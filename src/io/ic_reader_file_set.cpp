#include "cosmosim/io/ic_reader.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <functional>
#include <fstream>
#include <limits>
#include <numeric>
#include <optional>
#include <set>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/units.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>
#endif
#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace cosmosim::io {
namespace {

constexpr std::size_t kParticleTypeCount = 6U;
constexpr std::uint32_t kInvalidIndex = std::numeric_limits<std::uint32_t>::max();

[[nodiscard]] bool nearlyEqual(double lhs, double rhs) {
  return std::abs(lhs - rhs) <= 1.0e-10 * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

[[nodiscard]] IcSpeciesPolicy mapConfiguredPolicy(
    core::InitialConditionSpeciesPolicy policy,
    std::size_t type_index) {
  switch (policy) {
    case core::InitialConditionSpeciesPolicy::kReject:
      return IcSpeciesPolicy::kReject;
    case core::InitialConditionSpeciesPolicy::kDarkMatter:
      return type_index == 2U
          ? IcSpeciesPolicy::kCollisionlessFamily2AsDarkMatter
          : IcSpeciesPolicy::kCollisionlessFamily3AsDarkMatter;
    case core::InitialConditionSpeciesPolicy::kStar:
      return IcSpeciesPolicy::kStar;
    case core::InitialConditionSpeciesPolicy::kBlackHole:
      return IcSpeciesPolicy::kBlackHole;
    case core::InitialConditionSpeciesPolicy::kTracer:
      return IcSpeciesPolicy::kTracer;
  }
  throw std::invalid_argument("unknown configured IC species policy");
}

[[nodiscard]] std::uint32_t speciesTag(IcSpeciesPolicy policy) {
  switch (policy) {
    case IcSpeciesPolicy::kGas:
      return static_cast<std::uint32_t>(core::ParticleSpecies::kGas);
    case IcSpeciesPolicy::kDarkMatter:
    case IcSpeciesPolicy::kCollisionlessFamily2AsDarkMatter:
    case IcSpeciesPolicy::kCollisionlessFamily3AsDarkMatter:
      return static_cast<std::uint32_t>(core::ParticleSpecies::kDarkMatter);
    case IcSpeciesPolicy::kStar:
      return static_cast<std::uint32_t>(core::ParticleSpecies::kStar);
    case IcSpeciesPolicy::kBlackHole:
      return static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole);
    case IcSpeciesPolicy::kTracer:
      return static_cast<std::uint32_t>(core::ParticleSpecies::kTracer);
    case IcSpeciesPolicy::kReject:
      break;
  }
  throw std::runtime_error("attempted to materialize a rejected IC family");
}

class Sha256 {
 public:
  Sha256() { reset(); }
  void update(const std::uint8_t* data, std::size_t size) {
    for (std::size_t i = 0; i < size; ++i) {
      m_block[m_block_size++] = data[i];
      m_bit_count += 8U;
      if (m_block_size == 64U) {
        transform();
        m_block_size = 0U;
      }
    }
  }
  [[nodiscard]] std::array<std::uint8_t, 32> finish() {
    const std::uint64_t original_bits = m_bit_count;
    m_block[m_block_size++] = 0x80U;
    if (m_block_size > 56U) {
      while (m_block_size < 64U) m_block[m_block_size++] = 0U;
      transform();
      m_block_size = 0U;
    }
    while (m_block_size < 56U) m_block[m_block_size++] = 0U;
    for (int shift = 56; shift >= 0; shift -= 8) {
      m_block[m_block_size++] = static_cast<std::uint8_t>((original_bits >> shift) & 0xffU);
    }
    transform();
    std::array<std::uint8_t, 32> digest{};
    for (std::size_t i = 0; i < 8U; ++i) {
      digest[i * 4U + 0U] = static_cast<std::uint8_t>(m_state[i] >> 24U);
      digest[i * 4U + 1U] = static_cast<std::uint8_t>(m_state[i] >> 16U);
      digest[i * 4U + 2U] = static_cast<std::uint8_t>(m_state[i] >> 8U);
      digest[i * 4U + 3U] = static_cast<std::uint8_t>(m_state[i]);
    }
    return digest;
  }
 private:
  static constexpr std::array<std::uint32_t, 64> kRound{
      0x428a2f98U,0x71374491U,0xb5c0fbcfU,0xe9b5dba5U,0x3956c25bU,0x59f111f1U,0x923f82a4U,0xab1c5ed5U,
      0xd807aa98U,0x12835b01U,0x243185beU,0x550c7dc3U,0x72be5d74U,0x80deb1feU,0x9bdc06a7U,0xc19bf174U,
      0xe49b69c1U,0xefbe4786U,0x0fc19dc6U,0x240ca1ccU,0x2de92c6fU,0x4a7484aaU,0x5cb0a9dcU,0x76f988daU,
      0x983e5152U,0xa831c66dU,0xb00327c8U,0xbf597fc7U,0xc6e00bf3U,0xd5a79147U,0x06ca6351U,0x14292967U,
      0x27b70a85U,0x2e1b2138U,0x4d2c6dfcU,0x53380d13U,0x650a7354U,0x766a0abbU,0x81c2c92eU,0x92722c85U,
      0xa2bfe8a1U,0xa81a664bU,0xc24b8b70U,0xc76c51a3U,0xd192e819U,0xd6990624U,0xf40e3585U,0x106aa070U,
      0x19a4c116U,0x1e376c08U,0x2748774cU,0x34b0bcb5U,0x391c0cb3U,0x4ed8aa4aU,0x5b9cca4fU,0x682e6ff3U,
      0x748f82eeU,0x78a5636fU,0x84c87814U,0x8cc70208U,0x90befffaU,0xa4506cebU,0xbef9a3f7U,0xc67178f2U};
  static std::uint32_t rotate(std::uint32_t value, unsigned bits) {
    return (value >> bits) | (value << (32U - bits));
  }
  void reset() {
    m_state = {0x6a09e667U,0xbb67ae85U,0x3c6ef372U,0xa54ff53aU,0x510e527fU,0x9b05688cU,0x1f83d9abU,0x5be0cd19U};
    m_block_size = 0U; m_bit_count = 0U;
  }
  void transform() {
    std::array<std::uint32_t, 64> words{};
    for (std::size_t i = 0; i < 16U; ++i) {
      words[i] = (static_cast<std::uint32_t>(m_block[i*4U]) << 24U) |
          (static_cast<std::uint32_t>(m_block[i*4U+1U]) << 16U) |
          (static_cast<std::uint32_t>(m_block[i*4U+2U]) << 8U) |
          static_cast<std::uint32_t>(m_block[i*4U+3U]);
    }
    for (std::size_t i = 16U; i < 64U; ++i) {
      const std::uint32_t s0 = rotate(words[i-15U],7U) ^ rotate(words[i-15U],18U) ^ (words[i-15U] >> 3U);
      const std::uint32_t s1 = rotate(words[i-2U],17U) ^ rotate(words[i-2U],19U) ^ (words[i-2U] >> 10U);
      words[i] = words[i-16U] + s0 + words[i-7U] + s1;
    }
    std::uint32_t a=m_state[0],b=m_state[1],c=m_state[2],d=m_state[3],e=m_state[4],f=m_state[5],g=m_state[6],h=m_state[7];
    for (std::size_t i=0;i<64U;++i) {
      const std::uint32_t s1=rotate(e,6U)^rotate(e,11U)^rotate(e,25U);
      const std::uint32_t choice=(e&f)^((~e)&g);
      const std::uint32_t temp1=h+s1+choice+kRound[i]+words[i];
      const std::uint32_t s0=rotate(a,2U)^rotate(a,13U)^rotate(a,22U);
      const std::uint32_t majority=(a&b)^(a&c)^(b&c);
      const std::uint32_t temp2=s0+majority;
      h=g;g=f;f=e;e=d+temp1;d=c;c=b;b=a;a=temp1+temp2;
    }
    m_state[0]+=a;m_state[1]+=b;m_state[2]+=c;m_state[3]+=d;m_state[4]+=e;m_state[5]+=f;m_state[6]+=g;m_state[7]+=h;
  }
  std::array<std::uint32_t,8> m_state{};
  std::array<std::uint8_t,64> m_block{};
  std::size_t m_block_size=0U;
  std::uint64_t m_bit_count=0U;
};

[[nodiscard]] std::string sha256Hex(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("failed to open IC source for SHA-256: " + path.string());
  Sha256 hash; std::array<std::uint8_t,1U<<16U> buffer{};
  while (input) {
    input.read(reinterpret_cast<char*>(buffer.data()), static_cast<std::streamsize>(buffer.size()));
    const auto count=input.gcount(); if(count>0) hash.update(buffer.data(),static_cast<std::size_t>(count));
  }
  if(!input.eof()) throw std::runtime_error("failed while hashing IC source: " + path.string());
  static constexpr char kHex[]="0123456789abcdef"; const auto digest=hash.finish(); std::string out(64U,'0');
  for(std::size_t i=0;i<digest.size();++i){out[i*2U]=kHex[digest[i]>>4U];out[i*2U+1U]=kHex[digest[i]&0xfU];} return out;
}

[[nodiscard]] std::string sha256Hex(std::string_view value) {
  Sha256 hash;
  hash.update(
      reinterpret_cast<const std::uint8_t*>(value.data()), value.size());
  static constexpr char kHex[] = "0123456789abcdef";
  const auto digest = hash.finish();
  std::string out(64U, '0');
  for (std::size_t i = 0; i < digest.size(); ++i) {
    out[i * 2U] = kHex[digest[i] >> 4U];
    out[i * 2U + 1U] = kHex[digest[i] & 0xfU];
  }
  return out;
}

#if COSMOSIM_ENABLE_HDF5

class Hdf5Handle {
 public:
  explicit Hdf5Handle(hid_t handle=-1):m_handle(handle){}
  Hdf5Handle(const Hdf5Handle&)=delete; Hdf5Handle& operator=(const Hdf5Handle&)=delete;
  Hdf5Handle(Hdf5Handle&& other) noexcept:m_handle(other.m_handle){other.m_handle=-1;}
  ~Hdf5Handle(){close();}
  [[nodiscard]] hid_t get() const noexcept{return m_handle;}
  [[nodiscard]] bool valid() const noexcept{return m_handle>=0;}
 private:
  void close(){if(m_handle<0)return; switch(H5Iget_type(m_handle)){case H5I_FILE:H5Fclose(m_handle);break;case H5I_GROUP:H5Gclose(m_handle);break;case H5I_DATASET:H5Dclose(m_handle);break;case H5I_DATASPACE:H5Sclose(m_handle);break;case H5I_ATTR:H5Aclose(m_handle);break;case H5I_DATATYPE:H5Tclose(m_handle);break;default:break;}m_handle=-1;}
  hid_t m_handle=-1;
};

[[nodiscard]] bool pathExists(hid_t parent,std::string_view path){return H5Lexists(parent,std::string(path).c_str(),H5P_DEFAULT)>0;}
[[nodiscard]] bool attributeExists(hid_t parent,std::string_view name){return H5Aexists(parent,std::string(name).c_str())>0;}

void readAttributeU32x6(hid_t group,const char* name,std::array<std::uint32_t,6>& values,bool required=true){Hdf5Handle a(H5Aopen(group,name,H5P_DEFAULT));if(!a.valid()){if(required)throw std::runtime_error(std::string("missing Header/")+name);return;}if(H5Aread(a.get(),H5T_NATIVE_UINT32,values.data())<0)throw std::runtime_error(std::string("failed to read Header/")+name);}
void readAttributeF64x6(hid_t group,const char* name,std::array<double,6>& values){Hdf5Handle a(H5Aopen(group,name,H5P_DEFAULT));if(!a.valid()||H5Aread(a.get(),H5T_NATIVE_DOUBLE,values.data())<0)throw std::runtime_error(std::string("failed to read Header/")+name);}
void readAttributeF64(hid_t group,const char* name,double& value){Hdf5Handle a(H5Aopen(group,name,H5P_DEFAULT));if(!a.valid()||H5Aread(a.get(),H5T_NATIVE_DOUBLE,&value)<0)throw std::runtime_error(std::string("failed to read Header/")+name);}
void readAttributeU32(hid_t group,const char* name,std::uint32_t& value){Hdf5Handle a(H5Aopen(group,name,H5P_DEFAULT));if(!a.valid()||H5Aread(a.get(),H5T_NATIVE_UINT32,&value)<0)throw std::runtime_error(std::string("failed to read Header/")+name);}
[[nodiscard]] std::string readAttributeString(hid_t group,const char* name){Hdf5Handle a(H5Aopen(group,name,H5P_DEFAULT));if(!a.valid())throw std::runtime_error(std::string("missing Header/")+name);Hdf5Handle type(H5Aget_type(a.get()));if(!type.valid()||H5Tget_class(type.get())!=H5T_STRING)throw std::runtime_error(std::string("Header/")+name+" must be a string");if(H5Tis_variable_str(type.get())>0){char* raw=nullptr;if(H5Aread(a.get(),type.get(),&raw)<0)throw std::runtime_error(std::string("failed to read Header/")+name);std::string value=raw==nullptr?std::string{}:std::string(raw);if(raw!=nullptr)H5free_memory(raw);return value;}const std::size_t width=H5Tget_size(type.get());if(width==0U)throw std::runtime_error(std::string("Header/")+name+" has zero-width string type");std::vector<char> raw(width+1U,'\0');if(H5Aread(a.get(),type.get(),raw.data())<0)throw std::runtime_error(std::string("failed to read Header/")+name);const auto end=std::find(raw.begin(),raw.begin()+static_cast<std::ptrdiff_t>(width),'\0');return std::string(raw.begin(),end);}

[[nodiscard]] IcSchemaSummary readHeader(hid_t header) {
  IcSchemaSummary s; std::array<std::uint32_t,6> local{},low{}; readAttributeU32x6(header,"NumPart_ThisFile",local); readAttributeU32x6(header,"NumPart_Total",low); readAttributeU32x6(header,"NumPart_Total_HighWord",s.total_count_high_word,false); readAttributeF64x6(header,"MassTable",s.mass_table); readAttributeF64(header,"Time",s.scale_factor);readAttributeF64(header,"Redshift",s.redshift);readAttributeF64(header,"BoxSize",s.box_size);readAttributeF64(header,"Omega0",s.omega_matter);readAttributeF64(header,"OmegaLambda",s.omega_lambda);readAttributeF64(header,"HubbleParam",s.hubble_param);readAttributeU32(header,"NumFilesPerSnapshot",s.num_files_per_snapshot);
  for(std::size_t i=0;i<6;++i){s.count_by_type[i]=local[i];s.total_count_by_type[i]=static_cast<std::uint64_t>(low[i])|(static_cast<std::uint64_t>(s.total_count_high_word[i])<<32U);} return s;
}

[[nodiscard]] std::string headerAuditText(const IcSchemaSummary& s) {
  std::ostringstream out; out.precision(std::numeric_limits<double>::max_digits10); out<<"NumFilesPerSnapshot="<<s.num_files_per_snapshot<<";NumPart_ThisFile=";
  for(std::size_t i=0;i<6;++i)out<<(i?",":"")<<s.count_by_type[i];out<<";NumPart_Total=";for(std::size_t i=0;i<6;++i)out<<(i?",":"")<<s.total_count_by_type[i];out<<";NumPart_Total_HighWord=";for(std::size_t i=0;i<6;++i)out<<(i?",":"")<<s.total_count_high_word[i];out<<";MassTable=";for(std::size_t i=0;i<6;++i)out<<(i?",":"")<<s.mass_table[i];out<<";BoxSize="<<s.box_size<<";Time="<<s.scale_factor<<";Redshift="<<s.redshift<<";Omega0="<<s.omega_matter<<";OmegaLambda="<<s.omega_lambda<<";HubbleParam="<<s.hubble_param;return out.str();
}

struct TypeDescription { std::string name; IcScalarClass scalar_class=IcScalarClass::kFloatingPoint; std::uint8_t width=0; bool is_signed=false; IcByteOrder order=IcByteOrder::kNotApplicable; };
[[nodiscard]] TypeDescription describeType(hid_t type) {
  TypeDescription d; const H5T_class_t cls=H5Tget_class(type); const std::size_t width=H5Tget_size(type); if(width==0U||width>255U)throw std::runtime_error("unsupported HDF5 datatype width");d.width=static_cast<std::uint8_t>(width);
  if(cls==H5T_FLOAT){d.scalar_class=IcScalarClass::kFloatingPoint;d.name="float"+std::to_string(width*8U);d.is_signed=true;}else if(cls==H5T_INTEGER){d.scalar_class=IcScalarClass::kInteger;const H5T_sign_t sign=H5Tget_sign(type);d.is_signed=sign==H5T_SGN_2;d.name=std::string(d.is_signed?"int":"uint")+std::to_string(width*8U);}else throw std::runtime_error("IC datasets must use integer or floating-point scalar types");
  switch(H5Tget_order(type)){case H5T_ORDER_LE:d.order=IcByteOrder::kLittleEndian;break;case H5T_ORDER_BE:d.order=IcByteOrder::kBigEndian;break;case H5T_ORDER_NONE:d.order=IcByteOrder::kNotApplicable;break;default:d.order=IcByteOrder::kNative;break;} return d;
}

[[nodiscard]] std::vector<std::uint64_t> dataspaceDimensions(hid_t space) {
  const int rank=H5Sget_simple_extent_ndims(space);if(rank<0)throw std::runtime_error("failed to inspect HDF5 dataspace rank");if(rank==0)return {};std::vector<hsize_t> raw(static_cast<std::size_t>(rank));if(H5Sget_simple_extent_dims(space,raw.data(),nullptr)<0)throw std::runtime_error("failed to inspect HDF5 dimensions");std::vector<std::uint64_t> dims;dims.reserve(raw.size());for(hsize_t value:raw)dims.push_back(static_cast<std::uint64_t>(value));return dims;
}

[[nodiscard]] std::string selectAlias(hid_t group,const std::vector<std::string>& aliases,bool required,std::string_view canonical) {
  for(const auto& alias:aliases)if(pathExists(group,alias))return alias;
  if(required)throw std::runtime_error("required dataset missing for "+std::string(canonical));return {};
}

[[nodiscard]] std::vector<std::filesystem::path> discoverFiles(const std::filesystem::path& requested,std::uint32_t count) {
  if(count==0U)throw std::runtime_error("NumFilesPerSnapshot must be positive");if(count==1U)return {requested.lexically_normal()};
  const std::filesystem::path dir=requested.parent_path();const std::string filename=requested.filename().string();const std::size_t extension_pos=filename.rfind('.');if(extension_pos==std::string::npos)throw std::runtime_error("multifile IC path requires a .hdf5 or .h5 extension");const std::string extension=filename.substr(extension_pos);std::string prefix=filename.substr(0,extension_pos);const std::size_t index_dot=prefix.rfind('.');if(index_dot!=std::string::npos){const std::string maybe_index=prefix.substr(index_dot+1U);if(!maybe_index.empty()&&std::all_of(maybe_index.begin(),maybe_index.end(),[](char c){return c>='0'&&c<='9';}))prefix=prefix.substr(0,index_dot);}
  std::vector<std::filesystem::path> files;files.reserve(count);for(std::uint32_t i=0;i<count;++i){auto candidate=(dir/(prefix+"."+std::to_string(i)+extension)).lexically_normal();if(!std::filesystem::is_regular_file(candidate))throw std::runtime_error("missing multifile IC member: "+candidate.string());files.push_back(std::move(candidate));}return files;
}

struct Convention { core::UnitSystem source_units; IcCoordinateFrame frame=IcCoordinateFrame::kComoving; IcVelocityConvention velocity=IcVelocityConvention::kPhysicalPeculiar; };
[[nodiscard]] Convention conventionFor(IcDialect dialect,hid_t header) {
  if(dialect==IcDialect::kChuiCanonicalV1){const std::string schema_name=readAttributeString(header,"ChuiIcSchemaName");std::uint32_t schema_version=0U;readAttributeU32(header,"ChuiIcSchemaVersion",schema_version);const std::string coordinate_frame=readAttributeString(header,"ChuiCoordinateFrame");const std::string velocity_convention=readAttributeString(header,"ChuiVelocityConvention");const std::string manifest_digest=readAttributeString(header,"ConversionManifestSha256");if(schema_name!="chui_canonical_v1"||schema_version!=1U||coordinate_frame!="comoving"||velocity_convention!="physical_peculiar")throw std::runtime_error("canonical CHUÍ IC header has an unsupported schema or convention");if(manifest_digest.size()!=64U||!std::all_of(manifest_digest.begin(),manifest_digest.end(),[](char c){return (c>='0'&&c<='9')||(c>='a'&&c<='f');}))throw std::runtime_error("canonical CHUÍ IC ConversionManifestSha256 must be 64 lowercase hexadecimal characters");double length=0,mass=0,velocity=0;readAttributeF64(header,"ChuiLengthUnitToSI",length);readAttributeF64(header,"ChuiMassUnitToSI",mass);readAttributeF64(header,"ChuiVelocityUnitToSI",velocity);if(!std::isfinite(length)||!std::isfinite(mass)||!std::isfinite(velocity)||length<=0.0||mass<=0.0||velocity<=0.0)throw std::runtime_error("canonical CHUÍ IC unit scales must be finite and positive");core::UnitSystem units;units.length_si_per_code=length;units.mass_si_per_code=mass;units.velocity_si_per_code=velocity;return {.source_units=units,.frame=IcCoordinateFrame::kComoving,.velocity=IcVelocityConvention::kPhysicalPeculiar};}
  return {.source_units=core::makeUnitSystem("kpc","msun","km_s"),.frame=IcCoordinateFrame::kComoving,.velocity=IcVelocityConvention::kPhysicalPeculiar};
}

[[nodiscard]] IcFieldManifest inspectDataset(
    hid_t group,std::uint32_t file_index,std::string canonical_path,std::string selected_alias,
    const Convention& convention,const core::UnitSystem& target_units,IcFieldSemantics semantics,
    IcVelocityConvention velocity=IcVelocityConvention::kNotVelocity) {
  Hdf5Handle dataset(H5Dopen2(group,selected_alias.c_str(),H5P_DEFAULT));if(!dataset.valid())throw std::runtime_error("failed to open dataset "+canonical_path);Hdf5Handle type(H5Dget_type(dataset.get()));Hdf5Handle space(H5Dget_space(dataset.get()));if(!type.valid()||!space.valid())throw std::runtime_error("failed to inspect dataset "+canonical_path);const auto td=describeType(type.get());const auto dims=dataspaceDimensions(space.get());
  double base=1.0;std::string source_unit="dimensionless",target_unit="dimensionless";
  if(semantics==IcFieldSemantics::kCoordinate){base=convention.source_units.length_si_per_code;source_unit="source_length";target_unit="runtime_length";}else if(semantics==IcFieldSemantics::kVelocity){base=convention.source_units.velocity_si_per_code;source_unit="source_velocity";target_unit="runtime_velocity";}else if(semantics==IcFieldSemantics::kExtensive){base=convention.source_units.mass_si_per_code;source_unit="source_mass";target_unit="runtime_mass";}else if(semantics==IcFieldSemantics::kSpecific){base=convention.source_units.velocity_si_per_code*convention.source_units.velocity_si_per_code;source_unit="source_velocity_squared";target_unit="runtime_specific_energy";}else if(canonical_path.ends_with("/Density")){base=convention.source_units.mass_si_per_code/std::pow(convention.source_units.length_si_per_code,3.0);source_unit="source_mass/source_length^3";target_unit="runtime_density";}else if(canonical_path.ends_with("/BH_Mdot")){base=convention.source_units.mass_si_per_code/convention.source_units.timeSiPerCode();source_unit="source_mass/source_time";target_unit="runtime_mass/runtime_time";}
  return {.source_file_index=file_index,.dataset_path=std::move(canonical_path),.selected_alias=std::move(selected_alias),.scalar_type=td.name,.scalar_class=td.scalar_class,.byte_width=td.width,.is_signed=td.is_signed,.byte_order=td.order,.rank=static_cast<std::uint8_t>(dims.size()),.dimensions=dims,.record_count=dims.empty()?1U:dims.front(),.base_unit_to_si=base,.hubble_exponent=0.0,.scale_factor_exponent=0.0,.coordinate_frame=convention.frame,.velocity_convention=velocity,.semantics=semantics,.disposition=IcFieldDisposition::kConverted,.source_unit=std::move(source_unit),.target_unit=std::move(target_unit),.conversion_equation="target = stored * base_unit_to_si * h^hubble_exponent * a^scale_factor_exponent / target_si_per_code"};
}

[[nodiscard]] IcFieldManifest inspectHeaderAttribute(
    hid_t header,std::uint32_t file_index,std::string name,const Convention& convention,
    IcFieldSemantics semantics) {
  Hdf5Handle attr(H5Aopen(header,name.c_str(),H5P_DEFAULT));if(!attr.valid())throw std::runtime_error("failed to open Header/"+name);Hdf5Handle type(H5Aget_type(attr.get()));Hdf5Handle space(H5Aget_space(attr.get()));const auto td=describeType(type.get());auto dims=dataspaceDimensions(space.get());double base=1.0;std::string source_unit="dimensionless",target_unit="dimensionless";if(semantics==IcFieldSemantics::kCoordinate){base=convention.source_units.length_si_per_code;source_unit="source_length";target_unit="runtime_length";}else if(semantics==IcFieldSemantics::kExtensive){base=convention.source_units.mass_si_per_code;source_unit="source_mass";target_unit="runtime_mass";}
  return {.source_file_index=file_index,.dataset_path="/Header/"+name,.selected_alias=name,.scalar_type=td.name,.scalar_class=td.scalar_class,.byte_width=td.width,.is_signed=td.is_signed,.byte_order=td.order,.rank=static_cast<std::uint8_t>(dims.size()),.dimensions=dims,.record_count=dims.empty()?1U:dims.front(),.base_unit_to_si=base,.hubble_exponent=0.0,.scale_factor_exponent=0.0,.coordinate_frame=convention.frame,.velocity_convention=IcVelocityConvention::kNotVelocity,.semantics=semantics,.disposition=IcFieldDisposition::kConverted,.source_unit=std::move(source_unit),.target_unit=std::move(target_unit),.conversion_equation="target = stored * base_unit_to_si * h^hubble_exponent * a^scale_factor_exponent / target_si_per_code"};
}


void validateCrossFileSchema(const IcManifest& manifest) {
  if (manifest.num_files_per_snapshot <= 1U) {
    return;
  }
  const auto comparable_dimensions = [](const IcFieldManifest& lhs,
                                        const IcFieldManifest& rhs) {
    if (lhs.rank != rhs.rank || lhs.dimensions.size() != rhs.dimensions.size()) {
      return false;
    }
    if (lhs.dataset_path.starts_with("/Header/")) {
      return lhs.dimensions == rhs.dimensions;
    }
    return std::equal(
        lhs.dimensions.begin() + (lhs.dimensions.empty() ? 0 : 1),
        lhs.dimensions.end(),
        rhs.dimensions.begin() + (rhs.dimensions.empty() ? 0 : 1));
  };
  for (const IcFieldManifest& baseline : manifest.fields) {
    if (baseline.source_file_index != 0U) {
      continue;
    }
    for (std::uint32_t file_index = 1U;
         file_index < manifest.num_files_per_snapshot;
         ++file_index) {
      const auto candidate = std::find_if(
          manifest.fields.begin(), manifest.fields.end(),
          [&](const IcFieldManifest& field) {
            return field.source_file_index == file_index &&
                   field.dataset_path == baseline.dataset_path;
          });
      if (candidate == manifest.fields.end()) {
        throw std::runtime_error(
            "inconsistent source schema across IC files: missing " +
            baseline.dataset_path + " in file index " +
            std::to_string(file_index));
      }
      if (candidate->selected_alias != baseline.selected_alias ||
          candidate->scalar_type != baseline.scalar_type ||
          candidate->scalar_class != baseline.scalar_class ||
          candidate->byte_width != baseline.byte_width ||
          candidate->is_signed != baseline.is_signed ||
          candidate->byte_order != baseline.byte_order ||
          !comparable_dimensions(baseline, *candidate)) {
        throw std::runtime_error(
            "inconsistent source schema across IC files for " +
            baseline.dataset_path);
      }
    }
  }
  for (const IcFieldManifest& field : manifest.fields) {
    if (field.source_file_index == 0U) {
      continue;
    }
    const auto baseline = std::find_if(
        manifest.fields.begin(), manifest.fields.end(),
        [&](const IcFieldManifest& candidate) {
          return candidate.source_file_index == 0U &&
                 candidate.dataset_path == field.dataset_path;
        });
    if (baseline == manifest.fields.end()) {
      throw std::runtime_error(
          "inconsistent source schema across IC files: unexpected " +
          field.dataset_path + " in file index " +
          std::to_string(field.source_file_index));
    }
  }
}

struct Inspection { IcManifest manifest; std::vector<IcSchemaSummary> schemas; };

[[nodiscard]] Inspection inspectFileSet(
    const std::filesystem::path& requested,const core::SimulationConfig& config,
    const IcImportOptions& options) {
  Hdf5Handle first_file(H5Fopen(requested.string().c_str(),H5F_ACC_RDONLY,H5P_DEFAULT));if(!first_file.valid())throw std::runtime_error("failed to open IC file: "+requested.string());Hdf5Handle first_header(H5Gopen2(first_file.get(),"/Header",H5P_DEFAULT));if(!first_header.valid())throw std::runtime_error("IC file missing /Header group");const IcSchemaSummary first_schema=readHeader(first_header.get());
  std::vector<std::filesystem::path> files;
  if(options.manifest!=nullptr&&!options.manifest->source_files.empty())files=options.manifest->source_files;else files=discoverFiles(requested,first_schema.num_files_per_snapshot);
  if(files.size()!=first_schema.num_files_per_snapshot)throw std::runtime_error("manifest/discovery file count disagrees with NumFilesPerSnapshot");
  Inspection inspection;IcManifest& m=inspection.manifest;m.dialect=config.mode.ic_convention==core::InitialConditionConvention::kChuiCanonicalV1?IcDialect::kChuiCanonicalV1:IcDialect::kGadgetArepoBridgeV1;m.dialect_version="1";m.converter_version="chui_runtime_inspector_v1";m.source_files=files;m.num_files_per_snapshot=static_cast<std::uint32_t>(files.size());m.species_policy[2]=mapConfiguredPolicy(config.mode.ic_part_type2_policy,2U);m.species_policy[3]=mapConfiguredPolicy(config.mode.ic_part_type3_policy,3U);
  if(options.manifest!=nullptr){m.dialect=options.manifest->dialect;m.dialect_version=options.manifest->dialect_version;m.species_policy=options.manifest->species_policy;}if(m.dialect_version!="1")throw std::runtime_error("unsupported IC dialect version: "+m.dialect_version);if(config.mode.ic_convention==core::InitialConditionConvention::kChuiCanonicalV1&&m.dialect!=IcDialect::kChuiCanonicalV1)throw std::runtime_error("chui_canonical_v1 configuration requires a canonical CHUÍ manifest/input");if(config.mode.ic_convention==core::InitialConditionConvention::kGadgetArepoBridgeV1&&m.dialect!=IcDialect::kGadgetArepoBridgeV1)throw std::runtime_error("gadget_arepo_bridge_v1 configuration requires a bridge manifest/input");
  const core::UnitSystem target_units=core::makeUnitSystem(config.units.length_unit,config.units.mass_unit,config.units.velocity_unit);
  std::array<std::uint64_t,6> summed{};
  for(std::size_t file_index=0;file_index<files.size();++file_index){const auto& path=files[file_index];if(!std::filesystem::is_regular_file(path))throw std::runtime_error("IC source is not a regular file: "+path.string());Hdf5Handle file(H5Fopen(path.string().c_str(),H5F_ACC_RDONLY,H5P_DEFAULT));if(!file.valid())throw std::runtime_error("failed to open IC member: "+path.string());Hdf5Handle header(H5Gopen2(file.get(),"/Header",H5P_DEFAULT));if(!header.valid())throw std::runtime_error("IC member missing /Header: "+path.string());const IcSchemaSummary schema=readHeader(header.get());inspection.schemas.push_back(schema);
    if(schema.num_files_per_snapshot!=files.size()||schema.total_count_by_type!=first_schema.total_count_by_type||schema.total_count_high_word!=first_schema.total_count_high_word||schema.mass_table!=first_schema.mass_table||!nearlyEqual(schema.box_size,first_schema.box_size)||!nearlyEqual(schema.scale_factor,first_schema.scale_factor)||!nearlyEqual(schema.redshift,first_schema.redshift)||!nearlyEqual(schema.omega_matter,first_schema.omega_matter)||!nearlyEqual(schema.omega_lambda,first_schema.omega_lambda)||!nearlyEqual(schema.hubble_param,first_schema.hubble_param))throw std::runtime_error("inconsistent cosmology, box, epoch, mass table, totals, or NumFilesPerSnapshot across IC files");
    for(std::size_t type=0;type<6;++type){if(summed[type]>std::numeric_limits<std::uint64_t>::max()-schema.count_by_type[type])throw std::overflow_error("IC file-set particle count overflow");summed[type]+=schema.count_by_type[type];}
    m.num_part_this_file.push_back(schema.count_by_type);m.source_file_sizes_bytes.push_back(std::filesystem::file_size(path));const std::string hash=sha256Hex(path);m.source_sha256.push_back(hash);m.source_provenance_ids.push_back("sha256:"+hash);m.original_header_attributes.push_back(headerAuditText(schema));const Convention convention=conventionFor(m.dialect,header.get());m.fields.push_back(inspectHeaderAttribute(header.get(),static_cast<std::uint32_t>(file_index),"MassTable",convention,IcFieldSemantics::kExtensive));m.fields.push_back(inspectHeaderAttribute(header.get(),static_cast<std::uint32_t>(file_index),"BoxSize",convention,IcFieldSemantics::kCoordinate));
    for(std::size_t type=0;type<6;++type){const std::uint64_t count=schema.count_by_type[type];if(count==0U)continue;if(m.species_policy[type]==IcSpeciesPolicy::kReject)throw std::invalid_argument("populated PartType"+std::to_string(type)+" has explicit reject policy");const std::string group_path="/PartType"+std::to_string(type);Hdf5Handle group(H5Gopen2(file.get(),group_path.c_str(),H5P_DEFAULT));if(!group.valid())throw std::runtime_error("missing "+group_path+" in "+path.string());
      const auto add=[&](std::string canonical,const std::vector<std::string>& aliases,bool required,IcFieldSemantics semantics,IcVelocityConvention velocity=IcVelocityConvention::kNotVelocity){const std::string selected=selectAlias(group.get(),aliases,required,group_path+"/"+canonical);if(selected.empty())return false;auto field=inspectDataset(group.get(),static_cast<std::uint32_t>(file_index),group_path+"/"+canonical,selected,convention,target_units,semantics,velocity);if(field.record_count!=count)throw std::runtime_error("dataset record count disagrees with NumPart_ThisFile for "+field.dataset_path);if((canonical=="Coordinates"||canonical=="Velocities")&&(field.rank!=2U||field.dimensions.size()!=2U||field.dimensions[1]!=3U))throw std::runtime_error(field.dataset_path+" must have dimensions [N,3]");if(canonical!="Coordinates"&&canonical!="Velocities"&&field.rank!=1U)throw std::runtime_error(field.dataset_path+" must have rank 1");m.fields.push_back(std::move(field));return true;};
      add("Coordinates",{"Coordinates","Position","POS"},true,IcFieldSemantics::kCoordinate);add("Velocities",{"Velocities","Velocity","VEL"},options.require_velocities,IcFieldSemantics::kVelocity,convention.velocity);add("ParticleIDs",{"ParticleIDs","ParticleID","ID"},options.require_particle_ids,IcFieldSemantics::kIdentifier);const bool has_masses=add("Masses",{"Masses","Mass"},false,IcFieldSemantics::kExtensive);if(!has_masses&&schema.mass_table[type]<=0.0)throw std::runtime_error(group_path+" requires Masses because MassTable is zero");
      if(type==0U){if(!add("InternalEnergy",{"InternalEnergy","U","Internal_Energy"},false,IcFieldSemantics::kSpecific))m.defaulted_fields.push_back(group_path+"/InternalEnergy=zero");if(!add("Density",{"Density","Rho"},false,IcFieldSemantics::kIntensive))m.defaulted_fields.push_back(group_path+"/Density=zero");if(add("Metallicity",{"Metallicity","GFM_Metallicity"},false,IcFieldSemantics::kIntensive))m.dropped_fields.push_back(group_path+"/Metallicity: no gas metallicity lane");if(add("SmoothingLength",{"SmoothingLength","Hsml","Smoothing_Length"},false,IcFieldSemantics::kCoordinate))m.dropped_fields.push_back(group_path+"/SmoothingLength: no gas smoothing-length lane");}
      if(m.species_policy[type]==IcSpeciesPolicy::kStar){if(!add("StellarFormationTime",{"GFM_StellarFormationTime","StellarFormationTime","BirthTime"},false,IcFieldSemantics::kIntensive))m.defaulted_fields.push_back(group_path+"/StellarFormationTime=Header/Time");if(!add("InitialMass",{"GFM_InitialMass","InitialMass","BirthMass"},false,IcFieldSemantics::kExtensive))m.defaulted_fields.push_back(group_path+"/InitialMass=particle_mass");if(!add("Metallicity",{"GFM_Metallicity","Metallicity"},false,IcFieldSemantics::kIntensive))m.defaulted_fields.push_back(group_path+"/Metallicity=zero");}
      if(m.species_policy[type]==IcSpeciesPolicy::kBlackHole){add("BH_Mass",{"BH_Mass"},true,IcFieldSemantics::kExtensive);if(!add("BH_Mdot",{"BH_Mdot"},false,IcFieldSemantics::kIntensive))m.defaulted_fields.push_back(group_path+"/BH_Mdot=zero");}
      if(m.species_policy[type]==IcSpeciesPolicy::kTracer){add("ParentParticleIDs",{"ParentParticleIDs","TracerParentIDs"},true,IcFieldSemantics::kIdentifier);add("InjectionStep",{"InjectionStep"},true,IcFieldSemantics::kIdentifier);add("HostCellIndex",{"HostCellIndex"},true,IcFieldSemantics::kIdentifier);add("MassFractionOfHost",{"MassFractionOfHost"},true,IcFieldSemantics::kIntensive);add("LastHostMass",{"LastHostMass"},true,IcFieldSemantics::kExtensive);add("CumulativeExchangedMass",{"CumulativeExchangedMass"},true,IcFieldSemantics::kExtensive);}
    }
  }
  validateCrossFileSchema(m);
  if(summed!=first_schema.total_count_by_type)throw std::runtime_error("summed NumPart_ThisFile does not equal reconstructed 64-bit NumPart_Total");m.num_part_total=first_schema.total_count_by_type;m.num_part_total_high_word=first_schema.total_count_high_word;m.mass_table=first_schema.mass_table;m.box_size=first_schema.box_size;m.scale_factor=first_schema.scale_factor;m.redshift=first_schema.redshift;m.omega_matter=first_schema.omega_matter;m.omega_lambda=first_schema.omega_lambda;m.hubble_param=first_schema.hubble_param;m.converted_fields.reserve(m.fields.size());for(const auto& field:m.fields){m.converted_fields.push_back(field.dataset_path);if(std::find(m.conversion_equations.begin(),m.conversion_equations.end(),field.conversion_equation)==m.conversion_equations.end())m.conversion_equations.push_back(field.conversion_equation);} 
  if(options.manifest!=nullptr&&!options.manifest->source_sha256.empty()){const IcManifest& supplied=*options.manifest;if(supplied.num_files_per_snapshot!=m.num_files_per_snapshot||supplied.source_files!=m.source_files||supplied.source_sha256!=m.source_sha256||supplied.source_file_sizes_bytes!=m.source_file_sizes_bytes||supplied.num_part_this_file!=m.num_part_this_file||supplied.num_part_total!=m.num_part_total||supplied.num_part_total_high_word!=m.num_part_total_high_word||supplied.mass_table!=m.mass_table||supplied.species_policy!=m.species_policy||!nearlyEqual(supplied.box_size,m.box_size)||!nearlyEqual(supplied.scale_factor,m.scale_factor)||!nearlyEqual(supplied.redshift,m.redshift)||!nearlyEqual(supplied.omega_matter,m.omega_matter)||!nearlyEqual(supplied.omega_lambda,m.omega_lambda)||!nearlyEqual(supplied.hubble_param,m.hubble_param))throw std::runtime_error("supplied IC manifest provenance, scientific header, policies, or counts do not match inspected source files");std::vector<IcFieldManifest> observed=m.fields;if(observed.size()!=supplied.fields.size())throw std::runtime_error("supplied IC manifest field count does not match actual HDF5 schema");for(const auto& expected:supplied.fields){const auto it=std::find_if(observed.begin(),observed.end(),[&](const IcFieldManifest& actual){return actual.source_file_index==expected.source_file_index&&actual.dataset_path==expected.dataset_path&&actual.selected_alias==expected.selected_alias;});if(it==observed.end()||it->scalar_type!=expected.scalar_type||it->scalar_class!=expected.scalar_class||it->byte_width!=expected.byte_width||it->is_signed!=expected.is_signed||it->byte_order!=expected.byte_order||it->rank!=expected.rank||it->dimensions!=expected.dimensions||it->record_count!=expected.record_count)throw std::runtime_error("supplied IC manifest schema does not match actual HDF5 field "+expected.dataset_path);}m=supplied;}
  validateIcManifest(m);return inspection;
}

[[nodiscard]] const IcFieldManifest* findField(const IcManifest& manifest,std::size_t file_index,std::string_view path){const auto it=std::find_if(manifest.fields.begin(),manifest.fields.end(),[&](const IcFieldManifest& field){return field.source_file_index==file_index&&field.dataset_path==path;});return it==manifest.fields.end()?nullptr:&*it;}
[[nodiscard]] const IcFieldManifest& requireField(const IcManifest& manifest,std::size_t file_index,std::string_view path){const auto* field=findField(manifest,file_index,path);if(field==nullptr)throw std::runtime_error("manifest lacks inspected field "+std::string(path)+" for file "+std::to_string(file_index));return *field;}

void readChunkDouble(hid_t group,const IcFieldManifest& field,std::size_t start,std::size_t count,std::size_t components,std::vector<double>& out){Hdf5Handle dataset(H5Dopen2(group,field.selected_alias.c_str(),H5P_DEFAULT));Hdf5Handle file_space(H5Dget_space(dataset.get()));if(!dataset.valid()||!file_space.valid())throw std::runtime_error("failed to open/read "+field.dataset_path);const int rank=components==1U?1:2;std::array<hsize_t,2> offset{static_cast<hsize_t>(start),0U},extent{static_cast<hsize_t>(count),static_cast<hsize_t>(components)};if(H5Sselect_hyperslab(file_space.get(),H5S_SELECT_SET,offset.data(),nullptr,extent.data(),nullptr)<0)throw std::runtime_error("failed hyperslab for "+field.dataset_path);Hdf5Handle mem(H5Screate_simple(rank,extent.data(),nullptr));out.resize(count*components);if(H5Dread(dataset.get(),H5T_NATIVE_DOUBLE,mem.get(),file_space.get(),H5P_DEFAULT,out.data())<0)throw std::runtime_error("failed chunk read for "+field.dataset_path);}
void readChunkU64(hid_t group,const IcFieldManifest& field,std::size_t start,std::size_t count,std::vector<std::uint64_t>& out){Hdf5Handle dataset(H5Dopen2(group,field.selected_alias.c_str(),H5P_DEFAULT));Hdf5Handle file_space(H5Dget_space(dataset.get()));hsize_t offset[1]{static_cast<hsize_t>(start)},extent[1]{static_cast<hsize_t>(count)};if(H5Sselect_hyperslab(file_space.get(),H5S_SELECT_SET,offset,nullptr,extent,nullptr)<0)throw std::runtime_error("failed hyperslab for "+field.dataset_path);Hdf5Handle mem(H5Screate_simple(1,extent,nullptr));out.resize(count);if(H5Dread(dataset.get(),H5T_NATIVE_UINT64,mem.get(),file_space.get(),H5P_DEFAULT,out.data())<0)throw std::runtime_error("failed ID chunk read for "+field.dataset_path);}

void convertValues(std::vector<double>& values,const IcFieldManifest& field,const IcManifest& manifest,double target_si_per_code){double factor=icStoredToSiMultiplier(field,manifest.hubble_param,manifest.scale_factor)/target_si_per_code;if(field.semantics==IcFieldSemantics::kVelocity)factor*=icVelocityConventionMultiplier(field.velocity_convention,manifest.scale_factor);for(double& value:values){value*=factor;if(!std::isfinite(value))throw std::runtime_error("non-finite converted value in "+field.dataset_path);}}

struct ParticleRecord {
  std::uint64_t id=0;std::uint32_t species=0;double x=0,y=0,z=0,vx=0,vy=0,vz=0,mass=0;
  double gas_density=0,gas_internal_energy=0;double star_formation=0,star_birth_mass=0,star_metallicity=0;double bh_mass=0,bh_mdot=0;
  std::uint64_t tracer_parent=0,tracer_injection=0;std::uint32_t tracer_host=kInvalidIndex;double tracer_fraction=0,tracer_last_host_mass=0,tracer_exchanged_mass=0;
};

[[nodiscard]] std::uint64_t precedingRecordCount(const IcManifest& manifest,std::size_t file_index,std::size_t type_index){std::uint64_t count=0;for(std::size_t file=0;file<manifest.num_part_this_file.size();++file){for(std::size_t type=0;type<6;++type){if(file==file_index&&type==type_index)return count;if(count>std::numeric_limits<std::uint64_t>::max()-manifest.num_part_this_file[file][type])throw std::overflow_error("generated IC ID offset overflow");count+=manifest.num_part_this_file[file][type];}}throw std::logic_error("invalid file/type offset");}

[[nodiscard]] std::vector<ParticleRecord> readRecordChunk(const Inspection& inspection,std::size_t file_index,std::size_t type_index,std::size_t start,std::size_t count,const core::SimulationConfig& config,const IcImportOptions& options,IcImportCounters& counters){const IcManifest& manifest=inspection.manifest;const auto& path=manifest.source_files[file_index];Hdf5Handle file(H5Fopen(path.string().c_str(),H5F_ACC_RDONLY,H5P_DEFAULT));Hdf5Handle group(H5Gopen2(file.get(),("/PartType"+std::to_string(type_index)).c_str(),H5P_DEFAULT));if(!file.valid()||!group.valid())throw std::runtime_error("failed to open particle chunk source");const std::string prefix="/PartType"+std::to_string(type_index)+"/";const core::UnitSystem target=core::makeUnitSystem(config.units.length_unit,config.units.mass_unit,config.units.velocity_unit);
  std::vector<double> pos,vel,mass,gas_u,gas_rho,star_time,star_birth,star_metal,bh_mass,bh_mdot,tracer_fraction,tracer_last,tracer_exchange;std::vector<std::uint64_t> ids,tracer_parent,tracer_injection,tracer_host64;
  const auto& pf=requireField(manifest,file_index,prefix+"Coordinates");readChunkDouble(group.get(),pf,start,count,3U,pos);convertValues(pos,pf,manifest,target.length_si_per_code);if(pf.coordinate_frame==IcCoordinateFrame::kPhysical)for(double& value:pos)value=core::physicalToComovingLength(value,manifest.scale_factor);
  if(const auto* f=findField(manifest,file_index,prefix+"Velocities")){readChunkDouble(group.get(),*f,start,count,3U,vel);convertValues(vel,*f,manifest,target.velocity_si_per_code);}else vel.assign(count*3U,0.0);
  if(const auto* f=findField(manifest,file_index,prefix+"ParticleIDs"))readChunkU64(group.get(),*f,start,count,ids);else{ids.resize(count);const auto base=precedingRecordCount(manifest,file_index,type_index)+start+1U;for(std::size_t i=0;i<count;++i)ids[i]=base+i;}
  if(const auto* f=findField(manifest,file_index,prefix+"Masses")){readChunkDouble(group.get(),*f,start,count,1U,mass);convertValues(mass,*f,manifest,target.mass_si_per_code);}else{mass.assign(count,manifest.mass_table[type_index]);const auto& mass_header=requireField(manifest,file_index,"/Header/MassTable");convertValues(mass,mass_header,manifest,target.mass_si_per_code);}
  const IcSpeciesPolicy policy=manifest.species_policy[type_index];
  if(policy==IcSpeciesPolicy::kGas){if(const auto* f=findField(manifest,file_index,prefix+"InternalEnergy")){readChunkDouble(group.get(),*f,start,count,1U,gas_u);convertValues(gas_u,*f,manifest,target.velocity_si_per_code*target.velocity_si_per_code);}else gas_u.assign(count,0.0);if(const auto* f=findField(manifest,file_index,prefix+"Density")){readChunkDouble(group.get(),*f,start,count,1U,gas_rho);convertValues(gas_rho,*f,manifest,target.mass_si_per_code/std::pow(target.length_si_per_code,3.0));}else gas_rho.assign(count,0.0);}
  if(policy==IcSpeciesPolicy::kStar){if(const auto* f=findField(manifest,file_index,prefix+"StellarFormationTime"))readChunkDouble(group.get(),*f,start,count,1U,star_time);else star_time.assign(count,manifest.scale_factor);if(const auto* f=findField(manifest,file_index,prefix+"InitialMass")){readChunkDouble(group.get(),*f,start,count,1U,star_birth);convertValues(star_birth,*f,manifest,target.mass_si_per_code);}else star_birth=mass;if(const auto* f=findField(manifest,file_index,prefix+"Metallicity"))readChunkDouble(group.get(),*f,start,count,1U,star_metal);else star_metal.assign(count,0.0);}
  if(policy==IcSpeciesPolicy::kBlackHole){const auto& f=requireField(manifest,file_index,prefix+"BH_Mass");readChunkDouble(group.get(),f,start,count,1U,bh_mass);convertValues(bh_mass,f,manifest,target.mass_si_per_code);if(const auto* mdot=findField(manifest,file_index,prefix+"BH_Mdot")){readChunkDouble(group.get(),*mdot,start,count,1U,bh_mdot);convertValues(bh_mdot,*mdot,manifest,target.mass_si_per_code/target.timeSiPerCode());}else bh_mdot.assign(count,0.0);}
  if(policy==IcSpeciesPolicy::kTracer){readChunkU64(group.get(),requireField(manifest,file_index,prefix+"ParentParticleIDs"),start,count,tracer_parent);readChunkU64(group.get(),requireField(manifest,file_index,prefix+"InjectionStep"),start,count,tracer_injection);readChunkU64(group.get(),requireField(manifest,file_index,prefix+"HostCellIndex"),start,count,tracer_host64);const auto& ff=requireField(manifest,file_index,prefix+"MassFractionOfHost");readChunkDouble(group.get(),ff,start,count,1U,tracer_fraction);const auto& lf=requireField(manifest,file_index,prefix+"LastHostMass");readChunkDouble(group.get(),lf,start,count,1U,tracer_last);convertValues(tracer_last,lf,manifest,target.mass_si_per_code);const auto& ef=requireField(manifest,file_index,prefix+"CumulativeExchangedMass");readChunkDouble(group.get(),ef,start,count,1U,tracer_exchange);convertValues(tracer_exchange,ef,manifest,target.mass_si_per_code);}
  std::vector<ParticleRecord> records(count);for(std::size_t i=0;i<count;++i){ParticleRecord& r=records[i];r.id=ids[i];r.species=speciesTag(policy);r.x=pos[i*3U];r.y=pos[i*3U+1U];r.z=pos[i*3U+2U];r.vx=vel[i*3U];r.vy=vel[i*3U+1U];r.vz=vel[i*3U+2U];r.mass=mass[i];if(r.id==0U||!(r.mass>0.0)||!std::isfinite(r.mass))throw std::runtime_error("IC particle IDs must be nonzero and masses finite/positive");if(policy==IcSpeciesPolicy::kGas){r.gas_density=gas_rho[i];r.gas_internal_energy=gas_u[i];}if(policy==IcSpeciesPolicy::kStar){r.star_formation=star_time[i];r.star_birth_mass=star_birth[i];r.star_metallicity=star_metal[i];}if(policy==IcSpeciesPolicy::kBlackHole){r.bh_mass=bh_mass[i];r.bh_mdot=bh_mdot[i];if(!(r.bh_mass>0.0)||r.bh_mdot<0.0)throw std::runtime_error("black-hole sidecar values must be physical");}if(policy==IcSpeciesPolicy::kTracer){r.tracer_parent=tracer_parent[i];r.tracer_injection=tracer_injection[i];if(tracer_host64[i]>std::numeric_limits<std::uint32_t>::max())throw std::runtime_error("tracer host cell index exceeds uint32 range");r.tracer_host=static_cast<std::uint32_t>(tracer_host64[i]);r.tracer_fraction=tracer_fraction[i];r.tracer_last_host_mass=tracer_last[i];r.tracer_exchanged_mass=tracer_exchange[i];}}
  std::uint64_t bytes=0;for(const auto& field:manifest.fields)if(field.source_file_index==file_index&&field.dataset_path.starts_with(prefix))bytes+=static_cast<std::uint64_t>(count)*field.byte_width*(field.rank==2U?field.dimensions[1]:1U);counters.bytes_read+=bytes;counters.records_read+=count;counters.records_converted+=count;counters.peak_staging_bytes=std::max(counters.peak_staging_bytes,static_cast<std::uint64_t>(records.size()*sizeof(ParticleRecord)+pos.size()*sizeof(double)+vel.size()*sizeof(double)));return records;
}

void appendRecords(core::SimulationState& state,const std::vector<ParticleRecord>& records,std::uint32_t owner_rank){const std::size_t old_particles=state.particles.size();std::size_t gas_count=0,star_count=0,bh_count=0,tracer_count=0;for(const auto& r:records){if(r.species==static_cast<std::uint32_t>(core::ParticleSpecies::kGas))++gas_count;else if(r.species==static_cast<std::uint32_t>(core::ParticleSpecies::kStar))++star_count;else if(r.species==static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole))++bh_count;else if(r.species==static_cast<std::uint32_t>(core::ParticleSpecies::kTracer))++tracer_count;}
  const std::size_t old_gas=state.cells.size(),old_stars=state.star_particles.size(),old_bh=state.black_holes.size(),old_tracers=state.tracers.size();state.resizeParticles(old_particles+records.size());state.resizeCells(old_gas+gas_count);state.star_particles.resize(old_stars+star_count);state.black_holes.resize(old_bh+bh_count);state.tracers.resize(old_tracers+tracer_count);std::size_t gas_row=old_gas,star_row=old_stars,bh_row=old_bh,tracer_row=old_tracers;
  for(std::size_t i=0;i<records.size();++i){const auto& r=records[i];const std::size_t p=old_particles+i;state.particles.position_x_comoving[p]=r.x;state.particles.position_y_comoving[p]=r.y;state.particles.position_z_comoving[p]=r.z;state.particles.velocity_x_peculiar[p]=r.vx;state.particles.velocity_y_peculiar[p]=r.vy;state.particles.velocity_z_peculiar[p]=r.vz;state.particles.mass_code[p]=r.mass;state.particles.time_bin[p]=0U;state.particle_sidecar.particle_id[p]=r.id;state.particle_sidecar.species_tag[p]=r.species;state.particle_sidecar.owning_rank[p]=owner_rank;state.particle_sidecar.sfc_key[p]=0U;state.particle_sidecar.particle_flags[p]=0U;state.particle_sidecar.last_drift_time_code[p]=0.0;state.particle_sidecar.last_drift_scale_factor[p]=1.0;
    if(r.species==static_cast<std::uint32_t>(core::ParticleSpecies::kGas)){state.cells.center_x_comoving[gas_row]=r.x;state.cells.center_y_comoving[gas_row]=r.y;state.cells.center_z_comoving[gas_row]=r.z;state.cells.mass_code[gas_row]=r.mass;state.cells.time_bin[gas_row]=0U;state.cells.patch_index[gas_row]=0U;state.gas_cells.gas_cell_id[gas_row]=r.id;state.gas_cells.parent_particle_id[gas_row]=r.id;state.gas_cells.velocity_x_peculiar[gas_row]=r.vx;state.gas_cells.velocity_y_peculiar[gas_row]=r.vy;state.gas_cells.velocity_z_peculiar[gas_row]=r.vz;state.gas_cells.density_code[gas_row]=r.gas_density;state.gas_cells.internal_energy_code[gas_row]=r.gas_internal_energy;state.gas_cells.pressure_code[gas_row]=0.0;state.gas_cells.temperature_code[gas_row]=0.0;state.gas_cells.sound_speed_code[gas_row]=0.0;++gas_row;}
    else if(r.species==static_cast<std::uint32_t>(core::ParticleSpecies::kStar)){state.star_particles.particle_index[star_row]=static_cast<std::uint32_t>(p);state.star_particles.formation_scale_factor[star_row]=r.star_formation;state.star_particles.birth_mass_code[star_row]=r.star_birth_mass;state.star_particles.metallicity_mass_fraction[star_row]=r.star_metallicity;state.star_particles.stellar_age_years_last[star_row]=0.0;state.star_particles.stellar_returned_mass_cumulative_code[star_row]=0.0;state.star_particles.stellar_returned_metals_cumulative_code[star_row]=0.0;state.star_particles.stellar_feedback_energy_cumulative_erg[star_row]=0.0;for(std::size_t c=0;c<3U;++c){state.star_particles.stellar_returned_mass_channel_cumulative_code[c][star_row]=0.0;state.star_particles.stellar_returned_metals_channel_cumulative_code[c][star_row]=0.0;state.star_particles.stellar_feedback_energy_channel_cumulative_erg[c][star_row]=0.0;}++star_row;}
    else if(r.species==static_cast<std::uint32_t>(core::ParticleSpecies::kBlackHole)){state.black_holes.particle_index[bh_row]=static_cast<std::uint32_t>(p);state.black_holes.host_cell_index[bh_row]=kInvalidIndex;state.black_holes.subgrid_mass_code[bh_row]=r.bh_mass;state.black_holes.accretion_rate_code[bh_row]=r.bh_mdot;state.black_holes.feedback_energy_code[bh_row]=0.0;state.black_holes.eddington_ratio[bh_row]=0.0;state.black_holes.cumulative_accreted_mass_code[bh_row]=0.0;state.black_holes.cumulative_feedback_energy_code[bh_row]=0.0;state.black_holes.duty_cycle_active_time_code[bh_row]=0.0;state.black_holes.duty_cycle_total_time_code[bh_row]=0.0;++bh_row;}
    else if(r.species==static_cast<std::uint32_t>(core::ParticleSpecies::kTracer)){state.tracers.particle_index[tracer_row]=static_cast<std::uint32_t>(p);state.tracers.parent_particle_id[tracer_row]=r.tracer_parent;state.tracers.injection_step[tracer_row]=r.tracer_injection;state.tracers.host_cell_index[tracer_row]=r.tracer_host;state.tracers.mass_fraction_of_host[tracer_row]=r.tracer_fraction;state.tracers.last_host_mass_code[tracer_row]=r.tracer_last_host_mass;state.tracers.cumulative_exchanged_mass_code[tracer_row]=r.tracer_exchanged_mass;++tracer_row;}
  }
}

void finalizeImportedState(core::SimulationState& state,const IcManifest& manifest,const core::SimulationConfig& config){state.metadata.run_name=config.output.run_name;state.metadata.scale_factor=manifest.scale_factor;state.species.count_by_species={};for(const std::uint32_t species_tag:state.particle_sidecar.species_tag){if(species_tag>=state.species.count_by_species.size())throw std::runtime_error("IC import produced an out-of-range species tag");++state.species.count_by_species[species_tag];}state.rebuildSpeciesIndex();if(state.cells.size()>0U)state.refreshGasCellIdentityFromParticleOrder();if(!state.validateOwnershipInvariants())throw std::runtime_error("IC import produced invalid species/sidecar/ownership invariants");}

[[nodiscard]] double convertedBoxSizeCode(const IcManifest& manifest,const core::SimulationConfig& config){const core::UnitSystem target=core::makeUnitSystem(config.units.length_unit,config.units.mass_unit,config.units.velocity_unit);const auto& f=requireField(manifest,0U,"/Header/BoxSize");double box=manifest.box_size*icStoredToSiMultiplier(f,manifest.hubble_param,manifest.scale_factor)/target.length_si_per_code;if(f.coordinate_frame==IcCoordinateFrame::kPhysical)box=core::physicalToComovingLength(box,manifest.scale_factor);return box;}

void validateRuntimeCosmology(const IcManifest& manifest,const core::SimulationConfig& config){const double box_code=convertedBoxSizeCode(manifest,config);const core::UnitSystem target=core::makeUnitSystem(config.units.length_unit,config.units.mass_unit,config.units.velocity_unit);const core::UnitSystem mpc=core::makeUnitSystem("mpc","msun","km_s");const double box_mpc=box_code*target.length_si_per_code/mpc.length_si_per_code;if(!nearlyEqual(manifest.scale_factor,config.numerics.a_begin)||!nearlyEqual(manifest.omega_matter,config.cosmology.omega_matter)||!nearlyEqual(manifest.omega_lambda,config.cosmology.omega_lambda)||!nearlyEqual(manifest.hubble_param,config.cosmology.hubble_param)||!nearlyEqual(box_mpc,config.cosmology.box_size_x_mpc_comoving)||!nearlyEqual(box_mpc,config.cosmology.box_size_y_mpc_comoving)||!nearlyEqual(box_mpc,config.cosmology.box_size_z_mpc_comoving))throw std::runtime_error("IC manifest cosmology/BoxSize/start epoch does not match frozen runtime configuration");}

void validateSerialCountsAndIds(const core::SimulationState& state,const IcManifest& manifest){std::uint64_t expected=0;for(auto count:manifest.num_part_total){if(expected>std::numeric_limits<std::uint64_t>::max()-count)throw std::overflow_error("global particle count overflow");expected+=count;}if(state.particles.size()!=expected)throw std::runtime_error("IC import particle count mismatch");std::vector<std::uint64_t> ids(state.particle_sidecar.particle_id.begin(),state.particle_sidecar.particle_id.end());std::sort(ids.begin(),ids.end());if(std::adjacent_find(ids.begin(),ids.end())!=ids.end())throw std::runtime_error("IC import contains duplicate particle IDs");}

#endif  // COSMOSIM_ENABLE_HDF5

}  // namespace

IcReadResult readGadgetArepoHdf5Ic(const std::filesystem::path& ic_path,const core::SimulationConfig& config,const IcImportOptions& options) {
#if !COSMOSIM_ENABLE_HDF5
  static_cast<void>(ic_path);static_cast<void>(config);static_cast<void>(options);throw std::runtime_error("COSMOSIM_ENABLE_HDF5=OFF: GADGET/AREPO HDF5 IC reader unavailable in this build");
#else
  if(options.chunk_particle_count==0U)throw std::invalid_argument("chunk_particle_count must be positive");Inspection inspection=inspectFileSet(ic_path,config,options);validateRuntimeCosmology(inspection.manifest,config);IcReadResult result;result.report.manifest=inspection.manifest;result.report.schema=inspection.schemas.front();result.report.defaulted_fields=inspection.manifest.defaulted_fields;for(const auto& value:result.report.defaulted_fields){const auto eq=value.find('=');result.report.missing_optional_fields.push_back(value.substr(0,eq));}result.report.unsupported_fields=inspection.manifest.dropped_fields;
  for(std::size_t file=0;file<inspection.manifest.source_files.size();++file){result.report.counters.files_assigned+=1U;for(std::size_t type=0;type<6;++type){const std::size_t total=static_cast<std::size_t>(inspection.manifest.num_part_this_file[file][type]);for(std::size_t start=0;start<total;start+=options.chunk_particle_count){const std::size_t count=std::min(options.chunk_particle_count,total-start);auto records=readRecordChunk(inspection,file,type,start,count,config,options,result.report.counters);appendRecords(result.state,records,0U);++result.report.counters.chunks_assigned;result.report.counters.records_routed+=records.size();}}}
  validateSerialCountsAndIds(result.state,inspection.manifest);finalizeImportedState(result.state,inspection.manifest,config);result.report.counters.final_local_particle_count=result.state.particles.size();result.report.counters.final_local_gas_cell_count=result.state.cells.size();result.report.counters.final_local_star_count=result.state.star_particles.size();result.report.counters.final_local_black_hole_count=result.state.black_holes.size();result.report.counters.final_local_tracer_count=result.state.tracers.size();return result;
#endif
}


namespace {
#if COSMOSIM_ENABLE_HDF5 && COSMOSIM_ENABLE_MPI

constexpr std::size_t kWireRecordBytes = 168U;

void appendLe32(std::vector<std::uint8_t>& out, std::uint32_t value) {
  for (unsigned shift = 0U; shift < 32U; shift += 8U) {
    out.push_back(static_cast<std::uint8_t>((value >> shift) & 0xffU));
  }
}
void appendLe64(std::vector<std::uint8_t>& out, std::uint64_t value) {
  for (unsigned shift = 0U; shift < 64U; shift += 8U) {
    out.push_back(static_cast<std::uint8_t>((value >> shift) & 0xffU));
  }
}
void appendDouble(std::vector<std::uint8_t>& out, double value) {
  appendLe64(out, std::bit_cast<std::uint64_t>(value));
}
[[nodiscard]] std::uint32_t readLe32(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  if (offset + 4U > bytes.size()) throw std::runtime_error("truncated IC wire uint32");
  std::uint32_t value = 0U;
  for (unsigned shift = 0U; shift < 32U; shift += 8U) value |= static_cast<std::uint32_t>(bytes[offset++]) << shift;
  return value;
}
[[nodiscard]] std::uint64_t readLe64(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  if (offset + 8U > bytes.size()) throw std::runtime_error("truncated IC wire uint64");
  std::uint64_t value = 0U;
  for (unsigned shift = 0U; shift < 64U; shift += 8U) value |= static_cast<std::uint64_t>(bytes[offset++]) << shift;
  return value;
}
[[nodiscard]] double readDouble(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  return std::bit_cast<double>(readLe64(bytes, offset));
}
void serializeRecord(const ParticleRecord& record, std::vector<std::uint8_t>& out) {
  const std::size_t begin = out.size();
  appendLe64(out, record.id); appendLe32(out, record.species);
  appendDouble(out, record.x); appendDouble(out, record.y); appendDouble(out, record.z);
  appendDouble(out, record.vx); appendDouble(out, record.vy); appendDouble(out, record.vz); appendDouble(out, record.mass);
  appendDouble(out, record.gas_density); appendDouble(out, record.gas_internal_energy);
  appendDouble(out, record.star_formation); appendDouble(out, record.star_birth_mass); appendDouble(out, record.star_metallicity);
  appendDouble(out, record.bh_mass); appendDouble(out, record.bh_mdot);
  appendLe64(out, record.tracer_parent); appendLe64(out, record.tracer_injection); appendLe32(out, record.tracer_host);
  appendDouble(out, record.tracer_fraction); appendDouble(out, record.tracer_last_host_mass); appendDouble(out, record.tracer_exchanged_mass);
  if (out.size() - begin != kWireRecordBytes) throw std::logic_error("IC wire record byte contract drifted");
}
[[nodiscard]] ParticleRecord deserializeRecord(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  if (offset + kWireRecordBytes > bytes.size()) throw std::runtime_error("truncated IC wire record");
  ParticleRecord r; r.id=readLe64(bytes,offset);r.species=readLe32(bytes,offset);
  r.x=readDouble(bytes,offset);r.y=readDouble(bytes,offset);r.z=readDouble(bytes,offset);
  r.vx=readDouble(bytes,offset);r.vy=readDouble(bytes,offset);r.vz=readDouble(bytes,offset);r.mass=readDouble(bytes,offset);
  r.gas_density=readDouble(bytes,offset);r.gas_internal_energy=readDouble(bytes,offset);
  r.star_formation=readDouble(bytes,offset);r.star_birth_mass=readDouble(bytes,offset);r.star_metallicity=readDouble(bytes,offset);
  r.bh_mass=readDouble(bytes,offset);r.bh_mdot=readDouble(bytes,offset);
  r.tracer_parent=readLe64(bytes,offset);r.tracer_injection=readLe64(bytes,offset);r.tracer_host=readLe32(bytes,offset);
  r.tracer_fraction=readDouble(bytes,offset);r.tracer_last_host_mass=readDouble(bytes,offset);r.tracer_exchanged_mass=readDouble(bytes,offset);return r;
}

void populateSchemasFromManifest(Inspection& inspection) {
  inspection.schemas.clear();
  inspection.schemas.reserve(inspection.manifest.num_part_this_file.size());
  for (const auto& per_file_counts : inspection.manifest.num_part_this_file) {
    IcSchemaSummary schema;
    schema.count_by_type = per_file_counts;
    schema.total_count_by_type = inspection.manifest.num_part_total;
    schema.total_count_high_word = inspection.manifest.num_part_total_high_word;
    schema.mass_table = inspection.manifest.mass_table;
    schema.num_files_per_snapshot =
        inspection.manifest.num_files_per_snapshot;
    schema.box_size = inspection.manifest.box_size;
    schema.scale_factor = inspection.manifest.scale_factor;
    schema.redshift = inspection.manifest.redshift;
    schema.omega_matter = inspection.manifest.omega_matter;
    schema.omega_lambda = inspection.manifest.omega_lambda;
    schema.hubble_param = inspection.manifest.hubble_param;
    inspection.schemas.push_back(schema);
  }
}

[[nodiscard]] std::string collectiveFailureMessage(
    const parallel::MpiContext& mpi_context,
    const std::exception_ptr& local_failure,
    std::string_view phase) {
  const int local_failed = local_failure ? 1 : 0;
  int failed_count = 0;
  MPI_Allreduce(&local_failed, &failed_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (failed_count == 0) return {};
  int candidate = local_failure ? mpi_context.worldRank() : mpi_context.worldSize();
  int failure_rank = mpi_context.worldSize();
  MPI_Allreduce(&candidate, &failure_rank, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  std::string message;
  if (mpi_context.worldRank() == failure_rank) {
    try { std::rethrow_exception(local_failure); }
    catch (const std::exception& error) { message = error.what(); }
    catch (...) { message = "unknown non-standard exception"; }
  }
  std::uint64_t length = message.size();
  MPI_Bcast(&length, 1, MPI_UINT64_T, failure_rank, MPI_COMM_WORLD);
  message.resize(static_cast<std::size_t>(length));
  if (length > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) throw std::overflow_error("collective failure message exceeds MPI int limit");
  MPI_Bcast(message.data(), static_cast<int>(length), MPI_CHAR, failure_rank, MPI_COMM_WORLD);
  return std::string(phase) + " failed on rank " + std::to_string(failure_rank) + ": " + message;
}

[[nodiscard]] std::string broadcastRootString(
    const parallel::MpiContext& mpi_context,
    std::string root_value) {
  std::uint64_t length = mpi_context.isRoot() ? root_value.size() : 0U;
  MPI_Bcast(&length, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  if (length > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) throw std::overflow_error("IC manifest broadcast exceeds MPI int limit");
  root_value.resize(static_cast<std::size_t>(length));
  MPI_Bcast(root_value.data(), static_cast<int>(length), MPI_CHAR, 0, MPI_COMM_WORLD);
  return root_value;
}

[[nodiscard]] std::vector<std::uint8_t> alltoallBytes(
    const parallel::MpiContext& mpi_context,
    const std::vector<std::vector<std::uint8_t>>& per_rank,
    std::uint64_t& sent_bytes,
    std::uint64_t& received_bytes) {
  const int world_size=mpi_context.worldSize();std::vector<int> send_counts(static_cast<std::size_t>(world_size),0),recv_counts(static_cast<std::size_t>(world_size),0),send_displs(static_cast<std::size_t>(world_size),0),recv_displs(static_cast<std::size_t>(world_size),0);
  std::uint64_t send_total64=0U;for(int rank=0;rank<world_size;++rank){const auto size=per_rank[static_cast<std::size_t>(rank)].size();if(size>static_cast<std::size_t>(std::numeric_limits<int>::max())||send_total64>static_cast<std::uint64_t>(std::numeric_limits<int>::max())-size)throw std::overflow_error("bounded IC send payload exceeds MPI int count/displacement");send_displs[static_cast<std::size_t>(rank)]=static_cast<int>(send_total64);send_counts[static_cast<std::size_t>(rank)]=static_cast<int>(size);send_total64+=size;}
  MPI_Alltoall(send_counts.data(),1,MPI_INT,recv_counts.data(),1,MPI_INT,MPI_COMM_WORLD);std::uint64_t recv_total64=0U;for(int rank=0;rank<world_size;++rank){const int count=recv_counts[static_cast<std::size_t>(rank)];if(count<0||recv_total64>static_cast<std::uint64_t>(std::numeric_limits<int>::max())-static_cast<std::uint64_t>(count))throw std::overflow_error("bounded IC receive payload exceeds MPI int count/displacement");recv_displs[static_cast<std::size_t>(rank)]=static_cast<int>(recv_total64);recv_total64+=static_cast<std::uint64_t>(count);}
  std::vector<std::uint8_t> send(static_cast<std::size_t>(send_total64));for(int rank=0;rank<world_size;++rank)std::copy(per_rank[static_cast<std::size_t>(rank)].begin(),per_rank[static_cast<std::size_t>(rank)].end(),send.begin()+send_displs[static_cast<std::size_t>(rank)]);std::vector<std::uint8_t> recv(static_cast<std::size_t>(recv_total64));MPI_Alltoallv(send.data(),send_counts.data(),send_displs.data(),MPI_BYTE,recv.data(),recv_counts.data(),recv_displs.data(),MPI_BYTE,MPI_COMM_WORLD);sent_bytes+=send_total64;received_bytes+=recv_total64;return recv;
}

[[nodiscard]] int ownerForX(double x, double box_size, int world_size) {
  if (!std::isfinite(x) || !(box_size > 0.0) || x < 0.0 || x > box_size * (1.0 + 1.0e-12)) throw std::runtime_error("converted IC coordinate is outside the periodic box");
  if (x >= box_size) x = 0.0;
  const double fraction=x/box_size;const int owner=std::min(world_size-1,static_cast<int>(fraction*world_size));return std::max(0,owner);
}

[[nodiscard]] std::uint64_t exactDistributedIdAudit(
    const parallel::MpiContext& mpi_context,
    std::span<const std::uint64_t> local_ids,
    std::size_t batch_count,
    IcImportCounters& counters) {
  if(batch_count==0U)throw std::invalid_argument("IC ID audit batch count must be positive");const auto checked_sum=[](std::initializer_list<std::uint64_t> values){std::uint64_t total=0U;for(const std::uint64_t value:values){if(total>std::numeric_limits<std::uint64_t>::max()-value)throw std::overflow_error("IC staging-byte counter overflow");total+=value;}return total;};const std::uint64_t local_rounds=(local_ids.size()+batch_count-1U)/batch_count;std::uint64_t rounds=0U;MPI_Allreduce(&local_rounds,&rounds,1,MPI_UINT64_T,MPI_MAX,MPI_COMM_WORLD);std::vector<std::uint64_t> validator_ids;std::uint64_t peak_bytes=0U;
  for(std::uint64_t round=0;round<rounds;++round){const std::size_t begin=static_cast<std::size_t>(std::min<std::uint64_t>(round*batch_count,local_ids.size()));const std::size_t end=std::min(local_ids.size(),begin+batch_count);std::vector<std::vector<std::uint8_t>> buckets(static_cast<std::size_t>(mpi_context.worldSize()));for(std::size_t i=begin;i<end;++i){const std::uint64_t id=local_ids[i];const int validator=static_cast<int>((id^(id>>33U)^(id<<11U))%static_cast<std::uint64_t>(mpi_context.worldSize()));appendLe64(buckets[static_cast<std::size_t>(validator)],id);}std::uint64_t sent=0U,received_bytes=0U;const auto received=alltoallBytes(mpi_context,buckets,sent,received_bytes);counters.bytes_sent+=sent;counters.bytes_received+=received_bytes;const std::uint64_t retained=static_cast<std::uint64_t>(validator_ids.capacity())*sizeof(std::uint64_t);peak_bytes=std::max(peak_bytes,checked_sum({retained,sent,sent,received_bytes}));if(received.size()%8U!=0U)throw std::runtime_error("IC ID validation wire size is corrupt");std::size_t offset=0;while(offset<received.size())validator_ids.push_back(readLe64(received,offset));peak_bytes=std::max(peak_bytes,checked_sum({static_cast<std::uint64_t>(validator_ids.capacity())*sizeof(std::uint64_t),sent,received_bytes}));}
  std::sort(validator_ids.begin(),validator_ids.end());const int local_duplicate=std::adjacent_find(validator_ids.begin(),validator_ids.end())!=validator_ids.end()?1:0;int any_duplicate=0;MPI_Allreduce(&local_duplicate,&any_duplicate,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);if(any_duplicate!=0)throw std::runtime_error("distributed IC import contains duplicate particle IDs across files or ranks");return peak_bytes;
}

void validateDistributedTotals(
    const parallel::MpiContext& mpi_context,
    const core::SimulationState& state,
    const IcManifest& manifest,
    const std::array<double,5>& local_source_mass) {
  std::array<std::uint64_t,5> local_counts{};std::array<double,5> local_final_mass{};for(std::size_t i=0;i<state.particles.size();++i){const auto species=state.particle_sidecar.species_tag[i];if(species>=5U)throw std::runtime_error("invalid species tag after distributed IC routing");++local_counts[species];local_final_mass[species]+=state.particles.mass_code[i];}
  std::array<std::uint64_t,5> global_counts{};std::array<double,5> global_final_mass{},global_source_mass{};MPI_Allreduce(local_counts.data(),global_counts.data(),5,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);MPI_Allreduce(local_final_mass.data(),global_final_mass.data(),5,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);MPI_Allreduce(local_source_mass.data(),global_source_mass.data(),5,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  std::array<std::uint64_t,5> expected_counts{};for(std::size_t type=0;type<6;++type){if(manifest.num_part_total[type]>0U)expected_counts[speciesTag(manifest.species_policy[type])]+=manifest.num_part_total[type];}
  if(global_counts!=expected_counts)throw std::runtime_error("distributed IC global species counts do not match the manifest");for(std::size_t species=0;species<5U;++species){const double tolerance=1.0e-12*std::max({1.0,std::abs(global_source_mass[species]),std::abs(global_final_mass[species])});if(std::abs(global_source_mass[species]-global_final_mass[species])>tolerance)throw std::runtime_error("distributed IC global species mass total changed during routing");}
}

#endif  // COSMOSIM_ENABLE_HDF5 && COSMOSIM_ENABLE_MPI

}  // namespace

IcReadResult readDistributedGadgetArepoHdf5Ic(
    const std::filesystem::path& ic_path,
    const core::SimulationConfig& config,
    const parallel::MpiContext& mpi_context,
    const IcImportOptions& options) {
#if !COSMOSIM_ENABLE_HDF5
  static_cast<void>(ic_path);static_cast<void>(config);static_cast<void>(mpi_context);static_cast<void>(options);throw std::runtime_error("COSMOSIM_ENABLE_HDF5=OFF: distributed IC reader unavailable");
#elif !COSMOSIM_ENABLE_MPI
  if(mpi_context.worldSize()!=1)throw std::runtime_error("distributed IC reader requires an MPI-enabled build for world_size > 1");return readGadgetArepoHdf5Ic(ic_path,config,options);
#else
  if(!mpi_context.isEnabled()||mpi_context.worldSize()==1)return readGadgetArepoHdf5Ic(ic_path,config,options);if(options.chunk_particle_count==0U)throw std::invalid_argument("chunk_particle_count must be positive");
  Inspection root_inspection;std::exception_ptr inspection_failure;if(mpi_context.isRoot()){try{root_inspection=inspectFileSet(ic_path,config,options);validateRuntimeCosmology(root_inspection.manifest,config);}catch(...){inspection_failure=std::current_exception();}}
  const std::string failure=collectiveFailureMessage(mpi_context,inspection_failure,"IC file-set inspection");if(!failure.empty())throw std::runtime_error(failure);std::string manifest_json=mpi_context.isRoot()?serializeIcManifestJson(root_inspection.manifest):std::string{};manifest_json=broadcastRootString(mpi_context,std::move(manifest_json));Inspection inspection;std::exception_ptr manifest_failure;try{if(mpi_context.isRoot())inspection=std::move(root_inspection);else inspection.manifest=deserializeIcManifestJson(manifest_json);populateSchemasFromManifest(inspection);}catch(...){manifest_failure=std::current_exception();}const std::string manifest_error=collectiveFailureMessage(mpi_context,manifest_failure,"IC manifest broadcast/validation");if(!manifest_error.empty())throw std::runtime_error(manifest_error);
  const std::string local_manifest_json=serializeIcManifestJson(inspection.manifest);const std::string local_manifest_digest=sha256Hex(std::string_view(local_manifest_json));std::string root_manifest_digest=mpi_context.isRoot()?local_manifest_digest:std::string{};root_manifest_digest=broadcastRootString(mpi_context,std::move(root_manifest_digest));const int local_digest_mismatch=local_manifest_digest==root_manifest_digest?0:1;int any_digest_mismatch=0;MPI_Allreduce(&local_digest_mismatch,&any_digest_mismatch,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);if(any_digest_mismatch!=0)throw std::runtime_error("IC manifest SHA-256 digest is inconsistent across ranks");
  IcReadResult result;result.report.manifest=inspection.manifest;result.report.already_partitioned=true;result.report.schema.count_by_type=inspection.manifest.num_part_this_file.front();result.report.schema.total_count_by_type=inspection.manifest.num_part_total;result.report.schema.total_count_high_word=inspection.manifest.num_part_total_high_word;result.report.schema.mass_table=inspection.manifest.mass_table;result.report.schema.num_files_per_snapshot=inspection.manifest.num_files_per_snapshot;result.report.schema.box_size=inspection.manifest.box_size;result.report.schema.scale_factor=inspection.manifest.scale_factor;result.report.schema.redshift=inspection.manifest.redshift;result.report.schema.omega_matter=inspection.manifest.omega_matter;result.report.schema.omega_lambda=inspection.manifest.omega_lambda;result.report.schema.hubble_param=inspection.manifest.hubble_param;result.report.defaulted_fields=inspection.manifest.defaulted_fields;for(const auto& value:result.report.defaulted_fields){const auto eq=value.find('=');result.report.missing_optional_fields.push_back(value.substr(0,eq));}result.report.unsupported_fields=inspection.manifest.dropped_fields;
  const double box_size=convertedBoxSizeCode(inspection.manifest,config);const std::size_t distributed_chunk_particle_count=std::min(options.chunk_particle_count,config.mode.ic_staging_particle_count);if(distributed_chunk_particle_count==0U)throw std::invalid_argument("distributed IC staging particle count must be positive");std::uint64_t global_chunk_index=0U;std::set<std::size_t> assigned_files;std::array<double,5> local_source_mass{};
  for(std::size_t file=0;file<inspection.manifest.source_files.size();++file){for(std::size_t type=0;type<6;++type){const std::size_t total=static_cast<std::size_t>(inspection.manifest.num_part_this_file[file][type]);for(std::size_t start=0;start<total;start+=distributed_chunk_particle_count,++global_chunk_index){const std::size_t count=std::min(distributed_chunk_particle_count,total-start);const int reader_rank=static_cast<int>(global_chunk_index%static_cast<std::uint64_t>(mpi_context.worldSize()));std::vector<ParticleRecord> records;std::exception_ptr local_read_failure;if(mpi_context.worldRank()==reader_rank){try{records=readRecordChunk(inspection,file,type,start,count,config,options,result.report.counters);assigned_files.insert(file);++result.report.counters.chunks_assigned;for(const auto& r:records)local_source_mass[r.species]+=r.mass;}catch(...){local_read_failure=std::current_exception();}}
        const std::string read_failure=collectiveFailureMessage(mpi_context,local_read_failure,"IC chunk read/convert");if(!read_failure.empty())throw std::runtime_error(read_failure);std::vector<std::vector<std::uint8_t>> per_rank(static_cast<std::size_t>(mpi_context.worldSize()));if(mpi_context.worldRank()==reader_rank){for(const auto& record:records){const int owner=ownerForX(record.x,box_size,mpi_context.worldSize());serializeRecord(record,per_rank[static_cast<std::size_t>(owner)]);}result.report.counters.records_routed+=records.size();}
        std::uint64_t sent=0,received=0;const auto inbound_bytes=alltoallBytes(mpi_context,per_rank,sent,received);result.report.counters.bytes_sent+=sent;result.report.counters.bytes_received+=received;const std::uint64_t record_bytes=static_cast<std::uint64_t>(records.size())*sizeof(ParticleRecord);result.report.counters.peak_staging_bytes=std::max(result.report.counters.peak_staging_bytes,record_bytes+sent+sent+received);if(inbound_bytes.size()%kWireRecordBytes!=0U)throw std::runtime_error("distributed IC wire payload has invalid length");std::vector<ParticleRecord> inbound;inbound.reserve(inbound_bytes.size()/kWireRecordBytes);std::size_t offset=0;while(offset<inbound_bytes.size())inbound.push_back(deserializeRecord(inbound_bytes,offset));result.report.counters.peak_staging_bytes=std::max(result.report.counters.peak_staging_bytes,record_bytes+sent+received+static_cast<std::uint64_t>(inbound.capacity())*sizeof(ParticleRecord));appendRecords(result.state,inbound,static_cast<std::uint32_t>(mpi_context.worldRank()));
      }}}
  result.report.counters.files_assigned=assigned_files.size();finalizeImportedState(result.state,inspection.manifest,config);result.report.counters.peak_staging_bytes=std::max(result.report.counters.peak_staging_bytes,exactDistributedIdAudit(mpi_context,result.state.particle_sidecar.particle_id,config.mode.ic_staging_particle_count,result.report.counters));validateDistributedTotals(mpi_context,result.state,inspection.manifest,local_source_mass);result.report.counters.final_local_particle_count=result.state.particles.size();result.report.counters.final_local_gas_cell_count=result.state.cells.size();result.report.counters.final_local_star_count=result.state.star_particles.size();result.report.counters.final_local_black_hole_count=result.state.black_holes.size();result.report.counters.final_local_tracer_count=result.state.tracers.size();return result;
#endif
}

}  // namespace cosmosim::io
