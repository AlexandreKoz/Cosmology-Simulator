#include "cosmosim/io/ic_reader.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cstdlib>
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

[[nodiscard]] bool pathExists(hid_t parent, std::string_view path) {
  return H5Lexists(parent, std::string(path).c_str(), H5P_DEFAULT) > 0;
}
[[nodiscard]] bool attributeExists(hid_t parent, std::string_view name) {
  return H5Aexists(parent, std::string(name).c_str()) > 0;
}

struct TypeDescription {
  std::string name;
  IcScalarClass scalar_class = IcScalarClass::kFloatingPoint;
  std::uint8_t width = 0;
  bool is_signed = false;
  IcByteOrder order = IcByteOrder::kNotApplicable;
};

[[nodiscard]] TypeDescription describeType(hid_t type) {
  TypeDescription description;
  const H5T_class_t scalar_class = H5Tget_class(type);
  const std::size_t width = H5Tget_size(type);
  if (width == 0U || width > 255U) {
    throw std::runtime_error("unsupported HDF5 datatype width");
  }
  description.width = static_cast<std::uint8_t>(width);
  if (scalar_class == H5T_FLOAT) {
    description.scalar_class = IcScalarClass::kFloatingPoint;
    description.name = "float" + std::to_string(width * 8U);
    description.is_signed = true;
  } else if (scalar_class == H5T_INTEGER) {
    description.scalar_class = IcScalarClass::kInteger;
    const H5T_sign_t sign = H5Tget_sign(type);
    if (sign != H5T_SGN_NONE && sign != H5T_SGN_2) {
      throw std::runtime_error("unsupported HDF5 integer sign encoding");
    }
    description.is_signed = sign == H5T_SGN_2;
    description.name =
        std::string(description.is_signed ? "int" : "uint") +
        std::to_string(width * 8U);
  } else {
    throw std::runtime_error(
        "IC scalar data must use an integer or floating-point HDF5 type");
  }
  switch (H5Tget_order(type)) {
    case H5T_ORDER_LE:
      description.order = IcByteOrder::kLittleEndian;
      break;
    case H5T_ORDER_BE:
      description.order = IcByteOrder::kBigEndian;
      break;
    case H5T_ORDER_NONE:
      description.order = IcByteOrder::kNotApplicable;
      break;
    default:
      description.order = IcByteOrder::kNative;
      break;
  }
  return description;
}

[[nodiscard]] std::vector<std::uint64_t> attributeDimensions(hid_t attribute) {
  Hdf5Handle space(H5Aget_space(attribute));
  if (!space.valid()) {
    throw std::runtime_error("failed to inspect HDF5 attribute dataspace");
  }
  const int rank = H5Sget_simple_extent_ndims(space.get());
  if (rank < 0) {
    throw std::runtime_error("failed to inspect HDF5 attribute rank");
  }
  if (rank == 0) {
    return {};
  }
  std::vector<hsize_t> raw(static_cast<std::size_t>(rank));
  if (H5Sget_simple_extent_dims(space.get(), raw.data(), nullptr) < 0) {
    throw std::runtime_error("failed to inspect HDF5 attribute dimensions");
  }
  std::vector<std::uint64_t> dimensions;
  dimensions.reserve(raw.size());
  for (const hsize_t value : raw) {
    dimensions.push_back(static_cast<std::uint64_t>(value));
  }
  return dimensions;
}

[[nodiscard]] Hdf5Handle openValidatedAttribute(
    hid_t group,
    const char* name,
    bool required,
    IcScalarClass expected_class,
    std::span<const std::uint64_t> expected_dimensions,
    bool require_unsigned_integer = false) {
  Hdf5Handle attribute(H5Aopen(group, name, H5P_DEFAULT));
  if (!attribute.valid()) {
    if (required) {
      throw std::runtime_error(std::string("missing Header/") + name);
    }
    return attribute;
  }
  Hdf5Handle type(H5Aget_type(attribute.get()));
  if (!type.valid()) {
    throw std::runtime_error(std::string("failed to inspect Header/") + name);
  }
  const TypeDescription description = describeType(type.get());
  if (description.scalar_class != expected_class) {
    throw std::runtime_error(
        std::string("Header/") + name + " has an invalid datatype class");
  }
  if (require_unsigned_integer && description.is_signed) {
    throw std::runtime_error(
        std::string("Header/") + name + " must be unsigned integer data");
  }
  if (expected_class == IcScalarClass::kInteger && description.width > 4U) {
    throw std::runtime_error(
        std::string("Header/") + name +
        " must fit the uint32 header contract");
  }
  if (expected_class == IcScalarClass::kFloatingPoint &&
      description.width != 4U && description.width != 8U) {
    throw std::runtime_error(
        std::string("Header/") + name + " must use float32 or float64");
  }
  const std::vector<std::uint64_t> dimensions =
      attributeDimensions(attribute.get());
  if (!std::equal(
          dimensions.begin(), dimensions.end(), expected_dimensions.begin(),
          expected_dimensions.end())) {
    std::ostringstream message;
    message << "Header/" << name << " has shape [";
    for (std::size_t i = 0; i < dimensions.size(); ++i) {
      message << (i == 0 ? "" : ",") << dimensions[i];
    }
    message << "] but expected [";
    for (std::size_t i = 0; i < expected_dimensions.size(); ++i) {
      message << (i == 0 ? "" : ",") << expected_dimensions[i];
    }
    message << ']';
    throw std::runtime_error(message.str());
  }
  return attribute;
}

void readAttributeU32x6(
    hid_t group,
    const char* name,
    std::array<std::uint32_t, 6>& values,
    bool required = true) {
  static constexpr std::array<std::uint64_t, 1> kExpected{6U};
  Hdf5Handle attribute = openValidatedAttribute(
      group, name, required, IcScalarClass::kInteger, kExpected, true);
  if (!attribute.valid()) {
    return;
  }
  if (H5Aread(attribute.get(), H5T_NATIVE_UINT32, values.data()) < 0) {
    throw std::runtime_error(std::string("failed to read Header/") + name);
  }
}

void readAttributeF64x6(
    hid_t group,
    const char* name,
    std::array<double, 6>& values) {
  static constexpr std::array<std::uint64_t, 1> kExpected{6U};
  Hdf5Handle attribute = openValidatedAttribute(
      group, name, true, IcScalarClass::kFloatingPoint, kExpected);
  if (H5Aread(attribute.get(), H5T_NATIVE_DOUBLE, values.data()) < 0) {
    throw std::runtime_error(std::string("failed to read Header/") + name);
  }
}

void readAttributeF64(hid_t group, const char* name, double& value) {
  Hdf5Handle attribute = openValidatedAttribute(
      group, name, true, IcScalarClass::kFloatingPoint, {});
  if (H5Aread(attribute.get(), H5T_NATIVE_DOUBLE, &value) < 0) {
    throw std::runtime_error(std::string("failed to read Header/") + name);
  }
}

void readAttributeU32(hid_t group, const char* name, std::uint32_t& value) {
  Hdf5Handle attribute = openValidatedAttribute(
      group, name, true, IcScalarClass::kInteger, {}, true);
  if (H5Aread(attribute.get(), H5T_NATIVE_UINT32, &value) < 0) {
    throw std::runtime_error(std::string("failed to read Header/") + name);
  }
}

[[nodiscard]] std::string readAttributeString(
    hid_t group,
    const char* name) {
  Hdf5Handle attribute(H5Aopen(group, name, H5P_DEFAULT));
  if (!attribute.valid()) {
    throw std::runtime_error(std::string("missing Header/") + name);
  }
  if (!attributeDimensions(attribute.get()).empty()) {
    throw std::runtime_error(
        std::string("Header/") + name + " must be a scalar string");
  }
  Hdf5Handle type(H5Aget_type(attribute.get()));
  if (!type.valid() || H5Tget_class(type.get()) != H5T_STRING) {
    throw std::runtime_error(std::string("Header/") + name + " must be a string");
  }
  if (H5Tis_variable_str(type.get()) > 0) {
    char* raw = nullptr;
    if (H5Aread(attribute.get(), type.get(), &raw) < 0) {
      throw std::runtime_error(std::string("failed to read Header/") + name);
    }
    std::string value = raw == nullptr ? std::string{} : std::string(raw);
    if (raw != nullptr) {
      H5free_memory(raw);
    }
    return value;
  }
  const std::size_t width = H5Tget_size(type.get());
  if (width == 0U) {
    throw std::runtime_error(
        std::string("Header/") + name + " has zero-width string type");
  }
  std::vector<char> raw(width + 1U, '\0');
  if (H5Aread(attribute.get(), type.get(), raw.data()) < 0) {
    throw std::runtime_error(std::string("failed to read Header/") + name);
  }
  const auto end = std::find(
      raw.begin(), raw.begin() + static_cast<std::ptrdiff_t>(width), '\0');
  return std::string(raw.begin(), end);
}

[[nodiscard]] IcSchemaSummary readHeader(hid_t header) {
  IcSchemaSummary summary;
  std::array<std::uint32_t, 6> local{};
  std::array<std::uint32_t, 6> low{};
  readAttributeU32x6(header, "NumPart_ThisFile", local);
  readAttributeU32x6(header, "NumPart_Total", low);
  readAttributeU32x6(
      header, "NumPart_Total_HighWord", summary.total_count_high_word, false);
  readAttributeF64x6(header, "MassTable", summary.mass_table);
  readAttributeF64(header, "Time", summary.scale_factor);
  readAttributeF64(header, "Redshift", summary.redshift);
  readAttributeF64(header, "BoxSize", summary.box_size);
  readAttributeF64(header, "Omega0", summary.omega_matter);
  readAttributeF64(header, "OmegaLambda", summary.omega_lambda);
  readAttributeF64(header, "HubbleParam", summary.hubble_param);
  readAttributeU32(
      header, "NumFilesPerSnapshot", summary.num_files_per_snapshot);
  for (std::size_t i = 0; i < 6; ++i) {
    summary.count_by_type[i] = local[i];
    summary.total_count_by_type[i] =
        static_cast<std::uint64_t>(low[i]) |
        (static_cast<std::uint64_t>(summary.total_count_high_word[i]) << 32U);
  }
  return summary;
}

[[nodiscard]] std::string headerAuditText(const IcSchemaSummary& schema) {
  std::ostringstream output;
  output.precision(std::numeric_limits<double>::max_digits10);
  output << "NumFilesPerSnapshot=" << schema.num_files_per_snapshot
         << ";NumPart_ThisFile=";
  for (std::size_t i = 0; i < schema.count_by_type.size(); ++i) {
    output << (i == 0U ? "" : ",") << schema.count_by_type[i];
  }
  output << ";NumPart_Total=";
  for (std::size_t i = 0; i < schema.total_count_by_type.size(); ++i) {
    output << (i == 0U ? "" : ",") << schema.total_count_by_type[i];
  }
  output << ";NumPart_Total_HighWord=";
  for (std::size_t i = 0; i < schema.total_count_high_word.size(); ++i) {
    output << (i == 0U ? "" : ",") << schema.total_count_high_word[i];
  }
  output << ";MassTable=";
  for (std::size_t i = 0; i < schema.mass_table.size(); ++i) {
    output << (i == 0U ? "" : ",") << schema.mass_table[i];
  }
  output << ";BoxSize=" << schema.box_size
         << ";Time=" << schema.scale_factor
         << ";Redshift=" << schema.redshift
         << ";Omega0=" << schema.omega_matter
         << ";OmegaLambda=" << schema.omega_lambda
         << ";HubbleParam=" << schema.hubble_param;
  return output.str();
}

[[nodiscard]] std::vector<std::uint64_t> dataspaceDimensions(hid_t space) {
  const int rank=H5Sget_simple_extent_ndims(space);if(rank<0)throw std::runtime_error("failed to inspect HDF5 dataspace rank");if(rank==0)return {};std::vector<hsize_t> raw(static_cast<std::size_t>(rank));if(H5Sget_simple_extent_dims(space,raw.data(),nullptr)<0)throw std::runtime_error("failed to inspect HDF5 dimensions");std::vector<std::uint64_t> dims;dims.reserve(raw.size());for(hsize_t value:raw)dims.push_back(static_cast<std::uint64_t>(value));return dims;
}

herr_t collectLinkName(
    hid_t,
    const char* name,
    const H5L_info_t*,
    void* context) {
  auto* names = static_cast<std::vector<std::string>*>(context);
  names->emplace_back(name);
  return 0;
}

[[nodiscard]] std::vector<std::string> listLinkNames(hid_t group) {
  std::vector<std::string> names;
  hsize_t index = 0U;
  if (H5Literate(
          group, H5_INDEX_NAME, H5_ITER_NATIVE, &index, collectLinkName,
          &names) < 0) {
    throw std::runtime_error("failed to enumerate HDF5 group members");
  }
  return names;
}

[[nodiscard]] H5O_type_t objectType(hid_t group, std::string_view name) {
  H5O_info_t info{};
  if (H5Oget_info_by_name(
          group, std::string(name).c_str(), &info, H5O_INFO_BASIC, H5P_DEFAULT) < 0) {
    throw std::runtime_error(
        "failed to inspect HDF5 object " + std::string(name));
  }
  return info.type;
}

[[nodiscard]] std::string selectAlias(
    hid_t group,
    const std::vector<std::string>& aliases,
    bool required,
    std::string_view canonical) {
  std::vector<std::string> present;
  for (const auto& alias : aliases) {
    if (pathExists(group, alias)) {
      present.push_back(alias);
    }
  }
  if (present.size() > 1U) {
    std::ostringstream message;
    message << "ambiguous aliases for " << canonical << ": ";
    for (std::size_t i = 0; i < present.size(); ++i) {
      message << (i == 0U ? "" : ", ") << present[i];
    }
    throw std::runtime_error(message.str());
  }
  if (!present.empty()) {
    return present.front();
  }
  if (required) {
    throw std::runtime_error(
        "required dataset missing for " + std::string(canonical));
  }
  return {};
}

[[nodiscard]] std::vector<std::filesystem::path> discoverFiles(
    const std::filesystem::path& requested,
    std::uint32_t count) {
  if (count == 0U) {
    throw std::runtime_error("NumFilesPerSnapshot must be positive");
  }
  if (count == 1U) {
    return {requested.lexically_normal()};
  }

  const std::filesystem::path directory = requested.parent_path();
  const std::string filename = requested.filename().string();
  const std::size_t extension_position = filename.rfind('.');
  if (extension_position == std::string::npos) {
    throw std::runtime_error(
        "multifile IC path requires a .hdf5 or .h5 extension");
  }
  const std::string extension = filename.substr(extension_position);
  std::string prefix = filename.substr(0U, extension_position);
  const std::size_t index_separator = prefix.rfind('.');
  if (index_separator != std::string::npos) {
    const std::string maybe_index = prefix.substr(index_separator + 1U);
    if (!maybe_index.empty() &&
        std::all_of(
            maybe_index.begin(), maybe_index.end(),
            [](char value) { return value >= '0' && value <= '9'; })) {
      prefix = prefix.substr(0U, index_separator);
    }
  }

  std::vector<std::filesystem::path> files;
  files.reserve(count);
  for (std::uint32_t file_index = 0U; file_index < count; ++file_index) {
    auto candidate =
        (directory /
         (prefix + "." + std::to_string(file_index) + extension))
            .lexically_normal();
    if (!std::filesystem::is_regular_file(candidate)) {
      throw std::runtime_error(
          "missing multifile IC member: " + candidate.string());
    }
    files.push_back(std::move(candidate));
  }
  return files;
}

struct Convention {
  core::UnitSystem source_units;
  IcCoordinateFrame frame = IcCoordinateFrame::kComoving;
  IcVelocityConvention velocity = IcVelocityConvention::kPhysicalPeculiar;
  double length_hubble_exponent = 0.0;
  double length_scale_factor_exponent = 0.0;
  double mass_hubble_exponent = 0.0;
  double mass_scale_factor_exponent = 0.0;
  double velocity_hubble_exponent = 0.0;
  double velocity_scale_factor_exponent = 0.0;
};

[[nodiscard]] IcCoordinateFrame mapBridgeFrame(
    core::InitialConditionCoordinateFrame frame) {
  switch (frame) {
    case core::InitialConditionCoordinateFrame::kComoving:
      return IcCoordinateFrame::kComoving;
    case core::InitialConditionCoordinateFrame::kPhysical:
      return IcCoordinateFrame::kPhysical;
    case core::InitialConditionCoordinateFrame::kUnspecified:
      break;
  }
  throw std::invalid_argument("IC bridge coordinate frame is unspecified");
}

[[nodiscard]] IcVelocityConvention mapBridgeVelocityConvention(
    core::InitialConditionVelocityConvention convention) {
  switch (convention) {
    case core::InitialConditionVelocityConvention::kPhysicalPeculiar:
      return IcVelocityConvention::kPhysicalPeculiar;
    case core::InitialConditionVelocityConvention::kSqrtAScaledPeculiar:
      return IcVelocityConvention::kSqrtAScaledPeculiar;
    case core::InitialConditionVelocityConvention::kComovingCoordinateRate:
      return IcVelocityConvention::kComovingCoordinateRate;
    case core::InitialConditionVelocityConvention::kUnspecified:
      break;
  }
  throw std::invalid_argument("IC bridge velocity convention is unspecified");
}

[[nodiscard]] Convention conventionFor(
    IcDialect dialect,
    hid_t header,
    const core::SimulationConfig& config,
    bool supplied_manifest) {
  if (dialect == IcDialect::kChuiCanonicalV1) {
    const std::string schema_name =
        readAttributeString(header, "ChuiIcSchemaName");
    std::uint32_t schema_version = 0U;
    readAttributeU32(header, "ChuiIcSchemaVersion", schema_version);
    const std::string coordinate_frame =
        readAttributeString(header, "ChuiCoordinateFrame");
    const std::string velocity_convention =
        readAttributeString(header, "ChuiVelocityConvention");
    const std::string manifest_digest =
        readAttributeString(header, "ConversionManifestSha256");
    if (schema_name != "chui_canonical_v1" || schema_version != 1U ||
        coordinate_frame != "comoving" ||
        velocity_convention != "physical_peculiar") {
      throw std::runtime_error(
          "canonical CHUÍ IC header has an unsupported schema or convention");
    }
    if (manifest_digest.size() != 64U ||
        !std::all_of(
            manifest_digest.begin(), manifest_digest.end(), [](char c) {
              return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
            })) {
      throw std::runtime_error(
          "canonical CHUÍ IC ConversionManifestSha256 must be 64 lowercase hexadecimal characters");
    }
    double length = 0.0;
    double mass = 0.0;
    double velocity = 0.0;
    readAttributeF64(header, "ChuiLengthUnitToSI", length);
    readAttributeF64(header, "ChuiMassUnitToSI", mass);
    readAttributeF64(header, "ChuiVelocityUnitToSI", velocity);
    if (!std::isfinite(length) || !std::isfinite(mass) ||
        !std::isfinite(velocity) || length <= 0.0 || mass <= 0.0 ||
        velocity <= 0.0) {
      throw std::runtime_error(
          "canonical CHUÍ IC unit scales must be finite and positive");
    }
    core::UnitSystem units;
    units.length_si_per_code = length;
    units.mass_si_per_code = mass;
    units.velocity_si_per_code = velocity;
    return {
        .source_units = units,
        .frame = IcCoordinateFrame::kComoving,
        .velocity = IcVelocityConvention::kPhysicalPeculiar};
  }

  if (supplied_manifest) {
    // The supplied schema-v3 manifest is the authority for every field. These
    // neutral values are used only while re-inspecting HDF5 datatype/shape.
    core::UnitSystem units;
    units.length_si_per_code = 1.0;
    units.mass_si_per_code = 1.0;
    units.velocity_si_per_code = 1.0;
    return {.source_units = units};
  }

  core::UnitSystem units;
  units.length_si_per_code =
      config.mode.ic_bridge_source_length_unit_to_si;
  units.mass_si_per_code = config.mode.ic_bridge_source_mass_unit_to_si;
  units.velocity_si_per_code =
      config.mode.ic_bridge_source_velocity_unit_to_si;
  if (!(units.length_si_per_code > 0.0) ||
      !(units.mass_si_per_code > 0.0) ||
      !(units.velocity_si_per_code > 0.0) ||
      !std::isfinite(units.length_si_per_code) ||
      !std::isfinite(units.mass_si_per_code) ||
      !std::isfinite(units.velocity_si_per_code)) {
    throw std::invalid_argument(
        "gadget_arepo_bridge_v1 requires positive explicit source unit scales");
  }
  return {
      .source_units = units,
      .frame = mapBridgeFrame(config.mode.ic_bridge_coordinate_frame),
      .velocity = mapBridgeVelocityConvention(
          config.mode.ic_bridge_velocity_convention),
      .length_hubble_exponent =
          config.mode.ic_bridge_length_hubble_exponent,
      .length_scale_factor_exponent =
          config.mode.ic_bridge_length_scale_factor_exponent,
      .mass_hubble_exponent = config.mode.ic_bridge_mass_hubble_exponent,
      .mass_scale_factor_exponent =
          config.mode.ic_bridge_mass_scale_factor_exponent,
      .velocity_hubble_exponent =
          config.mode.ic_bridge_velocity_hubble_exponent,
      .velocity_scale_factor_exponent =
          config.mode.ic_bridge_velocity_scale_factor_exponent};
}

struct FieldConversionContract {
  double base_unit_to_si = 1.0;
  double hubble_exponent = 0.0;
  double scale_factor_exponent = 0.0;
  std::int8_t length_power = 0;
  std::int8_t mass_power = 0;
  std::int8_t time_power = 0;
  double frame_scale_factor_exponent = 0.0;
  std::uint8_t velocity_convention_power = 0;
  std::string source_unit = "dimensionless";
  std::string target_unit = "dimensionless";
};

[[nodiscard]] FieldConversionContract fieldConversionContract(
    std::string_view canonical_path,
    IcFieldSemantics semantics,
    const Convention& convention) {
  FieldConversionContract contract;
  const bool physical_frame =
      convention.frame == IcCoordinateFrame::kPhysical;
  if (canonical_path.ends_with("/Density")) {
    contract.base_unit_to_si =
        convention.source_units.mass_si_per_code /
        std::pow(convention.source_units.length_si_per_code, 3.0);
    contract.hubble_exponent = convention.mass_hubble_exponent -
        3.0 * convention.length_hubble_exponent;
    contract.scale_factor_exponent =
        convention.mass_scale_factor_exponent -
        3.0 * convention.length_scale_factor_exponent;
    contract.length_power = -3;
    contract.mass_power = 1;
    contract.frame_scale_factor_exponent = physical_frame ? 3.0 : 0.0;
    contract.source_unit = "source_mass/source_length^3";
    contract.target_unit = "runtime_mass/runtime_length^3";
    return contract;
  }
  if (canonical_path.ends_with("/BH_Mdot")) {
    contract.base_unit_to_si =
        convention.source_units.mass_si_per_code /
        convention.source_units.timeSiPerCode();
    contract.hubble_exponent = convention.mass_hubble_exponent +
        convention.velocity_hubble_exponent -
        convention.length_hubble_exponent;
    contract.scale_factor_exponent =
        convention.mass_scale_factor_exponent +
        convention.velocity_scale_factor_exponent -
        convention.length_scale_factor_exponent;
    contract.mass_power = 1;
    contract.time_power = -1;
    contract.frame_scale_factor_exponent = physical_frame ? 1.0 : 0.0;
    contract.velocity_convention_power = 1U;
    contract.source_unit = "source_mass/source_time";
    contract.target_unit = "runtime_mass/runtime_time";
    return contract;
  }
  switch (semantics) {
    case IcFieldSemantics::kCoordinate:
      contract.base_unit_to_si = convention.source_units.length_si_per_code;
      contract.hubble_exponent = convention.length_hubble_exponent;
      contract.scale_factor_exponent =
          convention.length_scale_factor_exponent;
      contract.length_power = 1;
      contract.frame_scale_factor_exponent = physical_frame ? -1.0 : 0.0;
      contract.source_unit = "source_length";
      contract.target_unit = "runtime_comoving_length";
      break;
    case IcFieldSemantics::kVelocity:
      contract.base_unit_to_si = convention.source_units.velocity_si_per_code;
      contract.hubble_exponent = convention.velocity_hubble_exponent;
      contract.scale_factor_exponent =
          convention.velocity_scale_factor_exponent;
      contract.length_power = 1;
      contract.time_power = -1;
      contract.velocity_convention_power = 1U;
      contract.source_unit = "source_velocity";
      contract.target_unit = "runtime_peculiar_velocity";
      break;
    case IcFieldSemantics::kExtensive:
      contract.base_unit_to_si = convention.source_units.mass_si_per_code;
      contract.hubble_exponent = convention.mass_hubble_exponent;
      contract.scale_factor_exponent = convention.mass_scale_factor_exponent;
      contract.mass_power = 1;
      contract.source_unit = "source_mass";
      contract.target_unit = "runtime_mass";
      break;
    case IcFieldSemantics::kSpecific:
      contract.base_unit_to_si =
          convention.source_units.velocity_si_per_code *
          convention.source_units.velocity_si_per_code;
      contract.hubble_exponent = 2.0 * convention.velocity_hubble_exponent;
      contract.scale_factor_exponent =
          2.0 * convention.velocity_scale_factor_exponent;
      contract.length_power = 2;
      contract.time_power = -2;
      contract.velocity_convention_power = 2U;
      contract.source_unit = "source_velocity^2";
      contract.target_unit = "runtime_specific_energy";
      break;
    case IcFieldSemantics::kIdentifier:
    case IcFieldSemantics::kIntensive:
      break;
  }
  return contract;
}

void validateDatasetSemanticType(
    const TypeDescription& description,
    IcFieldSemantics semantics,
    std::string_view canonical_path) {
  if (semantics == IcFieldSemantics::kIdentifier) {
    if (description.scalar_class != IcScalarClass::kInteger ||
        description.is_signed || description.width > sizeof(std::uint64_t)) {
      throw std::runtime_error(
          std::string(canonical_path) +
          " must use an unsigned integer datatype no wider than uint64");
    }
    return;
  }
  if (description.scalar_class != IcScalarClass::kFloatingPoint ||
      (description.width != 4U && description.width != 8U)) {
    throw std::runtime_error(
        std::string(canonical_path) +
        " must use float32 or float64 physical data");
  }
}

[[nodiscard]] IcFieldManifest inspectDataset(
    hid_t group,
    std::uint32_t file_index,
    std::string canonical_path,
    std::string selected_alias,
    const Convention& convention,
    IcFieldSemantics semantics,
    IcVelocityConvention velocity = IcVelocityConvention::kNotVelocity) {
  Hdf5Handle dataset(
      H5Dopen2(group, selected_alias.c_str(), H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error("failed to open dataset " + canonical_path);
  }
  Hdf5Handle type(H5Dget_type(dataset.get()));
  Hdf5Handle space(H5Dget_space(dataset.get()));
  if (!type.valid() || !space.valid()) {
    throw std::runtime_error("failed to inspect dataset " + canonical_path);
  }
  const TypeDescription type_description = describeType(type.get());
  validateDatasetSemanticType(type_description, semantics, canonical_path);
  const auto dimensions = dataspaceDimensions(space.get());
  const FieldConversionContract contract =
      fieldConversionContract(canonical_path, semantics, convention);
  return {
      .source_file_index = file_index,
      .dataset_path = std::move(canonical_path),
      .selected_alias = std::move(selected_alias),
      .scalar_type = type_description.name,
      .scalar_class = type_description.scalar_class,
      .byte_width = type_description.width,
      .is_signed = type_description.is_signed,
      .byte_order = type_description.order,
      .rank = static_cast<std::uint8_t>(dimensions.size()),
      .dimensions = dimensions,
      .record_count = dimensions.empty() ? 1U : dimensions.front(),
      .base_unit_to_si = contract.base_unit_to_si,
      .hubble_exponent = contract.hubble_exponent,
      .scale_factor_exponent = contract.scale_factor_exponent,
      .length_power = contract.length_power,
      .mass_power = contract.mass_power,
      .time_power = contract.time_power,
      .frame_scale_factor_exponent =
          contract.frame_scale_factor_exponent,
      .velocity_convention_power = contract.velocity_convention_power,
      .coordinate_frame = convention.frame,
      .velocity_convention = contract.velocity_convention_power > 0U
          ? convention.velocity
          : velocity,
      .semantics = semantics,
      .disposition = IcFieldDisposition::kConverted,
      .source_unit = contract.source_unit,
      .target_unit = contract.target_unit,
      .conversion_equation =
          "target = stored * base_unit_to_si * h^hubble_exponent * "
          "a^(scale_factor_exponent + frame_scale_factor_exponent) * "
          "velocity_convention_multiplier^velocity_convention_power / "
          "target_si_per_code(L^length_power M^mass_power T^time_power)"};
}

[[nodiscard]] IcFieldManifest inspectDroppedDataset(
    hid_t group,
    std::uint32_t file_index,
    std::string dataset_path,
    std::string selected_alias) {
  Hdf5Handle dataset(
      H5Dopen2(group, selected_alias.c_str(), H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error("failed to open dataset " + dataset_path);
  }
  Hdf5Handle type(H5Dget_type(dataset.get()));
  Hdf5Handle space(H5Dget_space(dataset.get()));
  if (!type.valid() || !space.valid()) {
    throw std::runtime_error("failed to inspect dataset " + dataset_path);
  }
  const TypeDescription description = describeType(type.get());
  if (description.scalar_class == IcScalarClass::kFloatingPoint &&
      description.width != 4U && description.width != 8U) {
    throw std::runtime_error(
        dataset_path + " uses an unsupported floating-point width");
  }
  if (description.scalar_class == IcScalarClass::kInteger &&
      description.width > sizeof(std::uint64_t)) {
    throw std::runtime_error(dataset_path + " uses an integer wider than uint64");
  }
  const std::vector<std::uint64_t> dimensions =
      dataspaceDimensions(space.get());
  return {
      .source_file_index = file_index,
      .dataset_path = std::move(dataset_path),
      .selected_alias = std::move(selected_alias),
      .scalar_type = description.name,
      .scalar_class = description.scalar_class,
      .byte_width = description.width,
      .is_signed = description.is_signed,
      .byte_order = description.order,
      .rank = static_cast<std::uint8_t>(dimensions.size()),
      .dimensions = dimensions,
      .record_count = dimensions.empty() ? 1U : dimensions.front(),
      .base_unit_to_si = 1.0,
      .hubble_exponent = 0.0,
      .scale_factor_exponent = 0.0,
      .length_power = 0,
      .mass_power = 0,
      .time_power = 0,
      .frame_scale_factor_exponent = 0.0,
      .velocity_convention_power = 0U,
      .coordinate_frame = IcCoordinateFrame::kComoving,
      .velocity_convention = IcVelocityConvention::kNotVelocity,
      .semantics = IcFieldSemantics::kIntensive,
      .disposition = IcFieldDisposition::kDropped,
      .source_unit = "unknown",
      .target_unit = "not_imported",
      .conversion_equation = "not converted: explicitly dropped"};
}

[[nodiscard]] IcFieldManifest inspectHeaderAttribute(
    hid_t header,
    std::uint32_t file_index,
    std::string name,
    const Convention& convention,
    IcFieldSemantics semantics) {
  Hdf5Handle attribute(H5Aopen(header, name.c_str(), H5P_DEFAULT));
  if (!attribute.valid()) {
    throw std::runtime_error("failed to open Header/" + name);
  }
  Hdf5Handle type(H5Aget_type(attribute.get()));
  Hdf5Handle space(H5Aget_space(attribute.get()));
  if (!type.valid() || !space.valid()) {
    throw std::runtime_error("failed to inspect Header/" + name);
  }
  const TypeDescription type_description = describeType(type.get());
  validateDatasetSemanticType(
      type_description, semantics, "/Header/" + name);
  const auto dimensions = dataspaceDimensions(space.get());
  const FieldConversionContract contract =
      fieldConversionContract("/Header/" + name, semantics, convention);
  return {
      .source_file_index = file_index,
      .dataset_path = "/Header/" + name,
      .selected_alias = name,
      .scalar_type = type_description.name,
      .scalar_class = type_description.scalar_class,
      .byte_width = type_description.width,
      .is_signed = type_description.is_signed,
      .byte_order = type_description.order,
      .rank = static_cast<std::uint8_t>(dimensions.size()),
      .dimensions = dimensions,
      .record_count = dimensions.empty() ? 1U : dimensions.front(),
      .base_unit_to_si = contract.base_unit_to_si,
      .hubble_exponent = contract.hubble_exponent,
      .scale_factor_exponent = contract.scale_factor_exponent,
      .length_power = contract.length_power,
      .mass_power = contract.mass_power,
      .time_power = contract.time_power,
      .frame_scale_factor_exponent =
          contract.frame_scale_factor_exponent,
      .velocity_convention_power = contract.velocity_convention_power,
      .coordinate_frame = convention.frame,
      .velocity_convention = IcVelocityConvention::kNotVelocity,
      .semantics = semantics,
      .disposition = IcFieldDisposition::kConverted,
      .source_unit = contract.source_unit,
      .target_unit = contract.target_unit,
      .conversion_equation =
          "target = stored * base_unit_to_si * h^hubble_exponent * "
          "a^(scale_factor_exponent + frame_scale_factor_exponent) / "
          "target_si_per_code(L^length_power M^mass_power T^time_power)"};
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

void checkedCounterAdd(
    std::uint64_t& destination,
    std::uint64_t value,
    std::string_view counter_name);

struct Inspection {
  IcManifest manifest;
  std::vector<IcSchemaSummary> schemas;
  IcImportCounters counters;
};

[[nodiscard]] bool isSupportedPartTypeName(
    std::string_view name,
    std::size_t& type_index) {
  static constexpr std::string_view kPrefix = "PartType";
  if (!name.starts_with(kPrefix)) {
    return false;
  }
  const std::string_view suffix = name.substr(kPrefix.size());
  if (suffix.empty() ||
      !std::all_of(suffix.begin(), suffix.end(), [](char value) {
        return value >= '0' && value <= '9';
      })) {
    return false;
  }
  std::uint64_t parsed = 0U;
  for (char value : suffix) {
    parsed = parsed * 10U + static_cast<std::uint64_t>(value - '0');
    if (parsed > std::numeric_limits<std::size_t>::max()) {
      throw std::overflow_error("PartType index overflow");
    }
  }
  type_index = static_cast<std::size_t>(parsed);
  return type_index < kParticleTypeCount;
}

[[nodiscard]] std::uint64_t logicalHeaderPayloadBytes(
    const IcSchemaSummary& schema) {
  static_cast<void>(schema);
  return 3U * 6U * sizeof(std::uint32_t) +
      6U * sizeof(double) + 6U * sizeof(double) + sizeof(std::uint32_t);
}

void recordRootSchemaDisposition(
    hid_t file,
    const IcSchemaSummary& schema,
    IcManifest& manifest) {
  for (const std::string& name : listLinkNames(file)) {
    if (name == "Header") {
      continue;
    }
    std::size_t type_index = 0U;
    if (isSupportedPartTypeName(name, type_index)) {
      if (schema.count_by_type[type_index] == 0U) {
        if (objectType(file, name) != H5O_TYPE_GROUP) {
          throw std::runtime_error("/" + name + " must be an HDF5 group");
        }
        Hdf5Handle group(H5Gopen2(file, name.c_str(), H5P_DEFAULT));
        if (!group.valid()) {
          throw std::runtime_error("failed to inspect /" + name);
        }
        if (!listLinkNames(group.get()).empty()) {
          throw std::runtime_error(
              "/" + name +
              " contains datasets although NumPart_ThisFile is zero");
        }
      }
      continue;
    }
    if (name.starts_with("PartType")) {
      throw std::runtime_error(
          "unsupported populated particle-family object /" + name);
    }
    const std::string disposition =
        "/" + name + ": auxiliary root object is not imported";
    manifest.preserved_auxiliary_fields.push_back(disposition);
    manifest.warnings.push_back(disposition);
  }
}

struct SourceFileInspection {
  std::filesystem::path path;
  IcSchemaSummary schema;
  std::uint64_t source_size_bytes = 0U;
  std::string source_sha256;
  std::string original_header_attributes;
  std::vector<IcFieldManifest> fields;
  std::vector<std::string> defaulted_fields;
  std::vector<std::string> dropped_fields;
  std::vector<std::string> preserved_auxiliary_fields;
  std::vector<std::string> warnings;
  IcImportCounters counters;
};

[[nodiscard]] SourceFileInspection inspectOneSourceFile(
    const std::filesystem::path& path,
    std::uint32_t file_index,
    IcDialect dialect,
    const std::array<IcSpeciesPolicy, kParticleTypeCount>& species_policy,
    const core::SimulationConfig& config,
    const IcImportOptions& options,
    bool has_authoritative_manifest) {
  if (!std::filesystem::is_regular_file(path)) {
    throw std::runtime_error(
        "IC source is not a regular file: " + path.string());
  }
  Hdf5Handle file(
      H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if (!file.valid()) {
    throw std::runtime_error("failed to open IC member: " + path.string());
  }
  Hdf5Handle header(H5Gopen2(file.get(), "/Header", H5P_DEFAULT));
  if (!header.valid()) {
    throw std::runtime_error("IC member missing /Header: " + path.string());
  }

  SourceFileInspection result;
  result.path = path;
  result.schema = readHeader(header.get());
  result.original_header_attributes = headerAuditText(result.schema);
  checkedCounterAdd(
      result.counters.metadata_bytes_read,
      logicalHeaderPayloadBytes(result.schema), "metadata_bytes_read");
  result.source_size_bytes = std::filesystem::file_size(path);
  result.source_sha256 = sha256Hex(path);
  checkedCounterAdd(
      result.counters.hash_bytes_read, result.source_size_bytes,
      "hash_bytes_read");

  IcManifest dispositions;
  recordRootSchemaDisposition(file.get(), result.schema, dispositions);
  result.preserved_auxiliary_fields =
      std::move(dispositions.preserved_auxiliary_fields);
  result.warnings = std::move(dispositions.warnings);

  const Convention convention = conventionFor(
      dialect, header.get(), config, has_authoritative_manifest);
  result.fields.push_back(inspectHeaderAttribute(
      header.get(), file_index, "MassTable", convention,
      IcFieldSemantics::kExtensive));
  result.fields.push_back(inspectHeaderAttribute(
      header.get(), file_index, "BoxSize", convention,
      IcFieldSemantics::kCoordinate));

  for (std::size_t type = 0; type < kParticleTypeCount; ++type) {
    const std::uint64_t count = result.schema.count_by_type[type];
    if (count == 0U) {
      continue;
    }
    if (species_policy[type] == IcSpeciesPolicy::kReject) {
      throw std::invalid_argument(
          "populated PartType" + std::to_string(type) +
          " has explicit reject policy");
    }
    const std::string group_path = "/PartType" + std::to_string(type);
    Hdf5Handle group(
        H5Gopen2(file.get(), group_path.c_str(), H5P_DEFAULT));
    if (!group.valid()) {
      throw std::runtime_error(
          "missing " + group_path + " in " + path.string());
    }
    std::set<std::string> handled_aliases;
    const auto add = [&](
                         std::string canonical,
                         const std::vector<std::string>& aliases,
                         bool required,
                         IcFieldSemantics semantics,
                         IcVelocityConvention velocity =
                             IcVelocityConvention::kNotVelocity,
                         IcFieldDisposition disposition =
                             IcFieldDisposition::kConverted) {
      const std::string selected = selectAlias(
          group.get(), aliases, required, group_path + "/" + canonical);
      if (selected.empty()) {
        return false;
      }
      handled_aliases.insert(selected);
      IcFieldManifest field = inspectDataset(
          group.get(), file_index, group_path + "/" + canonical, selected,
          convention, semantics, velocity);
      if (field.record_count != count) {
        throw std::runtime_error(
            "dataset record count disagrees with NumPart_ThisFile for " +
            field.dataset_path);
      }
      if ((canonical == "Coordinates" || canonical == "Velocities") &&
          (field.rank != 2U || field.dimensions.size() != 2U ||
           field.dimensions[1] != 3U)) {
        throw std::runtime_error(
            field.dataset_path + " must have dimensions [N,3]");
      }
      if (canonical != "Coordinates" && canonical != "Velocities" &&
          field.rank != 1U) {
        throw std::runtime_error(field.dataset_path + " must have rank 1");
      }
      field.disposition = disposition;
      result.fields.push_back(std::move(field));
      return true;
    };

    add("Coordinates", {"Coordinates", "Position", "POS"}, true,
        IcFieldSemantics::kCoordinate);
    add("Velocities", {"Velocities", "Velocity", "VEL"},
        options.require_velocities, IcFieldSemantics::kVelocity,
        convention.velocity);
    add("ParticleIDs", {"ParticleIDs", "ParticleID", "ID"},
        options.require_particle_ids, IcFieldSemantics::kIdentifier);
    const bool has_masses = add(
        "Masses", {"Masses", "Mass"}, false,
        IcFieldSemantics::kExtensive);
    if (!has_masses && result.schema.mass_table[type] <= 0.0) {
      throw std::runtime_error(
          group_path + " requires Masses because MassTable is zero");
    }

    if (type == 0U) {
      if (!add("InternalEnergy", {"InternalEnergy", "U", "Internal_Energy"},
               false, IcFieldSemantics::kSpecific)) {
        result.defaulted_fields.push_back(
            group_path + "/InternalEnergy=zero");
      }
      if (!add("Density", {"Density", "Rho"}, false,
               IcFieldSemantics::kIntensive)) {
        result.defaulted_fields.push_back(group_path + "/Density=zero");
      }
      if (add("Metallicity", {"Metallicity", "GFM_Metallicity"}, false,
              IcFieldSemantics::kIntensive,
              IcVelocityConvention::kNotVelocity,
              IcFieldDisposition::kDropped)) {
        result.dropped_fields.push_back(
            group_path + "/Metallicity: no gas metallicity lane");
      }
      if (add("SmoothingLength",
              {"SmoothingLength", "Hsml", "Smoothing_Length"}, false,
              IcFieldSemantics::kCoordinate,
              IcVelocityConvention::kNotVelocity,
              IcFieldDisposition::kDropped)) {
        result.dropped_fields.push_back(
            group_path + "/SmoothingLength: no gas smoothing-length lane");
      }
    }
    if (species_policy[type] == IcSpeciesPolicy::kStar) {
      if (!add("StellarFormationTime",
               {"GFM_StellarFormationTime", "StellarFormationTime",
                "BirthTime"},
               false, IcFieldSemantics::kIntensive)) {
        result.defaulted_fields.push_back(
            group_path + "/StellarFormationTime=Header/Time");
      }
      if (!add("InitialMass",
               {"GFM_InitialMass", "InitialMass", "BirthMass"}, false,
               IcFieldSemantics::kExtensive)) {
        result.defaulted_fields.push_back(
            group_path + "/InitialMass=particle_mass");
      }
      if (!add("Metallicity", {"GFM_Metallicity", "Metallicity"}, false,
               IcFieldSemantics::kIntensive)) {
        result.defaulted_fields.push_back(group_path + "/Metallicity=zero");
      }
    }
    if (species_policy[type] == IcSpeciesPolicy::kBlackHole) {
      add("BH_Mass", {"BH_Mass", "BlackHoleMass"}, true,
          IcFieldSemantics::kExtensive);
      if (!add("BH_Mdot", {"BH_Mdot", "BlackHoleAccretionRate"}, false,
               IcFieldSemantics::kIntensive)) {
        result.defaulted_fields.push_back(group_path + "/BH_Mdot=zero");
      }
    }
    if (species_policy[type] == IcSpeciesPolicy::kTracer) {
      add("ParentParticleIDs", {"ParentParticleIDs", "TracerParentIDs"},
          true, IcFieldSemantics::kIdentifier);
      add("InjectionStep", {"InjectionStep"}, true,
          IcFieldSemantics::kIdentifier);
      if (add(
              "HostCellIndex", {"HostCellIndex"}, false,
              IcFieldSemantics::kIdentifier,
              IcVelocityConvention::kNotVelocity,
              IcFieldDisposition::kDropped)) {
        result.dropped_fields.push_back(
            group_path +
            "/HostCellIndex: source-local row is remapped from ParentParticleIDs");
      }
      add("MassFractionOfHost", {"MassFractionOfHost"}, true,
          IcFieldSemantics::kIntensive);
      add("LastHostMass", {"LastHostMass"}, true,
          IcFieldSemantics::kExtensive);
      add("CumulativeExchangedMass", {"CumulativeExchangedMass"}, true,
          IcFieldSemantics::kExtensive);
    }

    for (const std::string& dataset_name : listLinkNames(group.get())) {
      if (handled_aliases.contains(dataset_name)) {
        continue;
      }
      if (objectType(group.get(), dataset_name) != H5O_TYPE_DATASET) {
        throw std::runtime_error(
            group_path + "/" + dataset_name +
            " is not a supported scalar/vector dataset");
      }
      result.fields.push_back(inspectDroppedDataset(
          group.get(), file_index, group_path + "/" + dataset_name,
          dataset_name));
      const std::string warning =
          group_path + "/" + dataset_name +
          ": unsupported dataset explicitly dropped";
      result.dropped_fields.push_back(warning);
      result.warnings.push_back(warning);
    }
  }

  result.counters.bytes_read = result.counters.metadata_bytes_read +
      result.counters.hash_bytes_read;
  return result;
}

[[nodiscard]] Inspection inspectFileSet(
    const std::filesystem::path& requested,
    const core::SimulationConfig& config,
    const IcImportOptions& options) {
  Hdf5Handle first_file(
      H5Fopen(requested.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if (!first_file.valid()) {
    throw std::runtime_error("failed to open IC file: " + requested.string());
  }
  Hdf5Handle first_header(
      H5Gopen2(first_file.get(), "/Header", H5P_DEFAULT));
  if (!first_header.valid()) {
    throw std::runtime_error("IC file missing /Header group");
  }
  const IcSchemaSummary first_schema = readHeader(first_header.get());

  std::vector<std::filesystem::path> files;
  const bool has_authoritative_manifest =
      options.manifest != nullptr && !options.manifest->source_sha256.empty();
  if (has_authoritative_manifest && !options.manifest->source_files.empty()) {
    files = options.manifest->source_files;
  } else {
    files = discoverFiles(requested, first_schema.num_files_per_snapshot);
  }
  if (files.size() != first_schema.num_files_per_snapshot) {
    throw std::runtime_error(
        "manifest/discovery file count disagrees with NumFilesPerSnapshot");
  }

  Inspection inspection;
  IcManifest& manifest = inspection.manifest;
  manifest.dialect =
      config.mode.ic_convention ==
              core::InitialConditionConvention::kChuiCanonicalV1
          ? IcDialect::kChuiCanonicalV1
          : IcDialect::kGadgetArepoBridgeV1;
  manifest.dialect_version = "1";
  manifest.converter_version = "chui_runtime_inspector_v3";
  manifest.source_files = files;
  manifest.num_files_per_snapshot =
      static_cast<std::uint32_t>(files.size());
  manifest.species_policy[2] =
      mapConfiguredPolicy(config.mode.ic_part_type2_policy, 2U);
  manifest.species_policy[3] =
      mapConfiguredPolicy(config.mode.ic_part_type3_policy, 3U);
  if (options.manifest != nullptr) {
    manifest.dialect = options.manifest->dialect;
    manifest.dialect_version = options.manifest->dialect_version;
    manifest.species_policy = options.manifest->species_policy;
  }
  if (manifest.dialect_version != "1") {
    throw std::runtime_error(
        "unsupported IC dialect version: " + manifest.dialect_version);
  }
  if (config.mode.ic_convention ==
          core::InitialConditionConvention::kChuiCanonicalV1 &&
      manifest.dialect != IcDialect::kChuiCanonicalV1) {
    throw std::runtime_error(
        "chui_canonical_v1 configuration requires a canonical CHUÍ "
        "manifest/input");
  }
  if (config.mode.ic_convention ==
          core::InitialConditionConvention::kGadgetArepoBridgeV1 &&
      manifest.dialect != IcDialect::kGadgetArepoBridgeV1) {
    throw std::runtime_error(
        "gadget_arepo_bridge_v1 configuration requires a bridge "
        "manifest/input");
  }

  std::array<std::uint64_t, 6> summed{};
  for (std::size_t file_index = 0; file_index < files.size(); ++file_index) {
    SourceFileInspection source = inspectOneSourceFile(
        files[file_index], static_cast<std::uint32_t>(file_index),
        manifest.dialect, manifest.species_policy, config, options,
        has_authoritative_manifest);
    const IcSchemaSummary& schema = source.schema;
    if (schema.num_files_per_snapshot != files.size() ||
        schema.total_count_by_type != first_schema.total_count_by_type ||
        schema.total_count_high_word != first_schema.total_count_high_word ||
        schema.mass_table != first_schema.mass_table ||
        !nearlyEqual(schema.box_size, first_schema.box_size) ||
        !nearlyEqual(schema.scale_factor, first_schema.scale_factor) ||
        !nearlyEqual(schema.redshift, first_schema.redshift) ||
        !nearlyEqual(schema.omega_matter, first_schema.omega_matter) ||
        !nearlyEqual(schema.omega_lambda, first_schema.omega_lambda) ||
        !nearlyEqual(schema.hubble_param, first_schema.hubble_param)) {
      throw std::runtime_error(
          "inconsistent cosmology, box, epoch, mass table, totals, or "
          "NumFilesPerSnapshot across IC files");
    }
    for (std::size_t type = 0; type < kParticleTypeCount; ++type) {
      if (summed[type] > std::numeric_limits<std::uint64_t>::max() -
              schema.count_by_type[type]) {
        throw std::overflow_error("IC file-set particle count overflow");
      }
      summed[type] += schema.count_by_type[type];
    }
    inspection.schemas.push_back(schema);
    manifest.num_part_this_file.push_back(schema.count_by_type);
    manifest.source_file_sizes_bytes.push_back(source.source_size_bytes);
    manifest.source_sha256.push_back(source.source_sha256);
    manifest.source_provenance_ids.push_back(
        "sha256:" + source.source_sha256);
    manifest.original_header_attributes.push_back(
        std::move(source.original_header_attributes));
    manifest.fields.insert(
        manifest.fields.end(),
        std::make_move_iterator(source.fields.begin()),
        std::make_move_iterator(source.fields.end()));
    manifest.defaulted_fields.insert(
        manifest.defaulted_fields.end(),
        std::make_move_iterator(source.defaulted_fields.begin()),
        std::make_move_iterator(source.defaulted_fields.end()));
    manifest.dropped_fields.insert(
        manifest.dropped_fields.end(),
        std::make_move_iterator(source.dropped_fields.begin()),
        std::make_move_iterator(source.dropped_fields.end()));
    manifest.preserved_auxiliary_fields.insert(
        manifest.preserved_auxiliary_fields.end(),
        std::make_move_iterator(source.preserved_auxiliary_fields.begin()),
        std::make_move_iterator(source.preserved_auxiliary_fields.end()));
    manifest.warnings.insert(
        manifest.warnings.end(),
        std::make_move_iterator(source.warnings.begin()),
        std::make_move_iterator(source.warnings.end()));
    checkedCounterAdd(
        inspection.counters.metadata_bytes_read,
        source.counters.metadata_bytes_read, "metadata_bytes_read");
    checkedCounterAdd(
        inspection.counters.hash_bytes_read,
        source.counters.hash_bytes_read, "hash_bytes_read");
  }

  validateCrossFileSchema(manifest);
  if (summed != first_schema.total_count_by_type) {
    throw std::runtime_error(
        "summed NumPart_ThisFile does not equal reconstructed 64-bit "
        "NumPart_Total");
  }
  manifest.num_part_total = first_schema.total_count_by_type;
  manifest.num_part_total_high_word = first_schema.total_count_high_word;
  manifest.mass_table = first_schema.mass_table;
  manifest.box_size = first_schema.box_size;
  manifest.scale_factor = first_schema.scale_factor;
  manifest.redshift = first_schema.redshift;
  manifest.omega_matter = first_schema.omega_matter;
  manifest.omega_lambda = first_schema.omega_lambda;
  manifest.hubble_param = first_schema.hubble_param;
  manifest.converted_fields.reserve(manifest.fields.size());
  for (const auto& field : manifest.fields) {
    if (field.disposition == IcFieldDisposition::kConverted) {
      manifest.converted_fields.push_back(field.dataset_path);
      if (std::find(
              manifest.conversion_equations.begin(),
              manifest.conversion_equations.end(), field.conversion_equation) ==
          manifest.conversion_equations.end()) {
        manifest.conversion_equations.push_back(field.conversion_equation);
      }
    }
  }

  if (has_authoritative_manifest) {
    const IcManifest& supplied = *options.manifest;
    if (supplied.num_files_per_snapshot !=
            manifest.num_files_per_snapshot ||
        supplied.source_files != manifest.source_files ||
        supplied.source_sha256 != manifest.source_sha256 ||
        supplied.source_file_sizes_bytes !=
            manifest.source_file_sizes_bytes ||
        supplied.num_part_this_file != manifest.num_part_this_file ||
        supplied.num_part_total != manifest.num_part_total ||
        supplied.num_part_total_high_word !=
            manifest.num_part_total_high_word ||
        supplied.mass_table != manifest.mass_table ||
        supplied.species_policy != manifest.species_policy ||
        !nearlyEqual(supplied.box_size, manifest.box_size) ||
        !nearlyEqual(supplied.scale_factor, manifest.scale_factor) ||
        !nearlyEqual(supplied.redshift, manifest.redshift) ||
        !nearlyEqual(supplied.omega_matter, manifest.omega_matter) ||
        !nearlyEqual(supplied.omega_lambda, manifest.omega_lambda) ||
        !nearlyEqual(supplied.hubble_param, manifest.hubble_param)) {
      throw std::runtime_error(
          "supplied IC manifest provenance, scientific header, policies, or "
          "counts do not match inspected source files");
    }
    if (manifest.fields.size() != supplied.fields.size()) {
      throw std::runtime_error(
          "supplied IC manifest field count does not match actual HDF5 "
          "schema");
    }
    for (const auto& expected : supplied.fields) {
      const auto observed = std::find_if(
          manifest.fields.begin(), manifest.fields.end(),
          [&](const IcFieldManifest& actual) {
            return actual.source_file_index == expected.source_file_index &&
                actual.dataset_path == expected.dataset_path &&
                actual.selected_alias == expected.selected_alias;
          });
      if (observed == manifest.fields.end() ||
          observed->scalar_type != expected.scalar_type ||
          observed->scalar_class != expected.scalar_class ||
          observed->byte_width != expected.byte_width ||
          observed->is_signed != expected.is_signed ||
          observed->byte_order != expected.byte_order ||
          observed->rank != expected.rank ||
          observed->dimensions != expected.dimensions ||
          observed->record_count != expected.record_count ||
          observed->disposition != expected.disposition) {
        throw std::runtime_error(
            "supplied IC manifest schema does not match actual HDF5 field " +
            expected.dataset_path);
      }
    }
    manifest = supplied;
  }
  validateIcManifest(manifest);
  inspection.counters.bytes_read =
      inspection.counters.metadata_bytes_read +
      inspection.counters.hash_bytes_read;
  return inspection;
}

[[nodiscard]] const IcFieldManifest* findField(const IcManifest& manifest,std::size_t file_index,std::string_view path){const auto it=std::find_if(manifest.fields.begin(),manifest.fields.end(),[&](const IcFieldManifest& field){return field.source_file_index==file_index&&field.dataset_path==path;});return it==manifest.fields.end()?nullptr:&*it;}
[[nodiscard]] const IcFieldManifest& requireField(const IcManifest& manifest,std::size_t file_index,std::string_view path){const auto* field=findField(manifest,file_index,path);if(field==nullptr)throw std::runtime_error("manifest lacks inspected field "+std::string(path)+" for file "+std::to_string(file_index));return *field;}

void checkedCounterAdd(
    std::uint64_t& destination,
    std::uint64_t value,
    std::string_view counter_name) {
  if (destination > std::numeric_limits<std::uint64_t>::max() - value) {
    throw std::overflow_error(
        "IC " + std::string(counter_name) + " counter overflow");
  }
  destination += value;
}

template <typename T>
[[nodiscard]] std::uint64_t vectorCapacityBytes(
    const std::vector<T>& values) {
  if (values.capacity() >
      std::numeric_limits<std::uint64_t>::max() / sizeof(T)) {
    throw std::overflow_error("IC vector capacity byte count overflow");
  }
  return static_cast<std::uint64_t>(values.capacity()) * sizeof(T);
}

void readChunkDouble(
    hid_t group,
    const IcFieldManifest& field,
    std::size_t start,
    std::size_t count,
    std::size_t components,
    std::vector<double>& out,
    IcImportCounters& counters) {
  Hdf5Handle dataset(
      H5Dopen2(group, field.selected_alias.c_str(), H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error("failed to open " + field.dataset_path);
  }
  Hdf5Handle file_space(H5Dget_space(dataset.get()));
  if (!file_space.valid()) {
    throw std::runtime_error("failed to inspect " + field.dataset_path);
  }
  const int rank = components == 1U ? 1 : 2;
  std::array<hsize_t, 2> offset{
      static_cast<hsize_t>(start), 0U};
  std::array<hsize_t, 2> extent{
      static_cast<hsize_t>(count), static_cast<hsize_t>(components)};
  if (H5Sselect_hyperslab(
          file_space.get(), H5S_SELECT_SET, offset.data(), nullptr,
          extent.data(), nullptr) < 0) {
    throw std::runtime_error("failed hyperslab for " + field.dataset_path);
  }
  Hdf5Handle mem(H5Screate_simple(rank, extent.data(), nullptr));
  if (!mem.valid()) {
    throw std::runtime_error(
        "failed to create memory dataspace for " + field.dataset_path);
  }
  out.resize(count * components);
  if (H5Dread(
          dataset.get(), H5T_NATIVE_DOUBLE, mem.get(), file_space.get(),
          H5P_DEFAULT, out.data()) < 0) {
    throw std::runtime_error("failed chunk read for " + field.dataset_path);
  }
  const std::uint64_t value_count =
      static_cast<std::uint64_t>(count) * components;
  checkedCounterAdd(
      counters.payload_bytes_read, value_count * field.byte_width,
      "payload_bytes_read");
  checkedCounterAdd(
      counters.converted_payload_bytes, value_count * sizeof(double),
      "converted_payload_bytes");
}

void readChunkU64(
    hid_t group,
    const IcFieldManifest& field,
    std::size_t start,
    std::size_t count,
    std::vector<std::uint64_t>& out,
    IcImportCounters& counters) {
  Hdf5Handle dataset(
      H5Dopen2(group, field.selected_alias.c_str(), H5P_DEFAULT));
  if (!dataset.valid()) {
    throw std::runtime_error("failed to open " + field.dataset_path);
  }
  Hdf5Handle file_space(H5Dget_space(dataset.get()));
  if (!file_space.valid()) {
    throw std::runtime_error("failed to inspect " + field.dataset_path);
  }
  hsize_t offset[1]{static_cast<hsize_t>(start)};
  hsize_t extent[1]{static_cast<hsize_t>(count)};
  if (H5Sselect_hyperslab(
          file_space.get(), H5S_SELECT_SET, offset, nullptr, extent,
          nullptr) < 0) {
    throw std::runtime_error("failed hyperslab for " + field.dataset_path);
  }
  Hdf5Handle mem(H5Screate_simple(1, extent, nullptr));
  out.resize(count);
  if (H5Dread(
          dataset.get(), H5T_NATIVE_UINT64, mem.get(), file_space.get(),
          H5P_DEFAULT, out.data()) < 0) {
    throw std::runtime_error("failed ID chunk read for " + field.dataset_path);
  }
  checkedCounterAdd(
      counters.payload_bytes_read,
      static_cast<std::uint64_t>(count) * field.byte_width,
      "payload_bytes_read");
  checkedCounterAdd(
      counters.converted_payload_bytes,
      static_cast<std::uint64_t>(count) * sizeof(std::uint64_t),
      "converted_payload_bytes");
}

void convertValues(
    std::vector<double>& values,
    const IcFieldManifest& field,
    const IcManifest& manifest,
    const core::UnitSystem& target_units) {
  const double factor =
      icFieldConversionMultiplier(field, manifest, target_units);
  for (double& value : values) {
    value *= factor;
    if (!std::isfinite(value)) {
      throw std::runtime_error(
          "non-finite converted value in " + field.dataset_path);
    }
  }
}

struct ParticleRecord {
  std::uint64_t id = 0;
  std::uint32_t species = 0;
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  double vx = 0.0;
  double vy = 0.0;
  double vz = 0.0;
  double mass = 0.0;
  double gas_density = 0.0;
  double gas_internal_energy = 0.0;
  double star_formation = 0.0;
  double star_birth_mass = 0.0;
  double star_metallicity = 0.0;
  double bh_mass = 0.0;
  double bh_mdot = 0.0;
  std::uint64_t tracer_parent = 0;
  std::uint64_t tracer_injection = 0;
  std::uint32_t tracer_host = kInvalidIndex;
  double tracer_fraction = 0.0;
  double tracer_last_host_mass = 0.0;
  double tracer_exchanged_mass = 0.0;
};

[[nodiscard]] std::uint64_t precedingRecordCount(
    const IcManifest& manifest,
    std::size_t file_index,
    std::size_t type_index) {
  std::uint64_t count = 0;
  for (std::size_t file = 0; file < manifest.num_part_this_file.size();
       ++file) {
    for (std::size_t type = 0; type < 6; ++type) {
      if (file == file_index && type == type_index) {
        return count;
      }
      if (count > std::numeric_limits<std::uint64_t>::max() -
                      manifest.num_part_this_file[file][type]) {
        throw std::overflow_error("generated IC ID offset overflow");
      }
      count += manifest.num_part_this_file[file][type];
    }
  }
  throw std::logic_error("invalid file/type offset");
}

[[nodiscard]] double convertedHeaderFieldCode(
    const IcManifest& manifest,
    const IcFieldManifest& field,
    double stored_value,
    const core::UnitSystem& target_units) {
  const double value = stored_value *
      icFieldConversionMultiplier(field, manifest, target_units);
  if (!std::isfinite(value)) {
    throw std::runtime_error(
        "non-finite converted Header value for " + field.dataset_path);
  }
  return value;
}

[[nodiscard]] double normalizePeriodicCoordinate(
    double value,
    double box_size,
    std::string_view axis_name) {
  if (!std::isfinite(value) || !(box_size > 0.0)) {
    throw std::runtime_error(
        "IC " + std::string(axis_name) + " coordinate is non-finite");
  }
  const double tolerance = 1.0e-12 * std::max(1.0, box_size);
  if (value < -tolerance || value > box_size + tolerance) {
    throw std::runtime_error(
        "IC " + std::string(axis_name) +
        " coordinate is outside the periodic box");
  }
  if (value < 0.0 || value >= box_size) {
    return 0.0;
  }
  return value;
}

void validateRecordScientificState(
    ParticleRecord& record,
    IcSpeciesPolicy policy,
    double box_size) {
  record.x = normalizePeriodicCoordinate(record.x, box_size, "x");
  record.y = normalizePeriodicCoordinate(record.y, box_size, "y");
  record.z = normalizePeriodicCoordinate(record.z, box_size, "z");
  if (!std::isfinite(record.vx) || !std::isfinite(record.vy) ||
      !std::isfinite(record.vz)) {
    throw std::runtime_error("IC velocity components must be finite");
  }
  if (record.id == 0U || !(record.mass > 0.0) ||
      !std::isfinite(record.mass)) {
    throw std::runtime_error(
        "IC particle IDs must be nonzero and masses finite/positive");
  }
  if (policy == IcSpeciesPolicy::kGas) {
    if (!std::isfinite(record.gas_density) || record.gas_density < 0.0 ||
        !std::isfinite(record.gas_internal_energy) ||
        record.gas_internal_energy < 0.0) {
      throw std::runtime_error(
          "gas density and internal energy must be finite and non-negative");
    }
  } else if (policy == IcSpeciesPolicy::kStar) {
    if (!std::isfinite(record.star_formation) ||
        record.star_formation < 0.0 ||
        !std::isfinite(record.star_birth_mass) ||
        !(record.star_birth_mass > 0.0) ||
        !std::isfinite(record.star_metallicity) ||
        record.star_metallicity < 0.0) {
      throw std::runtime_error("stellar IC sidecar values are invalid");
    }
  } else if (policy == IcSpeciesPolicy::kBlackHole) {
    if (!std::isfinite(record.bh_mass) || !(record.bh_mass > 0.0) ||
        !std::isfinite(record.bh_mdot) || record.bh_mdot < 0.0) {
      throw std::runtime_error("black-hole sidecar values must be physical");
    }
  } else if (policy == IcSpeciesPolicy::kTracer) {
    if (record.tracer_parent == 0U ||
        !std::isfinite(record.tracer_fraction) ||
        record.tracer_fraction < 0.0 || record.tracer_fraction > 1.0 ||
        !std::isfinite(record.tracer_last_host_mass) ||
        record.tracer_last_host_mass < 0.0 ||
        !std::isfinite(record.tracer_exchanged_mass)) {
      throw std::runtime_error("tracer IC sidecar values are invalid");
    }
  }
}

[[nodiscard]] std::vector<ParticleRecord> readRecordChunk(
    const Inspection& inspection,
    std::size_t file_index,
    std::size_t type_index,
    std::size_t start,
    std::size_t count,
    const core::SimulationConfig& config,
    const IcImportOptions& options,
    IcImportCounters& counters) {
  static_cast<void>(options);
  const IcManifest& manifest = inspection.manifest;
  const auto& path = manifest.source_files[file_index];
  Hdf5Handle file(
      H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  Hdf5Handle group(H5Gopen2(
      file.get(), ("/PartType" + std::to_string(type_index)).c_str(),
      H5P_DEFAULT));
  if (!file.valid() || !group.valid()) {
    throw std::runtime_error("failed to open particle chunk source");
  }
  const std::string prefix =
      "/PartType" + std::to_string(type_index) + "/";
  const core::UnitSystem target = core::makeUnitSystem(
      config.units.length_unit, config.units.mass_unit,
      config.units.velocity_unit);
  const double box_size = convertedHeaderFieldCode(
      manifest, requireField(manifest, file_index, "/Header/BoxSize"),
      manifest.box_size, target);

  std::vector<double> pos;
  std::vector<double> vel;
  std::vector<double> mass;
  std::vector<double> gas_u;
  std::vector<double> gas_rho;
  std::vector<double> star_time;
  std::vector<double> star_birth;
  std::vector<double> star_metal;
  std::vector<double> bh_mass;
  std::vector<double> bh_mdot;
  std::vector<double> tracer_fraction;
  std::vector<double> tracer_last;
  std::vector<double> tracer_exchange;
  std::vector<std::uint64_t> ids;
  std::vector<std::uint64_t> tracer_parent;
  std::vector<std::uint64_t> tracer_injection;
  std::vector<std::uint64_t> tracer_host64;

  const IcFieldManifest& position_field =
      requireField(manifest, file_index, prefix + "Coordinates");
  readChunkDouble(
      group.get(), position_field, start, count, 3U, pos, counters);
  convertValues(pos, position_field, manifest, target);

  if (const auto* field =
          findField(manifest, file_index, prefix + "Velocities")) {
    readChunkDouble(group.get(), *field, start, count, 3U, vel, counters);
    convertValues(vel, *field, manifest, target);
  } else {
    vel.assign(count * 3U, 0.0);
  }

  if (const auto* field =
          findField(manifest, file_index, prefix + "ParticleIDs")) {
    readChunkU64(group.get(), *field, start, count, ids, counters);
  } else {
    ids.resize(count);
    const std::uint64_t base =
        precedingRecordCount(manifest, file_index, type_index) + start + 1U;
    for (std::size_t i = 0; i < count; ++i) {
      ids[i] = base + i;
    }
  }

  if (const auto* field =
          findField(manifest, file_index, prefix + "Masses")) {
    readChunkDouble(group.get(), *field, start, count, 1U, mass, counters);
    convertValues(mass, *field, manifest, target);
  } else {
    mass.assign(count, manifest.mass_table[type_index]);
    convertValues(
        mass, requireField(manifest, file_index, "/Header/MassTable"),
        manifest, target);
  }

  const IcSpeciesPolicy policy = manifest.species_policy[type_index];
  if (policy == IcSpeciesPolicy::kGas) {
    if (const auto* field =
            findField(manifest, file_index, prefix + "InternalEnergy")) {
      readChunkDouble(
          group.get(), *field, start, count, 1U, gas_u, counters);
      convertValues(gas_u, *field, manifest, target);
    } else {
      gas_u.assign(count, 0.0);
    }
    if (const auto* field =
            findField(manifest, file_index, prefix + "Density")) {
      readChunkDouble(
          group.get(), *field, start, count, 1U, gas_rho, counters);
      convertValues(gas_rho, *field, manifest, target);
    } else {
      gas_rho.assign(count, 0.0);
    }
  }

  if (policy == IcSpeciesPolicy::kStar) {
    if (const auto* field = findField(
            manifest, file_index, prefix + "StellarFormationTime")) {
      readChunkDouble(
          group.get(), *field, start, count, 1U, star_time, counters);
    } else {
      star_time.assign(count, manifest.scale_factor);
    }
    if (const auto* field =
            findField(manifest, file_index, prefix + "InitialMass")) {
      readChunkDouble(
          group.get(), *field, start, count, 1U, star_birth, counters);
      convertValues(star_birth, *field, manifest, target);
    } else {
      star_birth = mass;
    }
    if (const auto* field =
            findField(manifest, file_index, prefix + "Metallicity")) {
      readChunkDouble(
          group.get(), *field, start, count, 1U, star_metal, counters);
    } else {
      star_metal.assign(count, 0.0);
    }
  }

  if (policy == IcSpeciesPolicy::kBlackHole) {
    const IcFieldManifest& field =
        requireField(manifest, file_index, prefix + "BH_Mass");
    readChunkDouble(
        group.get(), field, start, count, 1U, bh_mass, counters);
    convertValues(bh_mass, field, manifest, target);
    if (const auto* mdot =
            findField(manifest, file_index, prefix + "BH_Mdot")) {
      readChunkDouble(
          group.get(), *mdot, start, count, 1U, bh_mdot, counters);
      convertValues(bh_mdot, *mdot, manifest, target);
    } else {
      bh_mdot.assign(count, 0.0);
    }
  }

  if (policy == IcSpeciesPolicy::kTracer) {
    readChunkU64(
        group.get(),
        requireField(manifest, file_index, prefix + "ParentParticleIDs"),
        start, count, tracer_parent, counters);
    readChunkU64(
        group.get(),
        requireField(manifest, file_index, prefix + "InjectionStep"),
        start, count, tracer_injection, counters);
    tracer_host64.assign(count, kInvalidIndex);
    const auto& fraction_field =
        requireField(manifest, file_index, prefix + "MassFractionOfHost");
    readChunkDouble(
        group.get(), fraction_field, start, count, 1U, tracer_fraction,
        counters);
    const auto& last_field =
        requireField(manifest, file_index, prefix + "LastHostMass");
    readChunkDouble(
        group.get(), last_field, start, count, 1U, tracer_last, counters);
    convertValues(tracer_last, last_field, manifest, target);
    const auto& exchange_field = requireField(
        manifest, file_index, prefix + "CumulativeExchangedMass");
    readChunkDouble(
        group.get(), exchange_field, start, count, 1U, tracer_exchange,
        counters);
    convertValues(tracer_exchange, exchange_field, manifest, target);
  }

  std::vector<ParticleRecord> records(count);
  for (std::size_t i = 0; i < count; ++i) {
    ParticleRecord& record = records[i];
    record.id = ids[i];
    record.species = speciesTag(policy);
    record.x = pos[i * 3U];
    record.y = pos[i * 3U + 1U];
    record.z = pos[i * 3U + 2U];
    record.vx = vel[i * 3U];
    record.vy = vel[i * 3U + 1U];
    record.vz = vel[i * 3U + 2U];
    record.mass = mass[i];
    if (policy == IcSpeciesPolicy::kGas) {
      record.gas_density = gas_rho[i];
      record.gas_internal_energy = gas_u[i];
    } else if (policy == IcSpeciesPolicy::kStar) {
      record.star_formation = star_time[i];
      record.star_birth_mass = star_birth[i];
      record.star_metallicity = star_metal[i];
    } else if (policy == IcSpeciesPolicy::kBlackHole) {
      record.bh_mass = bh_mass[i];
      record.bh_mdot = bh_mdot[i];
    } else if (policy == IcSpeciesPolicy::kTracer) {
      record.tracer_parent = tracer_parent[i];
      record.tracer_injection = tracer_injection[i];
      if (tracer_host64[i] > std::numeric_limits<std::uint32_t>::max()) {
        throw std::runtime_error(
            "tracer host cell index exceeds uint32 range");
      }
      record.tracer_host = static_cast<std::uint32_t>(tracer_host64[i]);
      record.tracer_fraction = tracer_fraction[i];
      record.tracer_last_host_mass = tracer_last[i];
      record.tracer_exchanged_mass = tracer_exchange[i];
    }
    validateRecordScientificState(record, policy, box_size);
  }

  std::uint64_t staging_bytes = vectorCapacityBytes(records);
  const auto account_vector = [&](const auto& values) {
    checkedCounterAdd(
        staging_bytes, vectorCapacityBytes(values), "peak_staging_bytes");
  };
  account_vector(pos);
  account_vector(vel);
  account_vector(mass);
  account_vector(gas_u);
  account_vector(gas_rho);
  account_vector(star_time);
  account_vector(star_birth);
  account_vector(star_metal);
  account_vector(bh_mass);
  account_vector(bh_mdot);
  account_vector(tracer_fraction);
  account_vector(tracer_last);
  account_vector(tracer_exchange);
  account_vector(ids);
  account_vector(tracer_parent);
  account_vector(tracer_injection);
  account_vector(tracer_host64);
  counters.peak_staging_bytes =
      std::max(counters.peak_staging_bytes, staging_bytes);
  checkedCounterAdd(counters.records_read, count, "records_read");
  checkedCounterAdd(counters.records_converted, count, "records_converted");
  counters.bytes_read = counters.metadata_bytes_read +
      counters.hash_bytes_read + counters.payload_bytes_read;
  return records;
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

void finalizeImportedState(
    core::SimulationState& state,
    const IcManifest& manifest,
    const core::SimulationConfig& config) {
  state.metadata.run_name = config.output.run_name;
  state.metadata.scale_factor = manifest.scale_factor;
  state.species.count_by_species = {};
  for (const std::uint32_t species_tag :
       state.particle_sidecar.species_tag) {
    if (species_tag >= state.species.count_by_species.size()) {
      throw std::runtime_error(
          "IC import produced an out-of-range species tag");
    }
    ++state.species.count_by_species[species_tag];
  }
  state.rebuildSpeciesIndex();
  if (state.cells.size() > 0U) {
    state.refreshGasCellIdentityFromParticleOrder();
  }

  if (state.tracers.size() > 0U) {
    std::unordered_map<std::uint64_t, std::uint32_t> gas_row_by_parent;
    gas_row_by_parent.reserve(state.gas_cells.size());
    for (std::size_t row = 0U; row < state.gas_cells.size(); ++row) {
      const std::uint64_t parent_id =
          state.gas_cells.parent_particle_id[row];
      if (parent_id == 0U ||
          !gas_row_by_parent
               .emplace(parent_id, static_cast<std::uint32_t>(row))
               .second) {
        throw std::runtime_error(
            "IC gas parent identities must be nonzero and unique before "
            "tracer attachment");
      }
    }
    for (std::size_t row = 0U; row < state.tracers.size(); ++row) {
      const auto host = gas_row_by_parent.find(
          state.tracers.parent_particle_id[row]);
      if (host == gas_row_by_parent.end()) {
        throw std::runtime_error(
            "IC tracer parent must resolve to a gas cell on the same final "
            "owner rank");
      }
      state.tracers.host_cell_index[row] = host->second;
    }
  }

  if (!state.validateOwnershipInvariants()) {
    throw std::runtime_error(
        "IC import produced invalid species/sidecar/ownership invariants");
  }
}

[[nodiscard]] double convertedBoxSizeCode(
    const IcManifest& manifest,
    const core::SimulationConfig& config) {
  const core::UnitSystem target = core::makeUnitSystem(
      config.units.length_unit, config.units.mass_unit,
      config.units.velocity_unit);
  return convertedHeaderFieldCode(
      manifest, requireField(manifest, 0U, "/Header/BoxSize"),
      manifest.box_size, target);
}

void validateRuntimeCosmology(const IcManifest& manifest,const core::SimulationConfig& config){const double box_code=convertedBoxSizeCode(manifest,config);const core::UnitSystem target=core::makeUnitSystem(config.units.length_unit,config.units.mass_unit,config.units.velocity_unit);const core::UnitSystem mpc=core::makeUnitSystem("mpc","msun","km_s");const double box_mpc=box_code*target.length_si_per_code/mpc.length_si_per_code;if(!nearlyEqual(manifest.scale_factor,config.numerics.a_begin)||!nearlyEqual(manifest.omega_matter,config.cosmology.omega_matter)||!nearlyEqual(manifest.omega_lambda,config.cosmology.omega_lambda)||!nearlyEqual(manifest.hubble_param,config.cosmology.hubble_param)||!nearlyEqual(box_mpc,config.cosmology.box_size_x_mpc_comoving)||!nearlyEqual(box_mpc,config.cosmology.box_size_y_mpc_comoving)||!nearlyEqual(box_mpc,config.cosmology.box_size_z_mpc_comoving))throw std::runtime_error("IC manifest cosmology/BoxSize/start epoch does not match frozen runtime configuration");}

void validateSerialCountsAndIds(const core::SimulationState& state,const IcManifest& manifest){std::uint64_t expected=0;for(auto count:manifest.num_part_total){if(expected>std::numeric_limits<std::uint64_t>::max()-count)throw std::overflow_error("global particle count overflow");expected+=count;}if(state.particles.size()!=expected)throw std::runtime_error("IC import particle count mismatch");std::vector<std::uint64_t> ids(state.particle_sidecar.particle_id.begin(),state.particle_sidecar.particle_id.end());std::sort(ids.begin(),ids.end());if(std::adjacent_find(ids.begin(),ids.end())!=ids.end())throw std::runtime_error("IC import contains duplicate particle IDs");}

#endif  // COSMOSIM_ENABLE_HDF5

}  // namespace

std::string icSha256Hex(std::string_view value) {
  return sha256Hex(value);
}

std::string icSha256FileHex(const std::filesystem::path& input_path) {
  return sha256Hex(input_path);
}

IcReadResult readGadgetArepoHdf5Ic(
    const std::filesystem::path& ic_path,
    const core::SimulationConfig& config,
    const IcImportOptions& options) {
#if !COSMOSIM_ENABLE_HDF5
  static_cast<void>(ic_path);
  static_cast<void>(config);
  static_cast<void>(options);
  throw std::runtime_error(
      "COSMOSIM_ENABLE_HDF5=OFF: GADGET/AREPO HDF5 IC reader unavailable in this build");
#else
  if (options.chunk_particle_count == 0U) {
    throw std::invalid_argument("chunk_particle_count must be positive");
  }
  Inspection inspection = inspectFileSet(ic_path, config, options);
  if (options.validate_runtime_cosmology) {
    validateRuntimeCosmology(inspection.manifest, config);
  }

  IcReadResult result;
  result.report.counters = inspection.counters;
  result.report.manifest = inspection.manifest;
  result.report.schema = inspection.schemas.front();
  result.report.defaulted_fields = inspection.manifest.defaulted_fields;
  for (const auto& value : result.report.defaulted_fields) {
    const auto separator = value.find('=');
    result.report.missing_optional_fields.push_back(
        value.substr(0U, separator));
  }
  result.report.unsupported_fields = inspection.manifest.dropped_fields;

  for (std::size_t file_index = 0U;
       file_index < inspection.manifest.source_files.size(); ++file_index) {
    ++result.report.counters.files_assigned;
    for (std::size_t type_index = 0U; type_index < 6U; ++type_index) {
      const std::size_t total = static_cast<std::size_t>(
          inspection.manifest.num_part_this_file[file_index][type_index]);
      for (std::size_t start = 0U; start < total;
           start += options.chunk_particle_count) {
        const std::size_t count =
            std::min(options.chunk_particle_count, total - start);
        auto records = readRecordChunk(
            inspection, file_index, type_index, start, count, config, options,
            result.report.counters);
        appendRecords(result.state, records, 0U);
        ++result.report.counters.chunks_assigned;
        result.report.counters.records_routed += records.size();
      }
    }
  }

  validateSerialCountsAndIds(result.state, inspection.manifest);
  finalizeImportedState(result.state, inspection.manifest, config);
  result.report.counters.final_local_particle_count =
      result.state.particles.size();
  result.report.counters.final_local_gas_cell_count = result.state.cells.size();
  result.report.counters.final_local_star_count =
      result.state.star_particles.size();
  result.report.counters.final_local_black_hole_count =
      result.state.black_holes.size();
  result.report.counters.final_local_tracer_count = result.state.tracers.size();
  result.report.counters.bytes_read =
      result.report.counters.metadata_bytes_read +
      result.report.counters.hash_bytes_read +
      result.report.counters.payload_bytes_read;
  return result;
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
  MPI_Allreduce(
      &local_failed, &failed_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (failed_count == 0) {
    return {};
  }
  int candidate = local_failure ? mpi_context.worldRank()
                                : mpi_context.worldSize();
  int failure_rank = mpi_context.worldSize();
  MPI_Allreduce(
      &candidate, &failure_rank, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  static constexpr std::size_t kMaximumMessageBytes = 4095U;
  std::array<char, kMaximumMessageBytes + 1U> buffer{};
  std::uint32_t length = 0U;
  if (mpi_context.worldRank() == failure_rank) {
    std::string message;
    try {
      std::rethrow_exception(local_failure);
    } catch (const std::exception& error) {
      message = error.what();
    } catch (...) {
      message = "unknown non-standard exception";
    }
    length = static_cast<std::uint32_t>(
        std::min(message.size(), kMaximumMessageBytes));
    std::copy_n(message.data(), length, buffer.data());
  }
  MPI_Bcast(&length, 1, MPI_UINT32_T, failure_rank, MPI_COMM_WORLD);
  MPI_Bcast(
      buffer.data(), static_cast<int>(buffer.size()), MPI_CHAR, failure_rank,
      MPI_COMM_WORLD);
  return std::string(phase) + " failed on rank " +
      std::to_string(failure_rank) + ": " +
      std::string(buffer.data(), buffer.data() + length);
}

void injectIcTestFault(
    const parallel::MpiContext& mpi_context,
    std::string_view phase) {
#if COSMOSIM_ENABLE_TESTS
  const char* raw = std::getenv("COSMOSIM_IC_TEST_FAULT");
  if (raw == nullptr || *raw == '\0') {
    return;
  }
  const std::string specification(raw);
  const std::size_t separator = specification.rfind(':');
  if (separator == std::string::npos) {
    return;
  }
  const std::string configured_phase = specification.substr(0, separator);
  int configured_rank = -1;
  try {
    configured_rank = std::stoi(specification.substr(separator + 1U));
  } catch (...) {
    return;
  }
  if (configured_phase == phase &&
      configured_rank == mpi_context.worldRank()) {
    throw std::runtime_error(
        "test-only injected IC failure at phase " + std::string(phase));
  }
#else
  static_cast<void>(mpi_context);
  static_cast<void>(phase);
#endif
}


void mutateIcTestRoute(
    const parallel::MpiContext& mpi_context,
    std::vector<std::vector<std::uint8_t>>& per_rank) {
#if COSMOSIM_ENABLE_TESTS
  static bool mutation_applied = false;
  if (mutation_applied) {
    return;
  }
  const char* raw = std::getenv("COSMOSIM_IC_TEST_ROUTE_MUTATION");
  if (raw == nullptr || *raw == '\0') {
    return;
  }
  const std::string specification(raw);
  const std::size_t separator = specification.rfind(':');
  if (separator == std::string::npos) {
    return;
  }
  int configured_rank = -1;
  try {
    configured_rank = std::stoi(specification.substr(separator + 1U));
  } catch (...) {
    return;
  }
  if (configured_rank != mpi_context.worldRank()) {
    return;
  }
  const std::string operation = specification.substr(0U, separator);
  auto bucket = std::find_if(
      per_rank.begin(), per_rank.end(), [](const auto& candidate) {
        return candidate.size() >= kWireRecordBytes;
      });
  if (bucket == per_rank.end()) {
    return;
  }
  if (operation == "drop") {
    bucket->resize(bucket->size() - kWireRecordBytes);
  } else if (operation == "duplicate") {
    const std::size_t record_begin = bucket->size() - kWireRecordBytes;
    bucket->insert(
        bucket->end(), bucket->begin() + static_cast<std::ptrdiff_t>(record_begin),
        bucket->end());
  } else {
    return;
  }
  mutation_applied = true;
#else
  static_cast<void>(mpi_context);
  static_cast<void>(per_rank);
#endif
}

template <typename T, typename Function>
[[nodiscard]] T runCollectivePhase(
    const parallel::MpiContext& mpi_context,
    std::string_view phase,
    Function&& function) {
  std::optional<T> value;
  std::exception_ptr local_failure;
  try {
    value.emplace(std::invoke(std::forward<Function>(function)));
  } catch (...) {
    local_failure = std::current_exception();
  }
  const std::string failure =
      collectiveFailureMessage(mpi_context, local_failure, phase);
  if (!failure.empty()) {
    throw std::runtime_error(failure);
  }
  return std::move(*value);
}

template <typename Function>
void runCollectivePhaseVoid(
    const parallel::MpiContext& mpi_context,
    std::string_view phase,
    Function&& function) {
  std::exception_ptr local_failure;
  try {
    std::invoke(std::forward<Function>(function));
  } catch (...) {
    local_failure = std::current_exception();
  }
  const std::string failure =
      collectiveFailureMessage(mpi_context, local_failure, phase);
  if (!failure.empty()) {
    throw std::runtime_error(failure);
  }
}

[[nodiscard]] std::string broadcastRootString(
    const parallel::MpiContext& mpi_context,
    std::string root_value) {
  std::uint64_t length = mpi_context.isRoot() ? root_value.size() : 0U;
  MPI_Bcast(&length, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  runCollectivePhaseVoid(
      mpi_context, "root-string receive allocation", [&]() {
        if (length > std::numeric_limits<std::size_t>::max()) {
          throw std::overflow_error("root string exceeds local size_t");
        }
        root_value.resize(static_cast<std::size_t>(length));
      });
  static constexpr std::uint64_t kBroadcastChunk = 64U * 1024U * 1024U;
  for (std::uint64_t offset = 0U; offset < length;
       offset += kBroadcastChunk) {
    const int chunk = static_cast<int>(
        std::min(kBroadcastChunk, length - offset));
    MPI_Bcast(
        root_value.data() + static_cast<std::size_t>(offset), chunk, MPI_CHAR,
        0, MPI_COMM_WORLD);
  }
  return root_value;
}

struct AlltoallLayout {
  std::vector<int> send_counts;
  std::vector<int> receive_counts;
  std::vector<int> send_displacements;
  std::vector<int> receive_displacements;
  std::uint64_t send_total = 0U;
  std::uint64_t receive_total = 0U;
};

[[nodiscard]] std::vector<std::uint8_t> alltoallBytes(
    const parallel::MpiContext& mpi_context,
    const std::vector<std::vector<std::uint8_t>>& per_rank,
    std::uint64_t& sent_bytes,
    std::uint64_t& received_bytes,
    std::uint64_t* exchange_peak_bytes = nullptr) {
  const int world_size = mpi_context.worldSize();
  AlltoallLayout layout = runCollectivePhase<AlltoallLayout>(
      mpi_context, "IC all-to-all send-layout preparation", [&]() {
        if (per_rank.size() != static_cast<std::size_t>(world_size)) {
          throw std::invalid_argument(
              "IC all-to-all requires one send bucket per rank");
        }
        AlltoallLayout prepared;
        prepared.send_counts.assign(static_cast<std::size_t>(world_size), 0);
        prepared.receive_counts.assign(
            static_cast<std::size_t>(world_size), 0);
        prepared.send_displacements.assign(
            static_cast<std::size_t>(world_size), 0);
        prepared.receive_displacements.assign(
            static_cast<std::size_t>(world_size), 0);
        for (int rank = 0; rank < world_size; ++rank) {
          const std::size_t size = per_rank[static_cast<std::size_t>(rank)].size();
          if (size > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
              prepared.send_total >
                  static_cast<std::uint64_t>(std::numeric_limits<int>::max()) -
                      size) {
            throw std::overflow_error(
                "bounded IC send payload exceeds MPI int count/displacement");
          }
          prepared.send_displacements[static_cast<std::size_t>(rank)] =
              static_cast<int>(prepared.send_total);
          prepared.send_counts[static_cast<std::size_t>(rank)] =
              static_cast<int>(size);
          prepared.send_total += size;
        }
        return prepared;
      });

  MPI_Alltoall(
      layout.send_counts.data(), 1, MPI_INT, layout.receive_counts.data(), 1,
      MPI_INT, MPI_COMM_WORLD);

  layout = runCollectivePhase<AlltoallLayout>(
      mpi_context, "IC all-to-all receive-layout preparation", [&]() {
        for (int rank = 0; rank < world_size; ++rank) {
          const int count =
              layout.receive_counts[static_cast<std::size_t>(rank)];
          if (count < 0 ||
              layout.receive_total >
                  static_cast<std::uint64_t>(std::numeric_limits<int>::max()) -
                      static_cast<std::uint64_t>(count)) {
            throw std::overflow_error(
                "bounded IC receive payload exceeds MPI int "
                "count/displacement");
          }
          layout.receive_displacements[static_cast<std::size_t>(rank)] =
              static_cast<int>(layout.receive_total);
          layout.receive_total += static_cast<std::uint64_t>(count);
        }
        return std::move(layout);
      });

  struct ExchangeBuffers {
    std::vector<std::uint8_t> send;
    std::vector<std::uint8_t> receive;
  };
  ExchangeBuffers buffers = runCollectivePhase<ExchangeBuffers>(
      mpi_context, "IC all-to-all buffer allocation", [&]() {
        ExchangeBuffers prepared;
        prepared.send.resize(static_cast<std::size_t>(layout.send_total));
        prepared.receive.resize(
            static_cast<std::size_t>(layout.receive_total));
        for (int rank = 0; rank < world_size; ++rank) {
          const auto& bucket = per_rank[static_cast<std::size_t>(rank)];
          std::copy(
              bucket.begin(), bucket.end(),
              prepared.send.begin() +
                  layout.send_displacements[static_cast<std::size_t>(rank)]);
        }
        return prepared;
      });

  MPI_Alltoallv(
      buffers.send.data(), layout.send_counts.data(),
      layout.send_displacements.data(), MPI_BYTE, buffers.receive.data(),
      layout.receive_counts.data(), layout.receive_displacements.data(),
      MPI_BYTE, MPI_COMM_WORLD);
  sent_bytes += layout.send_total;
  received_bytes += layout.receive_total;
  if (exchange_peak_bytes != nullptr) {
    const std::uint64_t layout_bytes =
        static_cast<std::uint64_t>(layout.send_counts.capacity() +
                                   layout.receive_counts.capacity() +
                                   layout.send_displacements.capacity() +
                                   layout.receive_displacements.capacity()) *
        sizeof(int);
    *exchange_peak_bytes = layout_bytes +
        static_cast<std::uint64_t>(buffers.send.capacity()) +
        static_cast<std::uint64_t>(buffers.receive.capacity());
  }
  return std::move(buffers.receive);
}

[[nodiscard]] IcManifest makeTransferManifest(
    const SourceFileInspection& source,
    IcDialect dialect,
    const std::array<IcSpeciesPolicy, kParticleTypeCount>& species_policy) {
  IcManifest manifest;
  manifest.dialect = dialect;
  manifest.dialect_version = "1";
  manifest.converter_version = "chui_distributed_inspector_fragment_v3";
  manifest.source_files = {source.path};
  manifest.source_provenance_ids = {"sha256:" + source.source_sha256};
  manifest.source_file_sizes_bytes = {source.source_size_bytes};
  manifest.source_sha256 = {source.source_sha256};
  manifest.original_header_attributes = {source.original_header_attributes};
  manifest.num_files_per_snapshot = 1U;
  manifest.num_part_this_file = {source.schema.count_by_type};
  manifest.num_part_total = source.schema.count_by_type;
  for (std::size_t type = 0; type < kParticleTypeCount; ++type) {
    manifest.num_part_total_high_word[type] =
        manifest.num_part_total[type] >> 32U;
  }
  manifest.mass_table = source.schema.mass_table;
  manifest.box_size = source.schema.box_size;
  manifest.scale_factor = source.schema.scale_factor;
  manifest.redshift = source.schema.redshift;
  manifest.omega_matter = source.schema.omega_matter;
  manifest.omega_lambda = source.schema.omega_lambda;
  manifest.hubble_param = source.schema.hubble_param;
  manifest.species_policy = species_policy;
  manifest.fields = source.fields;
  for (IcFieldManifest& field : manifest.fields) {
    field.source_file_index = 0U;
    if (field.disposition == IcFieldDisposition::kConverted) {
      manifest.converted_fields.push_back(field.dataset_path);
      if (std::find(
              manifest.conversion_equations.begin(),
              manifest.conversion_equations.end(), field.conversion_equation) ==
          manifest.conversion_equations.end()) {
        manifest.conversion_equations.push_back(field.conversion_equation);
      }
    }
  }
  manifest.defaulted_fields = source.defaulted_fields;
  manifest.dropped_fields = source.dropped_fields;
  manifest.preserved_auxiliary_fields = source.preserved_auxiliary_fields;
  manifest.warnings = source.warnings;
  validateIcManifest(manifest);
  return manifest;
}

void appendSchemaSummary(
    std::vector<std::uint8_t>& output,
    const IcSchemaSummary& schema) {
  for (std::uint64_t value : schema.count_by_type) {
    appendLe64(output, value);
  }
  for (std::uint64_t value : schema.total_count_by_type) {
    appendLe64(output, value);
  }
  for (std::uint32_t value : schema.total_count_high_word) {
    appendLe32(output, value);
  }
  for (double value : schema.mass_table) {
    appendDouble(output, value);
  }
  appendLe32(output, schema.num_files_per_snapshot);
  appendDouble(output, schema.box_size);
  appendDouble(output, schema.scale_factor);
  appendDouble(output, schema.redshift);
  appendDouble(output, schema.omega_matter);
  appendDouble(output, schema.omega_lambda);
  appendDouble(output, schema.hubble_param);
}

[[nodiscard]] IcSchemaSummary readSchemaSummary(
    std::span<const std::uint8_t> input,
    std::size_t& offset) {
  IcSchemaSummary schema;
  for (std::uint64_t& value : schema.count_by_type) {
    value = readLe64(input, offset);
  }
  for (std::uint64_t& value : schema.total_count_by_type) {
    value = readLe64(input, offset);
  }
  for (std::uint32_t& value : schema.total_count_high_word) {
    value = readLe32(input, offset);
  }
  for (double& value : schema.mass_table) {
    value = readDouble(input, offset);
  }
  schema.num_files_per_snapshot = readLe32(input, offset);
  schema.box_size = readDouble(input, offset);
  schema.scale_factor = readDouble(input, offset);
  schema.redshift = readDouble(input, offset);
  schema.omega_matter = readDouble(input, offset);
  schema.omega_lambda = readDouble(input, offset);
  schema.hubble_param = readDouble(input, offset);
  return schema;
}

void appendLengthPrefixedString(
    std::vector<std::uint8_t>& output,
    std::string_view value) {
  appendLe64(output, value.size());
  output.insert(output.end(), value.begin(), value.end());
}

[[nodiscard]] std::string readLengthPrefixedString(
    std::span<const std::uint8_t> input,
    std::size_t& offset) {
  const std::uint64_t length = readLe64(input, offset);
  if (length > input.size() - offset) {
    throw std::runtime_error("truncated distributed IC metadata string");
  }
  std::string value(
      reinterpret_cast<const char*>(input.data() + offset),
      static_cast<std::size_t>(length));
  offset += static_cast<std::size_t>(length);
  return value;
}

[[nodiscard]] std::vector<std::uint8_t> gatherRootBytes(
    const parallel::MpiContext& mpi_context,
    const std::vector<std::uint8_t>& local_bytes) {
  const int local_count = runCollectivePhase<int>(
      mpi_context, "IC metadata gather local-count validation", [&]() {
        if (local_bytes.size() >
            static_cast<std::size_t>(std::numeric_limits<int>::max())) {
          throw std::overflow_error(
              "IC metadata fragment exceeds MPI int count");
        }
        return static_cast<int>(local_bytes.size());
      });
  std::vector<int> counts;
  if (mpi_context.isRoot()) {
    counts.resize(static_cast<std::size_t>(mpi_context.worldSize()));
  }
  MPI_Gather(
      &local_count, 1, MPI_INT,
      mpi_context.isRoot() ? counts.data() : nullptr, 1, MPI_INT, 0,
      MPI_COMM_WORLD);

  struct RootGatherLayout {
    std::vector<int> displacements;
    std::vector<std::uint8_t> bytes;
  };
  RootGatherLayout layout = runCollectivePhase<RootGatherLayout>(
      mpi_context, "IC metadata gather root-layout preparation", [&]() {
        RootGatherLayout prepared;
        if (!mpi_context.isRoot()) {
          return prepared;
        }
        prepared.displacements.resize(counts.size());
        std::uint64_t total = 0U;
        for (std::size_t rank = 0; rank < counts.size(); ++rank) {
          if (counts[rank] < 0 ||
              total >
                  static_cast<std::uint64_t>(std::numeric_limits<int>::max()) -
                      static_cast<std::uint64_t>(counts[rank])) {
            throw std::overflow_error(
                "IC metadata gather exceeds MPI int displacement range");
          }
          prepared.displacements[rank] = static_cast<int>(total);
          total += static_cast<std::uint64_t>(counts[rank]);
        }
        prepared.bytes.resize(static_cast<std::size_t>(total));
        return prepared;
      });

  MPI_Gatherv(
      local_bytes.empty() ? nullptr : local_bytes.data(), local_count, MPI_BYTE,
      mpi_context.isRoot() && !layout.bytes.empty() ? layout.bytes.data() : nullptr,
      mpi_context.isRoot() ? counts.data() : nullptr,
      mpi_context.isRoot() ? layout.displacements.data() : nullptr,
      MPI_BYTE, 0, MPI_COMM_WORLD);
  return std::move(layout.bytes);
}

[[nodiscard]] std::string encodeDiscoveredPaths(
    const std::vector<std::filesystem::path>& paths) {
  std::vector<std::uint8_t> bytes;
  appendLe64(bytes, paths.size());
  for (const auto& path : paths) {
    appendLengthPrefixedString(bytes, path.lexically_normal().generic_string());
  }
  return std::string(
      reinterpret_cast<const char*>(bytes.data()), bytes.size());
}

[[nodiscard]] std::vector<std::filesystem::path> decodeDiscoveredPaths(
    std::string_view encoded) {
  const auto* begin = reinterpret_cast<const std::uint8_t*>(encoded.data());
  std::span<const std::uint8_t> bytes(begin, encoded.size());
  std::size_t offset = 0U;
  const std::uint64_t count = readLe64(bytes, offset);
  if (count > std::numeric_limits<std::size_t>::max()) {
    throw std::overflow_error("discovered IC file count exceeds size_t");
  }
  std::vector<std::filesystem::path> paths;
  paths.reserve(static_cast<std::size_t>(count));
  for (std::uint64_t index = 0U; index < count; ++index) {
    paths.emplace_back(readLengthPrefixedString(bytes, offset));
  }
  if (offset != bytes.size()) {
    throw std::runtime_error("distributed IC path metadata has trailing bytes");
  }
  return paths;
}

void appendManifestLists(IcManifest& destination, IcManifest&& source) {
  destination.fields.insert(
      destination.fields.end(),
      std::make_move_iterator(source.fields.begin()),
      std::make_move_iterator(source.fields.end()));
  destination.defaulted_fields.insert(
      destination.defaulted_fields.end(),
      std::make_move_iterator(source.defaulted_fields.begin()),
      std::make_move_iterator(source.defaulted_fields.end()));
  destination.dropped_fields.insert(
      destination.dropped_fields.end(),
      std::make_move_iterator(source.dropped_fields.begin()),
      std::make_move_iterator(source.dropped_fields.end()));
  destination.preserved_auxiliary_fields.insert(
      destination.preserved_auxiliary_fields.end(),
      std::make_move_iterator(source.preserved_auxiliary_fields.begin()),
      std::make_move_iterator(source.preserved_auxiliary_fields.end()));
  destination.warnings.insert(
      destination.warnings.end(),
      std::make_move_iterator(source.warnings.begin()),
      std::make_move_iterator(source.warnings.end()));
}

[[nodiscard]] Inspection inspectFileSetDistributed(
    const std::filesystem::path& requested,
    const core::SimulationConfig& config,
    const parallel::MpiContext& mpi_context,
    const IcImportOptions& options) {
  const bool has_authoritative_manifest =
      options.manifest != nullptr && !options.manifest->source_sha256.empty();
  const IcDialect dialect =
      config.mode.ic_convention ==
              core::InitialConditionConvention::kChuiCanonicalV1
          ? IcDialect::kChuiCanonicalV1
          : IcDialect::kGadgetArepoBridgeV1;
  std::array<IcSpeciesPolicy, kParticleTypeCount> species_policy{
      IcSpeciesPolicy::kGas,
      IcSpeciesPolicy::kDarkMatter,
      mapConfiguredPolicy(config.mode.ic_part_type2_policy, 2U),
      mapConfiguredPolicy(config.mode.ic_part_type3_policy, 3U),
      IcSpeciesPolicy::kStar,
      IcSpeciesPolicy::kBlackHole};
  if (options.manifest != nullptr) {
    species_policy = options.manifest->species_policy;
  }

  std::string root_paths = runCollectivePhase<std::string>(
      mpi_context, "IC distributed file discovery", [&]() {
        if (!mpi_context.isRoot()) {
          return std::string{};
        }
        Hdf5Handle file(
            H5Fopen(requested.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
        if (!file.valid()) {
          throw std::runtime_error(
              "failed to open first IC source for distributed discovery");
        }
        Hdf5Handle header(H5Gopen2(file.get(), "/Header", H5P_DEFAULT));
        if (!header.valid()) {
          throw std::runtime_error("IC source is missing /Header");
        }
        std::uint32_t file_count = 0U;
        readAttributeU32(header.get(), "NumFilesPerSnapshot", file_count);
        std::vector<std::filesystem::path> paths;
        if (has_authoritative_manifest &&
            !options.manifest->source_files.empty()) {
          paths = options.manifest->source_files;
        } else {
          paths = discoverFiles(requested, file_count);
        }
        if (paths.size() != file_count) {
          throw std::runtime_error(
              "distributed IC discovery disagrees with "
              "NumFilesPerSnapshot");
        }
        return encodeDiscoveredPaths(paths);
      });
  root_paths = broadcastRootString(mpi_context, std::move(root_paths));
  const std::vector<std::filesystem::path> paths =
      runCollectivePhase<std::vector<std::filesystem::path>>(
          mpi_context, "IC distributed path metadata decoding", [&]() {
            return decodeDiscoveredPaths(root_paths);
          });

  struct LocalInspectionCollection {
    std::vector<std::pair<std::uint32_t, SourceFileInspection>> sources;
    IcImportCounters counters;
  };
  LocalInspectionCollection local =
      runCollectivePhase<LocalInspectionCollection>(
          mpi_context, "IC assigned-file inspection and hashing", [&]() {
            LocalInspectionCollection collection;
            for (std::size_t file_index = 0; file_index < paths.size();
                 ++file_index) {
              if (static_cast<int>(file_index %
                                   static_cast<std::size_t>(
                                       mpi_context.worldSize())) !=
                  mpi_context.worldRank()) {
                continue;
              }
              SourceFileInspection source = inspectOneSourceFile(
                  paths[file_index], 0U, dialect, species_policy, config,
                  options, has_authoritative_manifest);
              checkedCounterAdd(
                  collection.counters.metadata_bytes_read,
                  source.counters.metadata_bytes_read,
                  "metadata_bytes_read");
              checkedCounterAdd(
                  collection.counters.hash_bytes_read,
                  source.counters.hash_bytes_read, "hash_bytes_read");
              collection.sources.emplace_back(
                  static_cast<std::uint32_t>(file_index), std::move(source));
            }
            collection.counters.files_assigned = collection.sources.size();
            collection.counters.bytes_read =
                collection.counters.metadata_bytes_read +
                collection.counters.hash_bytes_read;
            return collection;
          });

  std::vector<std::uint8_t> local_metadata =
      runCollectivePhase<std::vector<std::uint8_t>>(
          mpi_context, "IC manifest-fragment serialization", [&]() {
            std::vector<std::uint8_t> bytes;
            appendLe64(bytes, local.sources.size());
            for (const auto& [file_index, source] : local.sources) {
              appendLe32(bytes, file_index);
              appendSchemaSummary(bytes, source.schema);
              const std::string json = serializeIcManifestJson(
                  makeTransferManifest(source, dialect, species_policy));
              appendLengthPrefixedString(bytes, json);
            }
            return bytes;
          });
  checkedCounterAdd(
      local.counters.manifest_metadata_bytes_communicated,
      local_metadata.size(), "manifest_metadata_bytes_communicated");
  const std::vector<std::uint8_t> gathered =
      gatherRootBytes(mpi_context, local_metadata);
  if (mpi_context.isRoot()) {
    checkedCounterAdd(
        local.counters.manifest_metadata_bytes_communicated, gathered.size(),
        "manifest_metadata_bytes_communicated");
  }

  Inspection root_inspection = runCollectivePhase<Inspection>(
      mpi_context, "IC distributed manifest assembly", [&]() {
        Inspection assembled;
        assembled.counters = local.counters;
        if (!mpi_context.isRoot()) {
          return assembled;
        }
        std::vector<std::optional<std::pair<IcSchemaSummary, IcManifest>>>
            fragments(paths.size());
        std::size_t offset = 0U;
        for (int rank = 0; rank < mpi_context.worldSize(); ++rank) {
          if (offset >= gathered.size()) {
            throw std::runtime_error(
                "distributed IC metadata is missing a rank fragment");
          }
          const std::uint64_t entry_count = readLe64(gathered, offset);
          for (std::uint64_t entry = 0U; entry < entry_count; ++entry) {
            const std::uint32_t file_index = readLe32(gathered, offset);
            if (file_index >= fragments.size() || fragments[file_index]) {
              throw std::runtime_error(
                  "distributed IC file fragment is duplicated or out of range");
            }
            IcSchemaSummary schema = readSchemaSummary(gathered, offset);
            IcManifest fragment = deserializeIcManifestJson(
                readLengthPrefixedString(gathered, offset));
            fragments[file_index].emplace(
                std::move(schema), std::move(fragment));
          }
        }
        if (offset != gathered.size()) {
          throw std::runtime_error(
              "distributed IC metadata has trailing bytes");
        }
        if (std::any_of(
                fragments.begin(), fragments.end(),
                [](const auto& fragment) { return !fragment.has_value(); })) {
          throw std::runtime_error(
              "distributed IC inspection did not cover every source file");
        }

        const IcSchemaSummary& baseline = fragments.front()->first;
        IcManifest& manifest = assembled.manifest;
        manifest.dialect = dialect;
        manifest.dialect_version = "1";
        manifest.converter_version = "chui_distributed_inspector_v3";
        manifest.num_files_per_snapshot =
            static_cast<std::uint32_t>(paths.size());
        manifest.source_files = paths;
        manifest.species_policy = species_policy;
        manifest.num_part_total = baseline.total_count_by_type;
        manifest.num_part_total_high_word = baseline.total_count_high_word;
        manifest.mass_table = baseline.mass_table;
        manifest.box_size = baseline.box_size;
        manifest.scale_factor = baseline.scale_factor;
        manifest.redshift = baseline.redshift;
        manifest.omega_matter = baseline.omega_matter;
        manifest.omega_lambda = baseline.omega_lambda;
        manifest.hubble_param = baseline.hubble_param;

        std::array<std::uint64_t, kParticleTypeCount> summed{};
        for (std::size_t file_index = 0; file_index < fragments.size();
             ++file_index) {
          IcSchemaSummary schema = fragments[file_index]->first;
          IcManifest fragment = std::move(fragments[file_index]->second);
          if (schema.num_files_per_snapshot != paths.size() ||
              schema.total_count_by_type != baseline.total_count_by_type ||
              schema.total_count_high_word != baseline.total_count_high_word ||
              schema.mass_table != baseline.mass_table ||
              !nearlyEqual(schema.box_size, baseline.box_size) ||
              !nearlyEqual(schema.scale_factor, baseline.scale_factor) ||
              !nearlyEqual(schema.redshift, baseline.redshift) ||
              !nearlyEqual(schema.omega_matter, baseline.omega_matter) ||
              !nearlyEqual(schema.omega_lambda, baseline.omega_lambda) ||
              !nearlyEqual(schema.hubble_param, baseline.hubble_param)) {
            throw std::runtime_error(
                "distributed IC source files disagree on scientific header "
                "or totals");
          }
          for (std::size_t type = 0; type < kParticleTypeCount; ++type) {
            if (summed[type] >
                std::numeric_limits<std::uint64_t>::max() -
                    schema.count_by_type[type]) {
              throw std::overflow_error(
                  "distributed IC particle-count sum overflow");
            }
            summed[type] += schema.count_by_type[type];
          }
          assembled.schemas.push_back(schema);
          manifest.num_part_this_file.push_back(schema.count_by_type);
          manifest.source_file_sizes_bytes.push_back(
              fragment.source_file_sizes_bytes.front());
          manifest.source_sha256.push_back(fragment.source_sha256.front());
          manifest.source_provenance_ids.push_back(
              fragment.source_provenance_ids.front());
          manifest.original_header_attributes.push_back(
              fragment.original_header_attributes.front());
          for (IcFieldManifest& field : fragment.fields) {
            field.source_file_index = static_cast<std::uint32_t>(file_index);
          }
          appendManifestLists(manifest, std::move(fragment));
        }
        if (summed != manifest.num_part_total) {
          throw std::runtime_error(
              "distributed IC per-file counts do not equal declared totals");
        }
        validateCrossFileSchema(manifest);
        for (const IcFieldManifest& field : manifest.fields) {
          if (field.disposition == IcFieldDisposition::kConverted) {
            manifest.converted_fields.push_back(field.dataset_path);
            if (std::find(
                    manifest.conversion_equations.begin(),
                    manifest.conversion_equations.end(),
                    field.conversion_equation) ==
                manifest.conversion_equations.end()) {
              manifest.conversion_equations.push_back(
                  field.conversion_equation);
            }
          }
        }
        if (has_authoritative_manifest) {
          const IcManifest& supplied = *options.manifest;
          if (supplied.source_files != manifest.source_files ||
              supplied.source_sha256 != manifest.source_sha256 ||
              supplied.source_file_sizes_bytes !=
                  manifest.source_file_sizes_bytes ||
              supplied.num_part_this_file != manifest.num_part_this_file ||
              supplied.num_part_total != manifest.num_part_total ||
              supplied.species_policy != manifest.species_policy ||
              supplied.fields.size() != manifest.fields.size()) {
            throw std::runtime_error(
                "supplied distributed IC manifest does not match inspected "
                "source provenance or schema");
          }
          for (const IcFieldManifest& expected : supplied.fields) {
            const auto observed = std::find_if(
                manifest.fields.begin(), manifest.fields.end(),
                [&](const IcFieldManifest& actual) {
                  return actual.source_file_index ==
                             expected.source_file_index &&
                      actual.dataset_path == expected.dataset_path &&
                      actual.selected_alias == expected.selected_alias;
                });
            if (observed == manifest.fields.end() ||
                observed->scalar_type != expected.scalar_type ||
                observed->scalar_class != expected.scalar_class ||
                observed->byte_width != expected.byte_width ||
                observed->is_signed != expected.is_signed ||
                observed->byte_order != expected.byte_order ||
                observed->rank != expected.rank ||
                observed->dimensions != expected.dimensions ||
                observed->disposition != expected.disposition) {
              throw std::runtime_error(
                  "supplied distributed IC manifest field does not match " +
                  expected.dataset_path);
            }
          }
          manifest = supplied;
        }
        validateIcManifest(manifest);
        return assembled;
      });

  std::string manifest_json = mpi_context.isRoot()
      ? serializeIcManifestJson(root_inspection.manifest)
      : std::string{};
  manifest_json = broadcastRootString(mpi_context, std::move(manifest_json));
  checkedCounterAdd(
      local.counters.manifest_metadata_bytes_communicated,
      manifest_json.size(), "manifest_metadata_bytes_communicated");
  Inspection inspection = runCollectivePhase<Inspection>(
      mpi_context, "IC distributed manifest decode", [&]() {
        Inspection decoded;
        decoded.counters = local.counters;
        decoded.manifest = deserializeIcManifestJson(manifest_json);
        populateSchemasFromManifest(decoded);
        return decoded;
      });
  return inspection;
}

[[nodiscard]] int ownerForX(double x, double box_size, int world_size) {
  if (!std::isfinite(x) || !(box_size > 0.0) || x < 0.0 || x > box_size * (1.0 + 1.0e-12)) throw std::runtime_error("converted IC coordinate is outside the periodic box");
  if (x >= box_size) x = 0.0;
  const double fraction=x/box_size;const int owner=std::min(world_size-1,static_cast<int>(fraction*world_size));return std::max(0,owner);
}

[[nodiscard]] std::uint64_t nestedByteCapacity(
    const std::vector<std::vector<std::uint8_t>>& buckets) {
  std::uint64_t bytes =
      static_cast<std::uint64_t>(buckets.capacity()) *
      sizeof(std::vector<std::uint8_t>);
  for (const auto& bucket : buckets) {
    if (bytes > std::numeric_limits<std::uint64_t>::max() -
            bucket.capacity()) {
      throw std::overflow_error("IC nested byte capacity overflow");
    }
    bytes += static_cast<std::uint64_t>(bucket.capacity());
  }
  return bytes;
}

void validateChunkCoverage(
    const parallel::MpiContext& mpi_context,
    std::size_t file_index,
    std::size_t type_index,
    std::size_t start,
    std::size_t count,
    int reader_rank) {
  const int local_completed =
      mpi_context.worldRank() == reader_rank ? 1 : 0;
  int global_completed = 0;
  MPI_Allreduce(
      &local_completed, &global_completed, 1, MPI_INT, MPI_SUM,
      MPI_COMM_WORLD);
  if (global_completed != 1) {
    throw std::runtime_error(
        "distributed IC chunk coverage is not exactly one reader");
  }
  const std::uint64_t token =
      (static_cast<std::uint64_t>(file_index) * 0x9e3779b97f4a7c15ULL) ^
      (static_cast<std::uint64_t>(type_index) * 0xbf58476d1ce4e5b9ULL) ^
      (static_cast<std::uint64_t>(start) * 0x94d049bb133111ebULL) ^
      static_cast<std::uint64_t>(count);
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  MPI_Allreduce(&token, &minimum, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&token, &maximum, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  if (minimum != maximum) {
    throw std::runtime_error(
        "distributed IC ranks disagree on the active chunk token");
  }
}

[[nodiscard]] std::uint64_t exactDistributedChunkReconciliation(
    const parallel::MpiContext& mpi_context,
    std::span<const ParticleRecord> source_records,
    std::span<const ParticleRecord> final_records,
    IcImportCounters& counters) {
  std::vector<std::vector<std::uint8_t>> buckets =
      runCollectivePhase<std::vector<std::vector<std::uint8_t>>>(
          mpi_context, "IC source-to-final audit bucket preparation", [&]() {
            injectIcTestFault(mpi_context, "source_final_reconciliation");
            std::vector<std::vector<std::uint8_t>> prepared(
                static_cast<std::size_t>(mpi_context.worldSize()));
            const auto append_entry = [&](std::uint64_t id,
                                          std::uint32_t source_count,
                                          std::uint32_t final_count) {
              const int validator = static_cast<int>(
                  (id ^ (id >> 33U) ^ (id << 11U)) %
                  static_cast<std::uint64_t>(mpi_context.worldSize()));
              auto& bucket = prepared[static_cast<std::size_t>(validator)];
              appendLe64(bucket, id);
              appendLe32(bucket, source_count);
              appendLe32(bucket, final_count);
            };
            for (const ParticleRecord& record : source_records) {
              append_entry(record.id, 1U, 0U);
            }
            for (const ParticleRecord& record : final_records) {
              append_entry(record.id, 0U, 1U);
            }
            return prepared;
          });

  std::uint64_t sent = 0U;
  std::uint64_t received_bytes = 0U;
  std::uint64_t exchange_peak = 0U;
  const std::vector<std::uint8_t> received = alltoallBytes(
      mpi_context, buckets, sent, received_bytes, &exchange_peak);
  checkedCounterAdd(counters.bytes_sent, sent, "bytes_sent");
  checkedCounterAdd(counters.bytes_received, received_bytes, "bytes_received");

  struct BalanceEntry {
    std::uint64_t id = 0U;
    std::uint32_t source_count = 0U;
    std::uint32_t final_count = 0U;
  };
  std::vector<BalanceEntry> entries =
      runCollectivePhase<std::vector<BalanceEntry>>(
          mpi_context, "IC source-to-final audit decode", [&]() {
            if (received.size() % 16U != 0U) {
              throw std::runtime_error(
                  "IC source-to-final audit wire size is corrupt");
            }
            std::vector<BalanceEntry> decoded;
            decoded.reserve(received.size() / 16U);
            std::size_t offset = 0U;
            while (offset < received.size()) {
              decoded.push_back(BalanceEntry{
                  .id = readLe64(received, offset),
                  .source_count = readLe32(received, offset),
                  .final_count = readLe32(received, offset)});
            }
            return decoded;
          });

  const int local_bad = runCollectivePhase<int>(
      mpi_context, "IC source-to-final audit reduction", [&]() {
        std::sort(
            entries.begin(), entries.end(),
            [](const BalanceEntry& lhs, const BalanceEntry& rhs) {
              return lhs.id < rhs.id;
            });
        std::size_t begin = 0U;
        while (begin < entries.size()) {
          std::size_t end = begin + 1U;
          std::uint64_t source_count = entries[begin].source_count;
          std::uint64_t final_count = entries[begin].final_count;
          while (end < entries.size() &&
                 entries[end].id == entries[begin].id) {
            source_count += entries[end].source_count;
            final_count += entries[end].final_count;
            ++end;
          }
          if (source_count != 1U || final_count != 1U) {
            return 1;
          }
          begin = end;
        }
        return 0;
      });
  int any_bad = 0;
  MPI_Allreduce(&local_bad, &any_bad, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (any_bad != 0) {
    throw std::runtime_error(
        "distributed IC source and final ID multisets do not reconcile "
        "exactly");
  }
  return std::max(
      nestedByteCapacity(buckets) + exchange_peak,
      static_cast<std::uint64_t>(received.capacity()) +
          static_cast<std::uint64_t>(entries.capacity()) *
              sizeof(BalanceEntry));
}

[[nodiscard]] std::uint64_t exactDistributedIdAudit(
    const parallel::MpiContext& mpi_context,
    std::span<const std::uint64_t> local_ids,
    std::size_t batch_count,
    IcImportCounters& counters) {
  const std::size_t validated_batch_count = runCollectivePhase<std::size_t>(
      mpi_context, "IC global ID audit configuration", [&]() {
        if (batch_count == 0U) {
          throw std::invalid_argument(
              "IC ID audit batch count must be positive");
        }
        return batch_count;
      });

  struct TemporaryAuditDirectory {
    std::filesystem::path path;
    ~TemporaryAuditDirectory() {
      std::error_code error;
      std::filesystem::remove_all(path, error);
    }
  };

  const std::filesystem::path audit_directory =
      runCollectivePhase<std::filesystem::path>(
          mpi_context, "IC global ID audit temporary storage", [&]() {
            const auto nonce = std::chrono::high_resolution_clock::now()
                                   .time_since_epoch()
                                   .count();
            std::filesystem::path path =
                std::filesystem::temp_directory_path() /
                ("cosmosim_ic_id_audit_" +
                 std::to_string(mpi_context.worldRank()) + "_" +
                 std::to_string(nonce));
            if (!std::filesystem::create_directories(path)) {
              throw std::runtime_error(
                  "failed to create temporary distributed IC ID-audit "
                  "directory");
            }
            return path;
          });
  TemporaryAuditDirectory audit_cleanup{audit_directory};

  const std::uint64_t local_rounds =
      (local_ids.size() + validated_batch_count - 1U) /
      validated_batch_count;
  std::uint64_t rounds = 0U;
  MPI_Allreduce(
      &local_rounds, &rounds, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);

  std::vector<std::optional<std::filesystem::path>> sorted_run_levels;
  std::uint64_t next_run_index = 0U;
  std::uint64_t peak_bytes = 0U;

  const auto pathStorageBytes = [&]() {
    std::uint64_t bytes =
        static_cast<std::uint64_t>(sorted_run_levels.capacity()) *
        sizeof(std::optional<std::filesystem::path>);
    for (const auto& path : sorted_run_levels) {
      if (path.has_value()) {
        checkedCounterAdd(
            bytes, path->native().capacity() * sizeof(std::filesystem::path::value_type),
            "ID-audit path storage");
      }
    }
    return bytes;
  };

  const auto writeSortedRun = [&](std::span<const std::uint64_t> ids) {
    const std::filesystem::path path =
        audit_directory / ("run_" + std::to_string(next_run_index++) + ".bin");
    std::ofstream output;
    output.rdbuf()->pubsetbuf(nullptr, 0);
    output.open(path, std::ios::binary | std::ios::trunc);
    if (!output) {
      throw std::runtime_error(
          "failed to create distributed IC ID-audit run");
    }
    for (const std::uint64_t id : ids) {
      output.write(
          reinterpret_cast<const char*>(&id),
          static_cast<std::streamsize>(sizeof(id)));
    }
    output.close();
    if (!output) {
      throw std::runtime_error(
          "failed to write distributed IC ID-audit run");
    }
    return path;
  };

  const auto mergeSortedRuns = [&](const std::filesystem::path& lhs_path,
                                   const std::filesystem::path& rhs_path) {
    const std::filesystem::path output_path =
        audit_directory / ("run_" + std::to_string(next_run_index++) + ".bin");
    std::ifstream lhs;
    std::ifstream rhs;
    std::ofstream output;
    lhs.rdbuf()->pubsetbuf(nullptr, 0);
    rhs.rdbuf()->pubsetbuf(nullptr, 0);
    output.rdbuf()->pubsetbuf(nullptr, 0);
    lhs.open(lhs_path, std::ios::binary);
    rhs.open(rhs_path, std::ios::binary);
    output.open(output_path, std::ios::binary | std::ios::trunc);
    if (!lhs || !rhs || !output) {
      throw std::runtime_error(
          "failed to open distributed IC ID-audit merge run");
    }

    std::uint64_t lhs_value = 0U;
    std::uint64_t rhs_value = 0U;
    bool has_lhs = static_cast<bool>(lhs.read(
        reinterpret_cast<char*>(&lhs_value),
        static_cast<std::streamsize>(sizeof(lhs_value))));
    bool has_rhs = static_cast<bool>(rhs.read(
        reinterpret_cast<char*>(&rhs_value),
        static_cast<std::streamsize>(sizeof(rhs_value))));
    std::optional<std::uint64_t> previous;
    bool duplicate = false;
    while (has_lhs || has_rhs) {
      const bool take_lhs = has_lhs && (!has_rhs || lhs_value <= rhs_value);
      const std::uint64_t value = take_lhs ? lhs_value : rhs_value;
      if (previous.has_value() && *previous == value) {
        duplicate = true;
        break;
      }
      previous = value;
      output.write(
          reinterpret_cast<const char*>(&value),
          static_cast<std::streamsize>(sizeof(value)));
      if (!output) {
        throw std::runtime_error(
            "failed to write distributed IC ID-audit merge run");
      }
      if (take_lhs) {
        has_lhs = static_cast<bool>(lhs.read(
            reinterpret_cast<char*>(&lhs_value),
            static_cast<std::streamsize>(sizeof(lhs_value))));
      } else {
        has_rhs = static_cast<bool>(rhs.read(
            reinterpret_cast<char*>(&rhs_value),
            static_cast<std::streamsize>(sizeof(rhs_value))));
      }
    }
    if (!lhs.eof() && lhs.fail()) {
      throw std::runtime_error(
          "failed to read distributed IC ID-audit left run");
    }
    if (!rhs.eof() && rhs.fail()) {
      throw std::runtime_error(
          "failed to read distributed IC ID-audit right run");
    }
    output.close();
    if (!duplicate && !output) {
      throw std::runtime_error(
          "failed to finalize distributed IC ID-audit merge run");
    }
    if (!duplicate) {
      std::error_code error;
      std::filesystem::remove(lhs_path, error);
      error.clear();
      std::filesystem::remove(rhs_path, error);
    }
    return std::pair{output_path, duplicate};
  };

  for (std::uint64_t round = 0U; round < rounds; ++round) {
    std::vector<std::vector<std::uint8_t>> buckets =
        runCollectivePhase<std::vector<std::vector<std::uint8_t>>>(
            mpi_context, "IC global ID audit bucket preparation", [&]() {
              const std::size_t begin = static_cast<std::size_t>(
                  std::min<std::uint64_t>(
                      round * validated_batch_count, local_ids.size()));
              const std::size_t end = std::min(
                  local_ids.size(), begin + validated_batch_count);
              std::vector<std::vector<std::uint8_t>> prepared(
                  static_cast<std::size_t>(mpi_context.worldSize()));
              for (std::size_t index = begin; index < end; ++index) {
                const std::uint64_t id = local_ids[index];
                const int validator = static_cast<int>(
                    (id ^ (id >> 33U) ^ (id << 11U)) %
                    static_cast<std::uint64_t>(mpi_context.worldSize()));
                appendLe64(
                    prepared[static_cast<std::size_t>(validator)], id);
              }
              return prepared;
            });
    std::uint64_t sent = 0U;
    std::uint64_t received_bytes = 0U;
    std::uint64_t exchange_peak = 0U;
    const std::vector<std::uint8_t> received = alltoallBytes(
        mpi_context, buckets, sent, received_bytes, &exchange_peak);
    checkedCounterAdd(counters.bytes_sent, sent, "bytes_sent");
    checkedCounterAdd(
        counters.bytes_received, received_bytes, "bytes_received");

    std::uint64_t round_vector_capacity = 0U;
    const int local_duplicate = runCollectivePhase<int>(
        mpi_context, "IC global ID audit sorted-run update", [&]() {
          if (received.size() % 8U != 0U) {
            throw std::runtime_error(
                "IC ID validation wire size is corrupt");
          }
          std::vector<std::uint64_t> round_ids;
          round_ids.reserve(received.size() / 8U);
          std::size_t offset = 0U;
          while (offset < received.size()) {
            round_ids.push_back(readLe64(received, offset));
          }
          round_vector_capacity =
              static_cast<std::uint64_t>(round_ids.capacity()) *
              sizeof(std::uint64_t);
          std::sort(round_ids.begin(), round_ids.end());
          if (std::adjacent_find(round_ids.begin(), round_ids.end()) !=
              round_ids.end()) {
            return 1;
          }
          if (round_ids.empty()) {
            return 0;
          }

          std::filesystem::path current = writeSortedRun(round_ids);
          std::size_t level = 0U;
          for (;;) {
            if (sorted_run_levels.size() <= level) {
              sorted_run_levels.resize(level + 1U);
            }
            if (!sorted_run_levels[level].has_value()) {
              sorted_run_levels[level] = std::move(current);
              break;
            }
            const auto [merged_path, duplicate] = mergeSortedRuns(
                *sorted_run_levels[level], current);
            sorted_run_levels[level].reset();
            if (duplicate) {
              return 1;
            }
            current = merged_path;
            ++level;
          }
          return 0;
        });
    int any_duplicate = 0;
    MPI_Allreduce(
        &local_duplicate, &any_duplicate, 1, MPI_INT, MPI_MAX,
        MPI_COMM_WORLD);
    if (any_duplicate != 0) {
      throw std::runtime_error(
          "distributed IC import contains duplicate particle IDs across "
          "files or ranks");
    }
    peak_bytes = std::max(
        peak_bytes,
        nestedByteCapacity(buckets) + exchange_peak + round_vector_capacity +
            pathStorageBytes());
  }

  const int local_cross_level_duplicate = runCollectivePhase<int>(
      mpi_context, "IC global ID audit final external merge", [&]() {
        std::optional<std::filesystem::path> merged;
        for (auto& level : sorted_run_levels) {
          if (!level.has_value()) {
            continue;
          }
          if (!merged.has_value()) {
            merged = std::move(*level);
            level.reset();
            continue;
          }
          const auto [merged_path, duplicate] =
              mergeSortedRuns(*merged, *level);
          level.reset();
          if (duplicate) {
            return 1;
          }
          merged = merged_path;
        }
        return 0;
      });
  int any_cross_level_duplicate = 0;
  MPI_Allreduce(
      &local_cross_level_duplicate, &any_cross_level_duplicate, 1, MPI_INT,
      MPI_MAX, MPI_COMM_WORLD);
  if (any_cross_level_duplicate != 0) {
    throw std::runtime_error(
        "distributed IC import contains duplicate particle IDs across files "
        "or ranks");
  }
  return std::max(peak_bytes, pathStorageBytes());
}

void validateDistributedTotals(
    const parallel::MpiContext& mpi_context,
    const core::SimulationState& state,
    const IcManifest& manifest,
    const std::array<double, 5>& local_source_mass) {
  struct LocalTotals {
    std::array<std::uint64_t, 5> counts{};
    std::array<double, 5> masses{};
  };
  const LocalTotals local = runCollectivePhase<LocalTotals>(
      mpi_context, "IC distributed local totals", [&]() {
        LocalTotals totals;
        std::array<double, 5> compensation{};
        for (std::size_t index = 0; index < state.particles.size(); ++index) {
          const std::uint32_t species =
              state.particle_sidecar.species_tag[index];
          if (species >= totals.counts.size()) {
            throw std::runtime_error(
                "invalid species tag after distributed IC routing");
          }
          ++totals.counts[species];
          const double value = state.particles.mass_code[index];
          if (!std::isfinite(value) || !(value > 0.0)) {
            throw std::runtime_error(
                "distributed IC final mass must be finite and positive");
          }
          const double corrected = value - compensation[species];
          const double updated = totals.masses[species] + corrected;
          compensation[species] =
              (updated - totals.masses[species]) - corrected;
          totals.masses[species] = updated;
        }
        return totals;
      });

  std::array<std::uint64_t, 5> global_counts{};
  std::array<double, 5> global_final_mass{};
  std::array<double, 5> global_source_mass{};
  MPI_Allreduce(
      local.counts.data(), global_counts.data(), 5, MPI_UINT64_T, MPI_SUM,
      MPI_COMM_WORLD);
  MPI_Allreduce(
      local.masses.data(), global_final_mass.data(), 5, MPI_DOUBLE, MPI_SUM,
      MPI_COMM_WORLD);
  MPI_Allreduce(
      local_source_mass.data(), global_source_mass.data(), 5, MPI_DOUBLE,
      MPI_SUM, MPI_COMM_WORLD);

  const std::array<std::uint64_t, 5> expected_counts =
      runCollectivePhase<std::array<std::uint64_t, 5>>(
          mpi_context, "IC distributed expected totals", [&]() {
            std::array<std::uint64_t, 5> expected{};
            for (std::size_t type = 0; type < kParticleTypeCount; ++type) {
              if (manifest.num_part_total[type] == 0U) {
                continue;
              }
              const std::uint32_t species =
                  speciesTag(manifest.species_policy[type]);
              if (expected[species] >
                  std::numeric_limits<std::uint64_t>::max() -
                      manifest.num_part_total[type]) {
                throw std::overflow_error(
                    "distributed IC expected species count overflow");
              }
              expected[species] += manifest.num_part_total[type];
            }
            return expected;
          });
  if (global_counts != expected_counts) {
    throw std::runtime_error(
        "distributed IC global species counts do not match the manifest");
  }
  for (std::size_t species = 0; species < expected_counts.size(); ++species) {
    const double tolerance = 1.0e-12 * std::max(
        {1.0, std::abs(global_source_mass[species]),
         std::abs(global_final_mass[species])});
    if (std::abs(
            global_source_mass[species] - global_final_mass[species]) >
        tolerance) {
      throw std::runtime_error(
          "distributed IC global species mass total changed during routing");
    }
  }
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
  if (!mpi_context.isEnabled() || mpi_context.worldSize() == 1) {
    return readGadgetArepoHdf5Ic(ic_path, config, options);
  }
  runCollectivePhaseVoid(
      mpi_context, "IC distributed configuration validation", [&]() {
        if (options.chunk_particle_count == 0U) {
          throw std::invalid_argument("chunk_particle_count must be positive");
        }
        if (config.mode.ic_staging_particle_count == 0U) {
          throw std::invalid_argument(
              "distributed IC staging particle count must be positive");
        }
      });

  Inspection inspection =
      inspectFileSetDistributed(ic_path, config, mpi_context, options);
  if (options.validate_runtime_cosmology) {
    runCollectivePhaseVoid(
        mpi_context, "IC distributed runtime-cosmology validation", [&]() {
          validateRuntimeCosmology(inspection.manifest, config);
        });
  }

  const std::string local_manifest_json =
      runCollectivePhase<std::string>(
          mpi_context, "IC local manifest serialization", [&]() {
            return serializeIcManifestJson(inspection.manifest);
          });
  const std::string local_manifest_digest = sha256Hex(
      std::string_view(local_manifest_json));
  std::string root_manifest_digest =
      mpi_context.isRoot() ? local_manifest_digest : std::string{};
  root_manifest_digest =
      broadcastRootString(mpi_context, std::move(root_manifest_digest));
  const int local_digest_mismatch =
      local_manifest_digest == root_manifest_digest ? 0 : 1;
  int any_digest_mismatch = 0;
  MPI_Allreduce(
      &local_digest_mismatch, &any_digest_mismatch, 1, MPI_INT, MPI_MAX,
      MPI_COMM_WORLD);
  if (any_digest_mismatch != 0) {
    throw std::runtime_error(
        "IC manifest SHA-256 digest is inconsistent across ranks");
  }

  IcReadResult result;
  result.report.counters = inspection.counters;
  result.report.manifest = inspection.manifest;
  result.report.already_partitioned = true;
  result.report.schema.count_by_type =
      inspection.manifest.num_part_this_file.front();
  result.report.schema.total_count_by_type =
      inspection.manifest.num_part_total;
  result.report.schema.total_count_high_word =
      inspection.manifest.num_part_total_high_word;
  result.report.schema.mass_table = inspection.manifest.mass_table;
  result.report.schema.num_files_per_snapshot =
      inspection.manifest.num_files_per_snapshot;
  result.report.schema.box_size = inspection.manifest.box_size;
  result.report.schema.scale_factor = inspection.manifest.scale_factor;
  result.report.schema.redshift = inspection.manifest.redshift;
  result.report.schema.omega_matter = inspection.manifest.omega_matter;
  result.report.schema.omega_lambda = inspection.manifest.omega_lambda;
  result.report.schema.hubble_param = inspection.manifest.hubble_param;
  result.report.defaulted_fields = inspection.manifest.defaulted_fields;
  for (const auto& value : result.report.defaulted_fields) {
    const auto equals = value.find('=');
    result.report.missing_optional_fields.push_back(
        value.substr(0, equals));
  }
  result.report.unsupported_fields = inspection.manifest.dropped_fields;

  struct RoutingConfiguration {
    double box_size = 0.0;
    std::size_t chunk_particle_count = 0U;
  };
  const RoutingConfiguration routing =
      runCollectivePhase<RoutingConfiguration>(
          mpi_context, "IC distributed routing configuration", [&]() {
            RoutingConfiguration values;
            values.box_size = convertedBoxSizeCode(
                inspection.manifest, config);
            values.chunk_particle_count = std::min(
                options.chunk_particle_count,
                config.mode.ic_staging_particle_count);
            if (values.chunk_particle_count == 0U) {
              throw std::invalid_argument(
                  "distributed IC staging particle count must be positive");
            }
            return values;
          });

  std::uint64_t global_chunk_index = 0U;
  std::set<std::size_t> assigned_files;
  std::array<double, 5> local_source_mass{};
  for (std::size_t file_index = 0;
       file_index < inspection.manifest.source_files.size(); ++file_index) {
    for (std::size_t type_index = 0; type_index < kParticleTypeCount;
         ++type_index) {
      const std::size_t total = static_cast<std::size_t>(
          inspection.manifest.num_part_this_file[file_index][type_index]);
      for (std::size_t start = 0; start < total;
           start += routing.chunk_particle_count, ++global_chunk_index) {
        const std::size_t count = std::min(
            routing.chunk_particle_count, total - start);
        const int reader_rank = static_cast<int>(
            global_chunk_index %
            static_cast<std::uint64_t>(mpi_context.worldSize()));

        std::vector<ParticleRecord> records =
            runCollectivePhase<std::vector<ParticleRecord>>(
                mpi_context, "IC chunk read and scientific conversion", [&]() {
                  if (mpi_context.worldRank() != reader_rank) {
                    return std::vector<ParticleRecord>{};
                  }
                  std::vector<ParticleRecord> local_records = readRecordChunk(
                      inspection, file_index, type_index, start, count, config,
                      options, result.report.counters);
                  assigned_files.insert(file_index);
                  ++result.report.counters.chunks_assigned;
                  for (const ParticleRecord& record : local_records) {
                    if (record.species >= local_source_mass.size()) {
                      throw std::runtime_error(
                          "source IC record has invalid species tag");
                    }
                    local_source_mass[record.species] += record.mass;
                  }
                  return local_records;
                });
        validateChunkCoverage(
            mpi_context, file_index, type_index, start, count, reader_rank);

        std::vector<std::vector<std::uint8_t>> per_rank =
            runCollectivePhase<std::vector<std::vector<std::uint8_t>>>(
                mpi_context, "IC owner calculation and serialization", [&]() {
                  injectIcTestFault(mpi_context, "owner_serialization");
                  std::vector<std::vector<std::uint8_t>> buckets(
                      static_cast<std::size_t>(mpi_context.worldSize()));
                  if (mpi_context.worldRank() == reader_rank) {
                    for (const ParticleRecord& record : records) {
                      const int owner = ownerForX(
                          record.x, routing.box_size,
                          mpi_context.worldSize());
                      serializeRecord(
                          record, buckets[static_cast<std::size_t>(owner)]);
                    }
                    checkedCounterAdd(
                        result.report.counters.records_routed,
                        records.size(), "records_routed");
                    mutateIcTestRoute(mpi_context, buckets);
                  }
                  return buckets;
                });
        const std::uint64_t serialized_capacity =
            runCollectivePhase<std::uint64_t>(
                mpi_context, "IC serialized-capacity accounting", [&]() {
                  return nestedByteCapacity(per_rank);
                });
        std::uint64_t serialized_size = 0U;
        for (const auto& bucket : per_rank) {
          checkedCounterAdd(
              serialized_size, bucket.size(), "serialized payload size");
        }
        checkedCounterAdd(
            result.report.counters.bytes_serialized, serialized_size,
            "bytes_serialized");

        runCollectivePhaseVoid(
            mpi_context, "IC send-layout fault-injection boundary", [&]() {
              injectIcTestFault(mpi_context, "send_layout");
            });
        std::uint64_t sent = 0U;
        std::uint64_t received_bytes = 0U;
        std::uint64_t exchange_peak = 0U;
        const std::vector<std::uint8_t> inbound_bytes = alltoallBytes(
            mpi_context, per_rank, sent, received_bytes, &exchange_peak);
        checkedCounterAdd(
            result.report.counters.bytes_sent, sent, "bytes_sent");
        checkedCounterAdd(
            result.report.counters.bytes_received, received_bytes,
            "bytes_received");

        std::vector<ParticleRecord> inbound =
            runCollectivePhase<std::vector<ParticleRecord>>(
                mpi_context, "IC wire validation and deserialization", [&]() {
                  injectIcTestFault(mpi_context, "payload_validation");
                  if (inbound_bytes.size() % kWireRecordBytes != 0U) {
                    throw std::runtime_error(
                        "distributed IC wire payload has invalid length");
                  }
                  std::vector<ParticleRecord> decoded;
                  decoded.reserve(inbound_bytes.size() / kWireRecordBytes);
                  std::size_t offset = 0U;
                  while (offset < inbound_bytes.size()) {
                    injectIcTestFault(mpi_context, "deserialization");
                    ParticleRecord record =
                        deserializeRecord(inbound_bytes, offset);
                    if (record.species >= 5U) {
                      throw std::runtime_error(
                          "distributed IC wire record has invalid species");
                    }
                    IcSpeciesPolicy policy = IcSpeciesPolicy::kDarkMatter;
                    if (record.species == static_cast<std::uint32_t>(
                                              core::ParticleSpecies::kGas)) {
                      policy = IcSpeciesPolicy::kGas;
                    } else if (record.species == static_cast<std::uint32_t>(
                                                     core::ParticleSpecies::kStar)) {
                      policy = IcSpeciesPolicy::kStar;
                    } else if (record.species == static_cast<std::uint32_t>(
                                                     core::ParticleSpecies::kBlackHole)) {
                      policy = IcSpeciesPolicy::kBlackHole;
                    } else if (record.species == static_cast<std::uint32_t>(
                                                     core::ParticleSpecies::kTracer)) {
                      policy = IcSpeciesPolicy::kTracer;
                    }
                    validateRecordScientificState(
                        record, policy, routing.box_size);
                    if (ownerForX(
                            record.x, routing.box_size,
                            mpi_context.worldSize()) !=
                        mpi_context.worldRank()) {
                      throw std::runtime_error(
                          "distributed IC record arrived at the wrong owner");
                    }
                    decoded.push_back(record);
                  }
                  return decoded;
                });

        const std::uint64_t reconciliation_peak =
            exactDistributedChunkReconciliation(
                mpi_context, records, inbound, result.report.counters);
        runCollectivePhaseVoid(
            mpi_context, "IC authoritative-state append", [&]() {
              injectIcTestFault(mpi_context, "sidecar_append");
              appendRecords(
                  result.state, inbound,
                  static_cast<std::uint32_t>(mpi_context.worldRank()));
            });

        const std::uint64_t routing_peak =
            vectorCapacityBytes(records) + serialized_capacity + exchange_peak +
            static_cast<std::uint64_t>(inbound_bytes.capacity()) +
            vectorCapacityBytes(inbound) + reconciliation_peak;
        result.report.counters.peak_staging_bytes = std::max(
            result.report.counters.peak_staging_bytes, routing_peak);
      }
    }
  }

  result.report.counters.files_assigned = assigned_files.size();
  runCollectivePhaseVoid(
      mpi_context, "IC final authoritative-state construction", [&]() {
        injectIcTestFault(mpi_context, "final_state");
        finalizeImportedState(result.state, inspection.manifest, config);
      });
  result.report.counters.peak_staging_bytes = std::max(
      result.report.counters.peak_staging_bytes,
      exactDistributedIdAudit(
          mpi_context, result.state.particle_sidecar.particle_id,
          config.mode.ic_staging_particle_count,
          result.report.counters));
  validateDistributedTotals(
      mpi_context, result.state, inspection.manifest, local_source_mass);

  result.report.counters.final_local_particle_count =
      result.state.particles.size();
  result.report.counters.final_local_gas_cell_count = result.state.cells.size();
  result.report.counters.final_local_star_count =
      result.state.star_particles.size();
  result.report.counters.final_local_black_hole_count =
      result.state.black_holes.size();
  result.report.counters.final_local_tracer_count = result.state.tracers.size();
  result.report.counters.bytes_read =
      result.report.counters.metadata_bytes_read +
      result.report.counters.hash_bytes_read +
      result.report.counters.payload_bytes_read;
  return result;
#endif
}

}  // namespace cosmosim::io
