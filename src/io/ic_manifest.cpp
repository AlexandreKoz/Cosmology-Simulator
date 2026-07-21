#include "cosmosim/io/ic_reader.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <variant>

namespace cosmosim::io {
namespace {

[[nodiscard]] std::string escapeJson(std::string_view value) {
  std::string result;
  result.reserve(value.size());
  for (const char ch : value) {
    switch (ch) {
      case '\\': result += "\\\\"; break;
      case '"': result += "\\\""; break;
      case '\n': result += "\\n"; break;
      case '\r': result += "\\r"; break;
      case '\t': result += "\\t"; break;
      default: result.push_back(ch); break;
    }
  }
  return result;
}

[[nodiscard]] std::string_view dialectName(IcDialect value) {
  switch (value) {
    case IcDialect::kChuiCanonicalV1: return "chui_canonical_v1";
    case IcDialect::kGadgetArepoBridgeV1: return "gadget_arepo_bridge_v1";
  }
  throw std::invalid_argument("unknown IC dialect");
}
[[nodiscard]] IcDialect parseDialect(std::string_view value) {
  if (value == "chui_canonical_v1") return IcDialect::kChuiCanonicalV1;
  if (value == "gadget_arepo_bridge_v1") return IcDialect::kGadgetArepoBridgeV1;
  throw std::invalid_argument("unsupported IC dialect: " + std::string(value));
}
[[nodiscard]] std::string_view frameName(IcCoordinateFrame value) {
  return value == IcCoordinateFrame::kComoving ? "comoving" : "physical";
}
[[nodiscard]] IcCoordinateFrame parseFrame(std::string_view value) {
  if (value == "comoving") return IcCoordinateFrame::kComoving;
  if (value == "physical") return IcCoordinateFrame::kPhysical;
  throw std::invalid_argument("unsupported IC coordinate frame: " + std::string(value));
}
[[nodiscard]] std::string_view velocityConventionName(IcVelocityConvention value) {
  switch (value) {
    case IcVelocityConvention::kNotVelocity: return "not_velocity";
    case IcVelocityConvention::kPhysicalPeculiar: return "physical_peculiar";
    case IcVelocityConvention::kSqrtAScaledPeculiar: return "sqrt_a_scaled_peculiar";
    case IcVelocityConvention::kComovingCoordinateRate: return "comoving_coordinate_rate";
  }
  throw std::invalid_argument("unknown IC velocity convention");
}
[[nodiscard]] IcVelocityConvention parseVelocityConvention(std::string_view value) {
  if (value == "not_velocity") return IcVelocityConvention::kNotVelocity;
  if (value == "physical_peculiar") return IcVelocityConvention::kPhysicalPeculiar;
  if (value == "sqrt_a_scaled_peculiar") return IcVelocityConvention::kSqrtAScaledPeculiar;
  if (value == "comoving_coordinate_rate") return IcVelocityConvention::kComovingCoordinateRate;
  throw std::invalid_argument("unsupported IC velocity convention: " + std::string(value));
}
[[nodiscard]] std::string_view semanticsName(IcFieldSemantics value) {
  switch (value) {
    case IcFieldSemantics::kCoordinate: return "coordinate";
    case IcFieldSemantics::kVelocity: return "velocity";
    case IcFieldSemantics::kExtensive: return "extensive";
    case IcFieldSemantics::kIntensive: return "intensive";
    case IcFieldSemantics::kSpecific: return "specific";
    case IcFieldSemantics::kIdentifier: return "identifier";
  }
  throw std::invalid_argument("unknown IC field semantics");
}
[[nodiscard]] IcFieldSemantics parseSemantics(std::string_view value) {
  if (value == "coordinate") return IcFieldSemantics::kCoordinate;
  if (value == "velocity") return IcFieldSemantics::kVelocity;
  if (value == "extensive") return IcFieldSemantics::kExtensive;
  if (value == "intensive") return IcFieldSemantics::kIntensive;
  if (value == "specific") return IcFieldSemantics::kSpecific;
  if (value == "identifier") return IcFieldSemantics::kIdentifier;
  throw std::invalid_argument("unsupported IC field semantics: " + std::string(value));
}
[[nodiscard]] std::string_view dispositionName(IcFieldDisposition value) {
  switch (value) {
    case IcFieldDisposition::kConverted: return "converted";
    case IcFieldDisposition::kPreserved: return "preserved";
    case IcFieldDisposition::kDefaulted: return "defaulted";
    case IcFieldDisposition::kDropped: return "dropped";
    case IcFieldDisposition::kRejected: return "rejected";
  }
  throw std::invalid_argument("unknown IC field disposition");
}
[[nodiscard]] IcFieldDisposition parseDisposition(std::string_view value) {
  if (value == "converted") return IcFieldDisposition::kConverted;
  if (value == "preserved") return IcFieldDisposition::kPreserved;
  if (value == "defaulted") return IcFieldDisposition::kDefaulted;
  if (value == "dropped") return IcFieldDisposition::kDropped;
  if (value == "rejected") return IcFieldDisposition::kRejected;
  throw std::invalid_argument("unsupported IC field disposition: " + std::string(value));
}
[[nodiscard]] std::string_view speciesPolicyName(IcSpeciesPolicy value) {
  switch (value) {
    case IcSpeciesPolicy::kReject: return "reject";
    case IcSpeciesPolicy::kGas: return "gas";
    case IcSpeciesPolicy::kDarkMatter: return "dark_matter";
    case IcSpeciesPolicy::kCollisionlessFamily2AsDarkMatter: return "family2_as_dark_matter";
    case IcSpeciesPolicy::kCollisionlessFamily3AsDarkMatter: return "family3_as_dark_matter";
    case IcSpeciesPolicy::kStar: return "star";
    case IcSpeciesPolicy::kBlackHole: return "black_hole";
    case IcSpeciesPolicy::kTracer: return "tracer";
  }
  throw std::invalid_argument("unknown IC species policy");
}
[[nodiscard]] IcSpeciesPolicy parseSpeciesPolicy(std::string_view value) {
  if (value == "reject") return IcSpeciesPolicy::kReject;
  if (value == "gas") return IcSpeciesPolicy::kGas;
  if (value == "dark_matter") return IcSpeciesPolicy::kDarkMatter;
  if (value == "family2_as_dark_matter") return IcSpeciesPolicy::kCollisionlessFamily2AsDarkMatter;
  if (value == "family3_as_dark_matter") return IcSpeciesPolicy::kCollisionlessFamily3AsDarkMatter;
  if (value == "star") return IcSpeciesPolicy::kStar;
  if (value == "black_hole") return IcSpeciesPolicy::kBlackHole;
  if (value == "tracer") return IcSpeciesPolicy::kTracer;
  throw std::invalid_argument("unsupported IC species policy: " + std::string(value));
}
[[nodiscard]] std::string_view scalarClassName(IcScalarClass value) {
  return value == IcScalarClass::kFloatingPoint ? "floating_point" : "integer";
}
[[nodiscard]] IcScalarClass parseScalarClass(std::string_view value) {
  if (value == "floating_point") return IcScalarClass::kFloatingPoint;
  if (value == "integer") return IcScalarClass::kInteger;
  throw std::invalid_argument("unsupported IC scalar class: " + std::string(value));
}
[[nodiscard]] std::string_view byteOrderName(IcByteOrder value) {
  switch (value) {
    case IcByteOrder::kNotApplicable: return "not_applicable";
    case IcByteOrder::kLittleEndian: return "little_endian";
    case IcByteOrder::kBigEndian: return "big_endian";
    case IcByteOrder::kNative: return "native";
  }
  throw std::invalid_argument("unknown IC byte order");
}
[[nodiscard]] IcByteOrder parseByteOrder(std::string_view value) {
  if (value == "not_applicable") return IcByteOrder::kNotApplicable;
  if (value == "little_endian") return IcByteOrder::kLittleEndian;
  if (value == "big_endian") return IcByteOrder::kBigEndian;
  if (value == "native") return IcByteOrder::kNative;
  throw std::invalid_argument("unsupported IC byte order: " + std::string(value));
}

struct JsonValue {
  enum class Kind { kNull, kBool, kNumber, kString, kArray, kObject };
  Kind kind = Kind::kNull;
  bool bool_value = false;
  std::string text;
  std::vector<JsonValue> array;
  std::map<std::string, JsonValue> object;
};

class JsonParser {
 public:
  explicit JsonParser(std::string_view input) : m_input(input) {}
  [[nodiscard]] JsonValue parse() {
    skipSpace();
    JsonValue value = parseValue();
    skipSpace();
    if (m_pos != m_input.size()) fail("trailing content");
    return value;
  }
 private:
  [[noreturn]] void fail(std::string_view message) const {
    throw std::invalid_argument("invalid IC manifest JSON at byte " + std::to_string(m_pos) + ": " + std::string(message));
  }
  void skipSpace() {
    while (m_pos < m_input.size() && (m_input[m_pos] == ' ' || m_input[m_pos] == '\n' || m_input[m_pos] == '\r' || m_input[m_pos] == '\t')) ++m_pos;
  }
  [[nodiscard]] bool consume(char expected) {
    skipSpace();
    if (m_pos < m_input.size() && m_input[m_pos] == expected) { ++m_pos; return true; }
    return false;
  }
  void require(char expected) { if (!consume(expected)) fail(std::string("expected '") + expected + "'"); }
  [[nodiscard]] JsonValue parseValue() {
    skipSpace();
    if (m_pos >= m_input.size()) fail("unexpected end of input");
    const char ch = m_input[m_pos];
    if (ch == '{') return parseObject();
    if (ch == '[') return parseArray();
    if (ch == '"') { JsonValue v; v.kind=JsonValue::Kind::kString; v.text=parseString(); return v; }
    if (ch == '-' || (ch >= '0' && ch <= '9')) return parseNumber();
    if (m_input.substr(m_pos, 4) == "true") { m_pos += 4; JsonValue v; v.kind=JsonValue::Kind::kBool; v.bool_value=true; return v; }
    if (m_input.substr(m_pos, 5) == "false") { m_pos += 5; JsonValue v; v.kind=JsonValue::Kind::kBool; return v; }
    if (m_input.substr(m_pos, 4) == "null") { m_pos += 4; return {}; }
    fail("unexpected token");
  }
  [[nodiscard]] JsonValue parseObject() {
    require('{'); JsonValue v; v.kind=JsonValue::Kind::kObject;
    skipSpace(); if (consume('}')) return v;
    while (true) {
      skipSpace(); if (m_pos >= m_input.size() || m_input[m_pos] != '"') fail("object key must be a string");
      std::string key=parseString(); require(':');
      if (!v.object.emplace(std::move(key), parseValue()).second) fail("duplicate object key");
      if (consume('}')) break; require(',');
    }
    return v;
  }
  [[nodiscard]] JsonValue parseArray() {
    require('['); JsonValue v; v.kind=JsonValue::Kind::kArray;
    skipSpace(); if (consume(']')) return v;
    while (true) { v.array.push_back(parseValue()); if (consume(']')) break; require(','); }
    return v;
  }
  [[nodiscard]] std::string parseString() {
    require('"'); std::string out;
    while (m_pos < m_input.size()) {
      char ch=m_input[m_pos++];
      if (ch=='"') return out;
      if (ch=='\\') {
        if (m_pos >= m_input.size()) fail("truncated escape");
        ch=m_input[m_pos++];
        switch (ch) {
          case '"': out.push_back('"'); break; case '\\': out.push_back('\\'); break;
          case '/': out.push_back('/'); break; case 'b': out.push_back('\b'); break;
          case 'f': out.push_back('\f'); break; case 'n': out.push_back('\n'); break;
          case 'r': out.push_back('\r'); break; case 't': out.push_back('\t'); break;
          default: fail("unsupported string escape");
        }
      } else { out.push_back(ch); }
    }
    fail("unterminated string");
  }
  [[nodiscard]] JsonValue parseNumber() {
    const std::size_t begin=m_pos;
    if (m_input[m_pos]=='-') ++m_pos;
    if (m_pos >= m_input.size() || m_input[m_pos]<'0' || m_input[m_pos]>'9') fail("invalid number");
    if (m_input[m_pos]=='0') ++m_pos; else while (m_pos<m_input.size() && m_input[m_pos]>='0' && m_input[m_pos]<='9') ++m_pos;
    if (m_pos<m_input.size() && m_input[m_pos]=='.') { ++m_pos; if (m_pos>=m_input.size() || m_input[m_pos]<'0' || m_input[m_pos]>'9') fail("invalid fraction"); while (m_pos<m_input.size() && m_input[m_pos]>='0' && m_input[m_pos]<='9') ++m_pos; }
    if (m_pos<m_input.size() && (m_input[m_pos]=='e' || m_input[m_pos]=='E')) { ++m_pos; if (m_pos<m_input.size() && (m_input[m_pos]=='+' || m_input[m_pos]=='-')) ++m_pos; if (m_pos>=m_input.size() || m_input[m_pos]<'0' || m_input[m_pos]>'9') fail("invalid exponent"); while (m_pos<m_input.size() && m_input[m_pos]>='0' && m_input[m_pos]<='9') ++m_pos; }
    JsonValue v; v.kind=JsonValue::Kind::kNumber; v.text=std::string(m_input.substr(begin,m_pos-begin)); return v;
  }
  std::string_view m_input; std::size_t m_pos=0;
};

[[nodiscard]] const JsonValue& member(const JsonValue& object, std::string_view name) {
  if (object.kind != JsonValue::Kind::kObject) throw std::invalid_argument("IC manifest root/member must be an object");
  const auto it=object.object.find(std::string(name));
  if (it==object.object.end()) throw std::invalid_argument("IC manifest missing required member: " + std::string(name));
  return it->second;
}
[[nodiscard]] std::string asString(const JsonValue& value, std::string_view name) {
  if (value.kind != JsonValue::Kind::kString) throw std::invalid_argument("IC manifest member must be string: " + std::string(name));
  return value.text;
}
[[nodiscard]] bool asBool(const JsonValue& value, std::string_view name) {
  if (value.kind != JsonValue::Kind::kBool) throw std::invalid_argument("IC manifest member must be bool: " + std::string(name));
  return value.bool_value;
}
[[nodiscard]] std::uint64_t asU64(const JsonValue& value, std::string_view name) {
  if (value.kind != JsonValue::Kind::kNumber || value.text.find_first_of("-.eE") != std::string::npos) throw std::invalid_argument("IC manifest member must be unsigned integer: " + std::string(name));
  std::uint64_t out=0; const auto result=std::from_chars(value.text.data(),value.text.data()+value.text.size(),out);
  if (result.ec!=std::errc{} || result.ptr!=value.text.data()+value.text.size()) throw std::invalid_argument("IC manifest unsigned integer out of range: " + std::string(name));
  return out;
}
[[nodiscard]] double asDouble(const JsonValue& value, std::string_view name) {
  if (value.kind != JsonValue::Kind::kNumber) throw std::invalid_argument("IC manifest member must be number: " + std::string(name));
  std::size_t used=0; double out=0.0; try { out=std::stod(value.text,&used); } catch (...) { throw std::invalid_argument("invalid IC manifest number: " + std::string(name)); }
  if (used!=value.text.size() || !std::isfinite(out)) throw std::invalid_argument("invalid IC manifest finite number: " + std::string(name));
  return out;
}
[[nodiscard]] const std::vector<JsonValue>& asArray(const JsonValue& value, std::string_view name) {
  if (value.kind != JsonValue::Kind::kArray) throw std::invalid_argument("IC manifest member must be array: " + std::string(name));
  return value.array;
}
[[nodiscard]] std::vector<std::string> stringArray(const JsonValue& value, std::string_view name) {
  std::vector<std::string> out; for (const auto& item:asArray(value,name)) out.push_back(asString(item,name)); return out;
}

template <typename T>
void writeNumericArray(std::ostringstream& out, const T& values) {
  out << '['; for (std::size_t i=0;i<values.size();++i) out << (i==0?"":", ") << values[i]; out << ']';
}
void writeStringArray(std::ostringstream& out, const std::vector<std::string>& values) {
  out << '['; for (std::size_t i=0;i<values.size();++i) out << (i==0?"":", ") << '"' << escapeJson(values[i]) << '"'; out << ']';
}

}  // namespace

void validateIcManifest(const IcManifest& manifest) {
  if (manifest.schema_name != "chui_ic_audit_manifest" || manifest.schema_version != 2U) throw std::invalid_argument("IcManifest requires schema chui_ic_audit_manifest version 2");
  if (manifest.converter_version.empty() || manifest.dialect_version.empty()) throw std::invalid_argument("IcManifest converter/dialect versions are required");
  if (manifest.source_files.empty()) throw std::invalid_argument("IcManifest requires at least one source file");
  const std::size_t file_count=manifest.source_files.size();
  if (manifest.num_files_per_snapshot!=file_count || manifest.num_part_this_file.size()!=file_count || manifest.source_provenance_ids.size()!=file_count || manifest.source_file_sizes_bytes.size()!=file_count || manifest.source_sha256.size()!=file_count || manifest.original_header_attributes.size()!=file_count) throw std::invalid_argument("IcManifest requires one count/hash/size/header record per source file");
  std::unordered_set<std::string> paths;
  for (std::size_t i=0;i<file_count;++i) {
    const std::string path=manifest.source_files[i].lexically_normal().generic_string();
    if (path.empty() || !paths.insert(path).second) throw std::invalid_argument("IcManifest source files must be non-empty and unique");
    if (manifest.source_file_sizes_bytes[i]==0U || manifest.source_sha256[i].size()!=64U || manifest.source_provenance_ids[i] != "sha256:" + manifest.source_sha256[i] || manifest.original_header_attributes[i].empty()) throw std::invalid_argument("IcManifest source hashes, sizes, provenance IDs, and original headers must be complete");
  }
  if (!(manifest.scale_factor>0.0) || !std::isfinite(manifest.scale_factor) || !std::isfinite(manifest.redshift) || !(manifest.hubble_param>0.0) || !std::isfinite(manifest.hubble_param) || !(manifest.box_size>0.0) || !std::isfinite(manifest.box_size) || !std::isfinite(manifest.omega_matter) || !std::isfinite(manifest.omega_lambda)) throw std::invalid_argument("IcManifest cosmology, BoxSize, h, and epoch must be finite and physical");
  const double expected_z=1.0/manifest.scale_factor-1.0;
  if (std::abs(expected_z-manifest.redshift)>1.0e-10*std::max(1.0,std::abs(expected_z))) throw std::invalid_argument("IcManifest Time/Redshift values are inconsistent");
  for (std::size_t type=0;type<6;++type) {
    if ((manifest.num_part_total[type]>>32U)!=manifest.num_part_total_high_word[type]) throw std::invalid_argument("IcManifest high word disagrees with 64-bit total");
    std::uint64_t sum=0; for (const auto& per_file:manifest.num_part_this_file) { if (sum>std::numeric_limits<std::uint64_t>::max()-per_file[type]) throw std::overflow_error("IcManifest particle count overflow"); sum+=per_file[type]; }
    if (sum!=manifest.num_part_total[type]) throw std::invalid_argument("IcManifest per-file counts do not sum to total");
    if (sum>0U && manifest.species_policy[type]==IcSpeciesPolicy::kReject) throw std::invalid_argument("IcManifest populated particle type has reject species policy");
    if (!std::isfinite(manifest.mass_table[type]) || manifest.mass_table[type]<0.0) throw std::invalid_argument("IcManifest MassTable values must be finite and non-negative");
  }
  std::unordered_set<std::string> field_keys;
  for (const IcFieldManifest& field:manifest.fields) {
    const std::string key=std::to_string(field.source_file_index)+":"+field.dataset_path;
    if (field.source_file_index>=file_count || field.dataset_path.empty() || field.selected_alias.empty() || field.scalar_type.empty() || !field_keys.insert(key).second) throw std::invalid_argument("IcManifest field source/path/type/alias must be complete and unique per file");
    const bool scalar_shape = field.rank==0U && field.dimensions.empty() && field.record_count==1U;
    const bool array_shape = field.rank>0U && field.rank==field.dimensions.size() && !field.dimensions.empty() && field.dimensions.front()==field.record_count;
    if (field.byte_width==0U || (!scalar_shape && !array_shape)) throw std::invalid_argument("IcManifest field datatype/rank/dimensions/count are inconsistent");
    if (!(field.base_unit_to_si>0.0) || !std::isfinite(field.base_unit_to_si) || !std::isfinite(field.hubble_exponent) || !std::isfinite(field.scale_factor_exponent) || field.source_unit.empty() || field.target_unit.empty() || field.conversion_equation.empty()) throw std::invalid_argument("IcManifest field conversion contract must be complete and finite");
    if (field.semantics==IcFieldSemantics::kVelocity && field.velocity_convention==IcVelocityConvention::kNotVelocity) throw std::invalid_argument("IcManifest velocity field requires a velocity convention");
  }
}

double icStoredToSiMultiplier(const IcFieldManifest& field,double hubble_param,double scale_factor) {
  if (!(field.base_unit_to_si>0.0) || !(hubble_param>0.0) || !(scale_factor>0.0)) throw std::invalid_argument("IC unit conversion requires positive base unit, h, and scale factor");
  return field.base_unit_to_si*std::pow(hubble_param,field.hubble_exponent)*std::pow(scale_factor,field.scale_factor_exponent);
}
double icVelocityConventionMultiplier(IcVelocityConvention convention,double scale_factor) {
  if (!(scale_factor>0.0)) throw std::invalid_argument("IC velocity conversion requires scale_factor > 0");
  switch (convention) {
    case IcVelocityConvention::kNotVelocity:
    case IcVelocityConvention::kPhysicalPeculiar: return 1.0;
    case IcVelocityConvention::kSqrtAScaledPeculiar: return 1.0/std::sqrt(scale_factor);
    case IcVelocityConvention::kComovingCoordinateRate: return scale_factor;
  }
  throw std::invalid_argument("unknown IC velocity convention");
}

std::string serializeIcManifestJson(const IcManifest& manifest) {
  validateIcManifest(manifest);
  std::ostringstream out; out << std::setprecision(std::numeric_limits<double>::max_digits10);
  out << "{\n  \"schema_name\": \"" << escapeJson(manifest.schema_name) << "\",\n  \"schema_version\": " << manifest.schema_version << ",\n  \"converter_version\": \"" << escapeJson(manifest.converter_version) << "\",\n  \"dialect\": \"" << dialectName(manifest.dialect) << "\",\n  \"dialect_version\": \"" << escapeJson(manifest.dialect_version) << "\",\n  \"num_files_per_snapshot\": " << manifest.num_files_per_snapshot << ",\n";
  out << "  \"source_files\": ["; for (std::size_t i=0;i<manifest.source_files.size();++i) out << (i?", ":"") << '"' << escapeJson(manifest.source_files[i].generic_string()) << '"'; out << "],\n  \"source_provenance_ids\": "; writeStringArray(out,manifest.source_provenance_ids);
  out << ",\n  \"source_file_sizes_bytes\": "; writeNumericArray(out,manifest.source_file_sizes_bytes);
  out << ",\n  \"source_sha256\": "; writeStringArray(out,manifest.source_sha256);
  out << ",\n  \"original_header_attributes\": "; writeStringArray(out,manifest.original_header_attributes);
  out << ",\n  \"num_part_this_file\": ["; for (std::size_t i=0;i<manifest.num_part_this_file.size();++i) { if(i) out<<", "; writeNumericArray(out,manifest.num_part_this_file[i]); } out << "],\n  \"num_part_total\": "; writeNumericArray(out,manifest.num_part_total);
  out << ",\n  \"num_part_total_high_word\": "; writeNumericArray(out,manifest.num_part_total_high_word);
  out << ",\n  \"mass_table\": "; writeNumericArray(out,manifest.mass_table);
  out << ",\n  \"box_size\": " << manifest.box_size << ",\n  \"time\": " << manifest.scale_factor << ",\n  \"redshift\": " << manifest.redshift << ",\n  \"omega_matter\": " << manifest.omega_matter << ",\n  \"omega_lambda\": " << manifest.omega_lambda << ",\n  \"hubble_param\": " << manifest.hubble_param << ",\n  \"species_policy\": [";
  for (std::size_t i=0;i<6;++i) out << (i?", ":"") << '"' << speciesPolicyName(manifest.species_policy[i]) << '"';
  out << "],\n  \"fields\": [";
  for (std::size_t i=0;i<manifest.fields.size();++i) {
    const auto& f=manifest.fields[i]; out << (i?",\n":"\n") << "    {\"source_file_index\": " << f.source_file_index << ", \"dataset_path\": \"" << escapeJson(f.dataset_path) << "\", \"selected_alias\": \"" << escapeJson(f.selected_alias) << "\", \"scalar_type\": \"" << escapeJson(f.scalar_type) << "\", \"scalar_class\": \"" << scalarClassName(f.scalar_class) << "\", \"byte_width\": " << static_cast<unsigned>(f.byte_width) << ", \"is_signed\": " << (f.is_signed?"true":"false") << ", \"byte_order\": \"" << byteOrderName(f.byte_order) << "\", \"rank\": " << static_cast<unsigned>(f.rank) << ", \"dimensions\": "; writeNumericArray(out,f.dimensions);
    out << ", \"record_count\": " << f.record_count << ", \"base_unit_to_si\": " << f.base_unit_to_si << ", \"hubble_exponent\": " << f.hubble_exponent << ", \"scale_factor_exponent\": " << f.scale_factor_exponent << ", \"coordinate_frame\": \"" << frameName(f.coordinate_frame) << "\", \"velocity_convention\": \"" << velocityConventionName(f.velocity_convention) << "\", \"semantics\": \"" << semanticsName(f.semantics) << "\", \"disposition\": \"" << dispositionName(f.disposition) << "\", \"source_unit\": \"" << escapeJson(f.source_unit) << "\", \"target_unit\": \"" << escapeJson(f.target_unit) << "\", \"conversion_equation\": \"" << escapeJson(f.conversion_equation) << "\"}";
  }
  out << (manifest.fields.empty()?"":"\n") << "  ],\n";
  const auto write_named_strings=[&](std::string_view name,const std::vector<std::string>& values,bool final=false){ out << "  \"" << name << "\": "; writeStringArray(out,values); out << (final?"\n":",\n"); };
  write_named_strings("defaulted_fields",manifest.defaulted_fields); write_named_strings("converted_fields",manifest.converted_fields); write_named_strings("dropped_fields",manifest.dropped_fields); write_named_strings("rejected_fields",manifest.rejected_fields); write_named_strings("preserved_auxiliary_fields",manifest.preserved_auxiliary_fields); write_named_strings("conversion_equations",manifest.conversion_equations); write_named_strings("warnings",manifest.warnings,true);
  out << "}\n"; return out.str();
}

IcManifest deserializeIcManifestJson(std::string_view json_text) {
  const JsonValue root=JsonParser(json_text).parse(); IcManifest m;
  m.schema_name=asString(member(root,"schema_name"),"schema_name"); m.schema_version=static_cast<std::uint32_t>(asU64(member(root,"schema_version"),"schema_version")); m.converter_version=asString(member(root,"converter_version"),"converter_version"); m.dialect=parseDialect(asString(member(root,"dialect"),"dialect")); m.dialect_version=asString(member(root,"dialect_version"),"dialect_version"); m.num_files_per_snapshot=static_cast<std::uint32_t>(asU64(member(root,"num_files_per_snapshot"),"num_files_per_snapshot"));
  for (const auto& v:asArray(member(root,"source_files"),"source_files")) m.source_files.emplace_back(asString(v,"source_files"));
  m.source_provenance_ids=stringArray(member(root,"source_provenance_ids"),"source_provenance_ids"); for (const auto& v:asArray(member(root,"source_file_sizes_bytes"),"source_file_sizes_bytes")) m.source_file_sizes_bytes.push_back(asU64(v,"source_file_sizes_bytes")); m.source_sha256=stringArray(member(root,"source_sha256"),"source_sha256"); m.original_header_attributes=stringArray(member(root,"original_header_attributes"),"original_header_attributes");
  for (const auto& row:asArray(member(root,"num_part_this_file"),"num_part_this_file")) { const auto& a=asArray(row,"num_part_this_file row"); if(a.size()!=6U) throw std::invalid_argument("num_part_this_file rows require six values"); std::array<std::uint64_t,6> values{}; for(std::size_t i=0;i<6;++i) values[i]=asU64(a[i],"num_part_this_file"); m.num_part_this_file.push_back(values); }
  const auto fill_u64_6=[&](std::string_view name,auto& target){ const auto& a=asArray(member(root,name),name); if(a.size()!=6U) throw std::invalid_argument(std::string(name)+" requires six values"); for(std::size_t i=0;i<6;++i) target[i]=static_cast<typename std::remove_reference_t<decltype(target)>::value_type>(asU64(a[i],name)); };
  fill_u64_6("num_part_total",m.num_part_total); fill_u64_6("num_part_total_high_word",m.num_part_total_high_word);
  { const auto& a=asArray(member(root,"mass_table"),"mass_table"); if(a.size()!=6U) throw std::invalid_argument("mass_table requires six values"); for(std::size_t i=0;i<6;++i) m.mass_table[i]=asDouble(a[i],"mass_table"); }
  m.box_size=asDouble(member(root,"box_size"),"box_size"); m.scale_factor=asDouble(member(root,"time"),"time"); m.redshift=asDouble(member(root,"redshift"),"redshift"); m.omega_matter=asDouble(member(root,"omega_matter"),"omega_matter"); m.omega_lambda=asDouble(member(root,"omega_lambda"),"omega_lambda"); m.hubble_param=asDouble(member(root,"hubble_param"),"hubble_param");
  { const auto& a=asArray(member(root,"species_policy"),"species_policy"); if(a.size()!=6U) throw std::invalid_argument("species_policy requires six values"); for(std::size_t i=0;i<6;++i) m.species_policy[i]=parseSpeciesPolicy(asString(a[i],"species_policy")); }
  for (const auto& item:asArray(member(root,"fields"),"fields")) {
    IcFieldManifest f; f.source_file_index=static_cast<std::uint32_t>(asU64(member(item,"source_file_index"),"source_file_index")); f.dataset_path=asString(member(item,"dataset_path"),"dataset_path"); f.selected_alias=asString(member(item,"selected_alias"),"selected_alias"); f.scalar_type=asString(member(item,"scalar_type"),"scalar_type"); f.scalar_class=parseScalarClass(asString(member(item,"scalar_class"),"scalar_class")); f.byte_width=static_cast<std::uint8_t>(asU64(member(item,"byte_width"),"byte_width")); f.is_signed=asBool(member(item,"is_signed"),"is_signed"); f.byte_order=parseByteOrder(asString(member(item,"byte_order"),"byte_order")); f.rank=static_cast<std::uint8_t>(asU64(member(item,"rank"),"rank")); for(const auto& d:asArray(member(item,"dimensions"),"dimensions")) f.dimensions.push_back(asU64(d,"dimensions")); f.record_count=asU64(member(item,"record_count"),"record_count"); f.base_unit_to_si=asDouble(member(item,"base_unit_to_si"),"base_unit_to_si"); f.hubble_exponent=asDouble(member(item,"hubble_exponent"),"hubble_exponent"); f.scale_factor_exponent=asDouble(member(item,"scale_factor_exponent"),"scale_factor_exponent"); f.coordinate_frame=parseFrame(asString(member(item,"coordinate_frame"),"coordinate_frame")); f.velocity_convention=parseVelocityConvention(asString(member(item,"velocity_convention"),"velocity_convention")); f.semantics=parseSemantics(asString(member(item,"semantics"),"semantics")); f.disposition=parseDisposition(asString(member(item,"disposition"),"disposition")); f.source_unit=asString(member(item,"source_unit"),"source_unit"); f.target_unit=asString(member(item,"target_unit"),"target_unit"); f.conversion_equation=asString(member(item,"conversion_equation"),"conversion_equation"); m.fields.push_back(std::move(f));
  }
  m.defaulted_fields=stringArray(member(root,"defaulted_fields"),"defaulted_fields"); m.converted_fields=stringArray(member(root,"converted_fields"),"converted_fields"); m.dropped_fields=stringArray(member(root,"dropped_fields"),"dropped_fields"); m.rejected_fields=stringArray(member(root,"rejected_fields"),"rejected_fields"); m.preserved_auxiliary_fields=stringArray(member(root,"preserved_auxiliary_fields"),"preserved_auxiliary_fields"); m.conversion_equations=stringArray(member(root,"conversion_equations"),"conversion_equations"); m.warnings=stringArray(member(root,"warnings"),"warnings"); validateIcManifest(m); return m;
}

IcManifest readIcManifestJson(const std::filesystem::path& input_path) {
  std::ifstream input(input_path,std::ios::binary); if(!input) throw std::runtime_error("failed to open IC manifest: "+input_path.string()); std::ostringstream text; text<<input.rdbuf(); if(!input.eof() && input.fail()) throw std::runtime_error("failed to read IC manifest: "+input_path.string()); return deserializeIcManifestJson(text.str());
}
void writeIcManifestJson(const IcManifest& manifest,const std::filesystem::path& output_path) {
  const auto temporary=output_path.string()+".part"; std::ofstream output(temporary,std::ios::binary|std::ios::trunc); if(!output) throw std::runtime_error("failed to open IC manifest output: "+temporary); output<<serializeIcManifestJson(manifest); output.close(); if(!output) throw std::runtime_error("failed to write IC manifest output: "+temporary); std::error_code error; std::filesystem::rename(temporary,output_path,error); if(error) { std::filesystem::remove(output_path,error); error.clear(); std::filesystem::rename(temporary,output_path,error); } if(error) throw std::runtime_error("failed to atomically finalize IC manifest: "+error.message());
}

}  // namespace cosmosim::io
