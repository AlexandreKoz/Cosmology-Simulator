#include <sstream>
#include <string>

#include "bench/reporting/bench_report.hpp"

int main() {
  cosmosim::bench::BenchmarkReporter reporter("bench_smoke");
  auto execution = cosmosim::bench::defaultExecutionConfig(1, 2);
  cosmosim::bench::addExecutionFields(reporter, execution);
  reporter.addField("sample_counter", 42);
  cosmosim::bench::addBandwidthFields(reporter, 1024, 0.5);

  std::ostringstream out;
  reporter.write(out);
  const std::string rendered = out.str();

  if (rendered.find("bench_smoke") == std::string::npos) {
    return 1;
  }
  if (rendered.find("sample_counter=42") == std::string::npos) {
    return 1;
  }
  if (rendered.find("effective_bandwidth_gb_s=") == std::string::npos) {
    return 1;
  }
  return 0;
}
