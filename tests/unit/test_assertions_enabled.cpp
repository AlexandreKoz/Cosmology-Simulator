#ifdef NDEBUG
#error "CosmoSim test targets must keep standard assertions enabled"
#endif

#include <cassert>

int main() {
  assert(true);
  return 0;
}
