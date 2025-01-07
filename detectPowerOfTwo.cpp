#include <iostream>
#include <bit>    //for count leading zeroes
#include <limits> //for number of bits in type, needed for count leading zeroes.

template<typename Unsigned> auto floorlog2(Unsigned x) {
  return std::numeric_limits<decltype(x)>::digits + ~std::countl_zero(x);
}

int main(int argc, char *argv[]) {
  unsigned bit = argc > 1 ? atoi(argv[1]) : 8;
  unsigned limit = 1 << bit;
  unsigned mask = limit - 1;

  for (unsigned i = 0; i < limit; i++) {
    std::cout << std::hex;
    std::cout << i << '\t';
    unsigned accum = (-i ^ i);
    std::cout << accum << '\t';
    accum += i << 1;
    std::cout << accum << '\t';
    bool byXOR = 0 == accum;
    std::cout << byXOR << '\t';
    std::cout << (0 == (mask & ((i ^ -i) + (i << 1)))) << '\t';
    auto clz1 = floorlog2(i);
    auto clz2 = floorlog2(i - 1);
    bool byCLZ = clz1 != clz2;
    std::cout << byCLZ << '\t';
    bool byCLZ2 = i == (1 << clz1);
    std::cout << byCLZ2 << '\t';
    if (byXOR != byCLZ) {
      std::cout << std::dec << int(-i) << '\t';
    }
    std::cout << std::endl;
  }
  return 0;
}
