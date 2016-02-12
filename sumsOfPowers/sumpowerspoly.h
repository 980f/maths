#ifndef SUMPOWERSPOLY_H
#define SUMPOWERSPOLY_H

#include <vector>  //used for caching polys, a necessary thing
using namespace std;

#include "ratioofproducts.h"
/** may eventually use a fancier infinite precision rational number class. */
typedef RatioOfProducts Rational;

class SumPowersPoly
{
    /** we need to cache because: 1: the content of each poly is const and expensive to construct; 2: we need to construct lower degrees to get to higher degrees and we don't want to force the user to know that, and we also don't want to discard all of that work since usually all the lower powers are needed.*/
    static vector<SumPowersPoly *> cache;
    /** retain construction arg verbatim */
    const int power;
    /** cache 1/(power+1) */
    const Rational firsterm;
    /** cache number of terms beyond first two, is size of coeffs but we don't want to have to constantly remember that */
    const int numCoeffs;
    /** the non-trivial coefficients, i.e. those that are not embedded inline */
    vector<const Rational> coeffs;
private: //only the factory knows when one should be created.
    /** the constructor is very non-trivial, it may create many other objects of this class before creating the requested one.*/
    SumPowersPoly(unsigned power);
public:
    /** because all methods are const and we cache all created objects we don't offer a constructor, just a factory that returns references to shared objects.
    This may cause many constructions and reallocations so you might want to get a factory object for the largest power that you might use early in your application's execution.*/
    static const SumPowersPoly& factory(unsigned power);
  /** @returns sum of integers raised to the constructed power @param n highest integer.
      using int instead of unsigned to detect erroneous use, instead of computing a number that would overflow. For negative values the sum will be 0, a fairly silent failure.
      return int instead of long because we will likely have internal overflows before we get close to a sum that requires it, unless we use longs througout which we won't
      Also not using float for the return because we are only using floats initially and the number is theoretically integral.
      If @param higherpower is not null then n^(power+1) will be written to it. This is rarely usefull but nice for initial debug of code.*/
  int sum(int n,int *higherpower=nullptr) const;
};

#endif // SUMPOWERSPOLY_H
