#include "sumpowerspoly.h"

/** number of powers for which the polynomial is builtin*/
static const int intrinsicCoeffs=2;

SumPowersPoly::SumPowersPoly(unsigned power):
  power(power),
  firsterm(1,power+1),
  numCoeffs((power-2)/2)
{
  coeffs.reserve(numCoeffs);
  //we will call the factory and as such with recurse into this constructor.
  if(power>1){

  } else {
 //one coefficient which is builtin (firsterm)
 //2nd coefficient is builtin (1/2)
  }
}

const SumPowersPoly &SumPowersPoly::factory(unsigned power)
{
  //if we already have constructed one return it
  if(power< cache.size()){
    return *cache[power];
  }
  //else make a new one first making all the lower powers for it
  for(int pow=0;pow<power;++pow){
    factory(pow);
  }
  cache.push_back(new SumPowersPoly(power));
  return *cache.back();
}

int SumPowersPoly::sum(int n, int *higherpower) const
{
    if(power<0){
        return -1;//hope this will blow hard
    }
    if(n<0){
        return 0;//probably a trivial user error, this is probably a good response for that
    }
    if(power==0){
        return n;
    }
    if(n==0){
        return 0;
    }
#if 0 //don't special case these until the code is fully debugged, they are handy for showing up coding errors.
    if(n==1){
        return 1;
    }
    if(n==2){
        return 1+2<<power;
    }
#endif
    //////////////////
    int nsquared=n*n;//FUE.

    if(power==1){//special case makes rest of code easier to follow
        if(higherpower){
            *higherpower=nsquared;
        }
        return (n+nsquared)>>1;//ignoring firsterm as we know it is == 1/2.
    }

    Rational summer;
    int poweredup=(power&1)?nsquared:n;//odd powers have a higher power of n in each term than their sibling even one
//     if we have
    for(int x=numCoeffs;x-->0;){

        //summer+=coeffs[x]*poweredup;;
        if(x){//all but last iteration boost the power
            poweredup*=nsquared;
        }
    }
    poweredup *= n;
    summer+= poweredup/2;
    poweredup *= n;
    summer+= poweredup * firsterm;

    if(higherpower){
        *higherpower=poweredup;
    }
    return summer;
}
