#include "ratioofproducts.h"



////////////
/// \brief RatioOfProducts::RatioOfProducts
///

RatioOfProducts::RatioOfProducts()
{

}

void RatioOfProducts::operator *=(const RatioOfProducts &other)
{
  numerator*=other.numerator;
  denominator *= other.denominator;
  //todo: remove common terms.
}

////////////////////////
/// \brief RatioOfProducts::Product::Product
///

void RatioOfProducts::Product::removeSplit(const RatioOfProducts::SignPower &split)
{
  pow2-=split.pow2; //remove 2's
  if(split.isNegative){//remove sign
    toggleSign();\
  }
}

RatioOfProducts::Product::Product():SignPower()
{
  //vector should initialize as empty.
}

int RatioOfProducts::Product::asInt() const
{
  if(zeroed){
    return 0;
  }
  //////////////
  int accumulator=1;
  for(unsigned x:factors){
    accumulator*=x; //ignoring overflows
  }
  return applyTo(accumulator);
}

RatioOfProducts::Product &RatioOfProducts::Product::operator *=(int term)
{
  if(zeroed){
    //do nothing, once zeroed it stays zero.
  } else {
    if(term==0){
      zeroed=true;
      //clean for debug, should be ignored.
      isNegative=false;
      pow2=0;
    } else {
      PsuedoFloat split(term);
      if(split.magnitude>1){
        factors.push_back(split.magnitude);
      }
      times(split);//sign and pow2
    }
  }
  return *this;
}

bool RatioOfProducts::Product::removeTerm(int term)
{
  if(term==0){
    //refuse rather than zero
    return false;
  }
  PsuedoFloat split(term);
  if(split.pow2>pow2){//not enough 2's.
    return false;
  }
  //but don't divide by 2's yet.
  for(int x=factors.size();x-->0;){
    if(factors[x] % split.magnitude == 0){
      factors[x]/=split.magnitude;
      removeSplit(split);
      if(factors[x]==1){
        factors.erase(factors.begin()+x);
      }
      return true;
    }
  }
  //until we do prime factorizations we can't profitably remove term from multiple separate factors.
  return false;
}

RatioOfProducts::Product &RatioOfProducts::Product::timesFactorial(int n)
{
  while(n>2){
    operator *=(n--);
  }
  if(n==2){
    ++pow2;
  }
  return *this;
}


RatioOfProducts &RatioOfProducts::timesRatio(int n, int d)
{
  //try divide first, extends high range
  if(!numerator.removeTerm(d)){
    denominator*=d;
  }
  numerator*=n;
  return *this;
}


RatioOfProducts::PsuedoFloat RatioOfProducts::PsuedoFloat::gcd(const RatioOfProducts::PsuedoFloat &one,const RatioOfProducts::PsuedoFloat &other)
{
  PsuedoFloat gcd;
  gcd.pow2= min(one.pow2,other.pow2);//relies upon other code keeping pow2 >=0.
  gcd.magnitude=1;
  //we now ignore the other power of 2.

  unsigned a=one.magnitude;
  unsigned b=other.magnitude;

  while(a!=b && a>1 && b>1) {//compares to 1 are significant performance optimizations
    if(a>b){
      a-=b;
      a>>=CTZ(a);
    } else {
      b-=a;
      b>>=CTZ(b);
    }
  }
  if(a==b){
    gcd.magnitude=a;
  }
  //else one of them is 1 and that means the magnitudes were coprime and our init value of gcd= 1^(lesser power) is the proper answer.
  return gcd;
}
