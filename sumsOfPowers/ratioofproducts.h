#ifndef RATIOOFPRODUCTS_H
#define RATIOOFPRODUCTS_H

/** a simple infinite precision rational number that does not handle addition nor does it guarantee maximum possible range.
It was created to serve the coefficients of the sums of powers of the integers */

#include <vector>
using namespace std;


class RatioOfProducts
{
    /** shared part of number splitter and product */
    class SignPower {
    protected: //hiding so we can slip in a per-platform asm implementation.
      /** power of 2 == number of trailing zeroes */
      static int CTZ(const int term) {
        if(term==0){
          return 0;
        }
        int pow2=0;
        while(0==(term&1<<pow2))++pow2;
        return pow2;
      }
    public:
      bool isNegative;
      int pow2;

      int split(int term){
          isNegative=term<0;
          if(isNegative){
              term=-term;
          }
          pow2=CTZ(term);
          return term>>pow2;
      }

      int applyTo(unsigned magnitude)const{
          int term=magnitude<<pow2; //might overflow!
          return isNegative?-term:term; //user beware on 0x800....
      }

      inline void toggleSign(){
          isNegative=!isNegative;
      }

      SignPower(int init=0):isNegative(false),pow2(0){
          if(init){
              split(init);
          }
      }

      void times(const SignPower & other){
          if(other.isNegative){
              toggleSign();
          }
          pow2+=other.pow2;
      }

    };

    /** a slightly floaty number */
    class PsuedoFloat: public SignPower {
      public:
        unsigned magnitude;

        PsuedoFloat(int init=0):SignPower(init),magnitude(0){
            magnitude=(isNegative?-init:init)>>pow2;
        }

        //note: default copy constructor is fine.

        /** core transformation: */
        PsuedoFloat &operator =(int term){
            magnitude=split(term);
            return *this;
        }

        PsuedoFloat & operator *=(const PsuedoFloat &other){
            magnitude*=other.magnitude;
            times(other);
            return *this;
        }

        int asInt() const {
            return applyTo(magnitude);
        }

        /** @returns whether effective value of this is not zero */
        operator bool() const {
          return magnitude!=0;//can ignore pow2, asInt will return 0 if magnitude is 0
        }

        /** @returns gcd of one and other */
        static PsuedoFloat gcd(const PsuedoFloat &one, const PsuedoFloat &other);
    };


    class Product: public SignPower  {
        /** once any factor is zero the product is zero forever more */
        bool zeroed;
        /** the unlimited number of odd factors */
        vector<unsigned> factors;

        void removeSplit(const SignPower &split);
    public:
        Product();
        int asInt()const;
        Product & operator *=(int term);
        Product & operator *=(const Product &term);
        /** @returns whether term could be removed, if so then it has been removed */
        bool removeTerm(int term);

        /** multiply by factorial of @param n, @returns this*/
        Product& timesFactorial(int n);
    };

    Product numerator;
    Product denominator;

public:
  /** creates 1/1 */
  RatioOfProducts();
  /** handy for creating const ratioes */
  RatioOfProducts(int n,int d);

  /** try to reduce number of terms by 'cancelling' terms common to num and denom. */
  void reduce(int howhardtotry=0);

  /** add factors to top and bottom, removing pairs when possible */
  void operator *=(const RatioOfProducts &other);
  /** add with some attempt at keeping the denominator from being outrageously long */
  void operator +=(const RatioOfProducts &other);


  void operator *=(int term);
  void operator /=(int term);

  /** multiply by factorial of @param n, @returns this*/
  RatioOfProducts& timesRatio(int n,int d);

  /** multiply by factorial of @param n, @returns this*/
  RatioOfProducts& timesFactorial(int n);
  /** divide by factorial of @param n, @returns this */
  RatioOfProducts& overFactorial(int n);

  double ratio()const;
  /** @returns signed integer that is somewhere close to the theoretical value */
  int truncated()const;
  /** @returns product of the signs, may be spurious if denominator became zero */
  bool isNegative()const;
};

#endif // RATIOOFPRODUCTS_H
