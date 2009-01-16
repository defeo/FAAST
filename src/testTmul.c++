/* Courtesy of E. Schost */

#include "Tmul.hpp"
#include <assert.h>
#include <NTL/ZZ_pXFactoring.h>

// to check
zz_p cnv(ZZ_p a){
  return to_zz_p(to_long(a._ZZ_p__rep));
}
// to check
zz_pX cnv(ZZ_pX& A){
  zz_pX a;
  for (long i = 0; i <= deg(A); i++)
    SetCoeff(a, i, cnv(coeff(A, i)));
  return a;
}
// to check
GF2 cnv2(ZZ_p a){
  return to_GF2(to_long(a._ZZ_p__rep));
}
// to check
GF2X cnv2(ZZ_pX& A){
  GF2X a;
  for (long i = 0; i <= deg(A); i++)
    SetCoeff(a, i, cnv2(coeff(A, i)));
  return a;
}

int main(int argc, char** argv){
  zz_p::init(2);
  ZZ_p::init(to_ZZ(2));

  ZZ_pX F, A, B, X;
  ZZ_pXModulus Fmod;
  ZZ_pXMultiplier Bmul;

  zz_pX f, a, b, x;
  zz_pXModulus fmod;
  zz_pXMultiplier bmul;

  GF2X f2, a2, b2, x2;
  GF2XModulus f2mod;
  GF2XTransMultiplier b2mul;
  
  const long d = 200;
  int col = 1;

  while(1){
    // builds a polynomial in F_2[x] as a ZZ_pX
    F=BuildIrred_ZZ_pX(d);
    // and the modulus
    Fmod = ZZ_pXModulus(F);
    
    // converts F to zz_pX
    f=cnv(F);
    // get the modulus
    fmod = zz_pXModulus(f);
    
    // converts F to GF2X
    f2=cnv2(F);
    // get the modulus
    f2mod = GF2XModulus(f2);
    
    // random polys + multiplier
    A=random_ZZ_pX(d);
    B=random_ZZ_pX(d);
    build(Bmul, B, Fmod);
    
    // same over zz_pX
    a=cnv(A);
    b=cnv(B);
    build(bmul, b, fmod);
    
    // same over GF2X
    a2=cnv2(A);
    b2=cnv2(B);
    build(b2mul, b2, f2mod);
    
    // does transposed product
    TransMulMod(X, A, Bmul, Fmod);
    TransMulMod(x, a, bmul, fmod);
    TransMulMod(x2, a2, b2mul, f2mod);
    
    // check the results
    assert (cnv(X)==x);
   	assert (cnv2(X) == x2);
    cout << ".";
    if (col == 0) cout << endl;
    cout.flush();
    col++; col %= 80;
  }

  return 0;
}
