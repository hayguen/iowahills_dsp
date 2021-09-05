/*
 This software is part of iowahills_dsp, a set of DSP routines under MIT License.
 2016 By Daniel Klostermann, Iowa Hills Software, LLC  IowaHills.com
 Copyright (c) 2021  Hayati Ayguen <h_ayguen@web.de>
 All rights reserved.

This complex math code is here to help users get started with the IowaHills Filter Code Kit.
This code was not written by a software professional and has not been properly vetted.
It is recommended that this code be replaced with the compilers complex math library.
*/

#ifndef CplxD_HEADER
#define CplxD_HEADER
//---------------------------------------------------------------------------
#include <cmath>


namespace iowahills {
  const double CPLXDMATH_ZERO_TEST = 1.0E-50;   // This is used to test for 0.
  const double CPLXDMATH_ZERO_RESULT = 1.0E50;  // This is returned by the log and log10 functions to avoid a math exception.
}


class CplxD  // complex double
{
  private:   // No private variables.

  public:
  double re; // In a typical compiler math library, these are declared as private, but this requires getter and setter functions. Public access to the real and imag parts is much faster, which is important for FFT routines, for example.
  double im;

  // Constructors
  CplxD();               // The default that gets called when vars are declared.
  CplxD(double, double);

  // Assignment
  void operator = (double X);

  // These operators modify r and i, then return *this
  CplxD & operator += (const CplxD& X);   // CplxD += CplxD
  CplxD & operator += (const double& X);  // CplxD += double
  CplxD & operator -= (const CplxD& X);   // etc
  CplxD & operator -= (const double& X);
  CplxD & operator *= (const CplxD& X);
  CplxD & operator *= (const double& X);
  CplxD & operator /= (const CplxD& X);
  CplxD & operator /= (const double& X);
  bool    operator == (const CplxD& X);
  bool    operator == (const double& X);
  bool    operator != (const CplxD& X);
  bool    operator != (const double& X);


  // These operators require 2 args so must be declared as friends.
  // They return a CplxD value (not *this)
  friend CplxD operator + (const CplxD& X, const CplxD& Y);
  friend CplxD operator + (const double& A, const CplxD& X);
  friend CplxD operator + (const CplxD& X, const double& A);

  friend CplxD operator - (const CplxD& X);   // negative sign, declared as a friend to avoid ambiguity with subtraction
  friend CplxD operator - (const CplxD& X, const CplxD& Y);
  friend CplxD operator - (const double& A, const CplxD& X);
  friend CplxD operator - (const CplxD& X, const double& A);

  friend CplxD operator * (const CplxD& X, const CplxD& Y);
  friend CplxD operator * (const double& A, const CplxD& X);
  friend CplxD operator * (const CplxD& X, const double& A);

  friend CplxD operator / (const CplxD& X, const CplxD& Y);
  friend CplxD operator / (const double& A, const CplxD& X);
  friend CplxD operator / (const CplxD& X, const double& A);

  // std::complex<> has no member conj() or cabs()
  //CplxD  conj() const;   // X.conj()
  //double cabs() const;   // X.cabs()

  // If re and im were private, these would need to be friend functions.
  // put the following global functions into the namespace iowahills
#if 0
  CplxD exp(const CplxD& X);
  CplxD log(const CplxD& X);
  CplxD log10(const CplxD& X);
  CplxD sqrt(const CplxD& X);
  CplxD cos(const CplxD& X);
  CplxD sin(const CplxD& X);
  CplxD tan(const CplxD& X);
  CplxD acos(const CplxD& X);
  CplxD asin(const CplxD& X);
  CplxD atan(const CplxD& X);
  CplxD cosh(const CplxD& X);
  CplxD sinh(const CplxD& X);
  CplxD tanh(const CplxD& X);
  CplxD acosh(const CplxD& X);
  CplxD asinh(const CplxD& X);
  CplxD atanh(const CplxD& X);

  CplxD pow(const CplxD& X, const double& Y);
  CplxD pow(const double& X, const CplxD& Y);
  CplxD pow(const CplxD& X, const CplxD& Y);

  CplxD  conj(const CplxD& X);    // conj(X)
  double abs(const CplxD& X);     // abs(X)  // -> std::complex<> is abs() - not cabs()
  CplxD  Zero(const CplxD& X, const double& Epsilon); // Sets the insignificant part of X to zero.

  double dB(const CplxD& X);           // 20 * log10(|X|)
  double DegreeAngle(const CplxD& X);  // DegreeAngle returns values 180 to -180
  double RadianAngle(const CplxD& X);  // RadianAngle returns radian values Pi to -Pi
  void QuadFormula(double *P, CplxD *Root0, CplxD *Root1);  // Quadratic formula
#endif

};  // end class CplxD



// Note: These function definitions need to be in this header file (not a cpp file)
// in order for the linker to use the inline functions.

// Costructors --------------------------------------------------------------------------------
// Constructors don't return a value.

// Default Constructor. Only called when a var is declared w/o an init.
  inline CplxD::CplxD()
   {
    //re = 0.0;  // Use these to init the vars, including arrays, to 0.
    //im = 0.0;  // Not using these helps find the use of un set vars.
   }

  // Constructor  Used for a declaration with an init:  CplxD A = CplxD(1.0, 0.0);
  // or a simple assignment: Y = CplxD(1.0, 2.0);
  inline CplxD::CplxD(double RealPart, double ImagPart)
   {
    re = RealPart;
    im = ImagPart;
   }


// Assignments -------------------------------------------------------------------

/*
  Copy assignment. Is implicitly defined for all classes.  x = CplxD(1.0, 0.0);
  would be declared as  void operator =(CplxD);
  void CplxD::operator = (CplxD X)
   {
    re = X.re;
    im = X.im;
   }
*/

// Assignment overload  As in X = 0.0;
  inline void CplxD::operator = (double X)
   {
    re = X;
    im = 0.0;
   }

// Addition -------------------------------------------------------------------

// Addition CplxD + CplxD
  inline CplxD operator + (const CplxD& X, const CplxD& Y)
   {
    return CplxD(X.re + Y.re, X.im + Y.im);
   }

// Addition double + CplxD
  inline CplxD operator + (const double& a, const CplxD& Y)
   {
    return CplxD(a + Y.re, Y.im);
   }

// Addition CplxD + double
  inline CplxD operator + (const CplxD& X, const double& a)
   {
    return CplxD(X.re + a, X.im);
   }

// Subtraction -------------------------------------------------------------------

//  CplxD - CplxD
  inline CplxD operator - (const CplxD& X, const CplxD& Y)
   {
    return CplxD(X.re - Y.re, X.im - Y.im);
   }

//  double - CplxD
  inline CplxD operator - (const double& a, const CplxD& X)
   {
    return CplxD(a - X.re, -X.im);
   }

//  CplxD - double
  inline CplxD operator - (const CplxD& X, const double& a)
   {
    return CplxD(X.re - a, X.im);
   }

// Multiplication -------------------------------------------------------------------

//  CplxD * CplxD
  inline CplxD operator * (const CplxD& X, const CplxD& Y)
   {
    return CplxD(X.re * Y.re - X.im * Y.im, X.re * Y.im + X.im * Y.re);
   }

//  CplxD * double
  inline CplxD operator * (const CplxD& X, const double& a)
   {
    return CplxD(X.re * a, X.im * a);
   }

// double * CplxD
  inline CplxD operator * (const double& a, const CplxD& X)
   {
    return CplxD(X.re * a, X.im * a);
   }

// Division -------------------------------------------------------------------

//  CplxD / CplxD
  inline CplxD operator / (const CplxD& X, const CplxD& Y)
   {
    double Mag = Y.re * Y.re + Y.im * Y.im;
    if(Mag < iowahills::CPLXDMATH_ZERO_TEST)
        Mag = iowahills::CPLXDMATH_ZERO_TEST; // to avoid / 0
    double re = ((X.re * Y.re) + (X.im * Y.im)) / Mag;
    double im = ((X.im * Y.re) - (X.re * Y.im)) / Mag;
    return CplxD(re, im);
   }


// CplxD / double
  inline CplxD operator / (const CplxD& X, const double& B)
   {
    if(fabs(B) < iowahills::CPLXDMATH_ZERO_TEST)
        return CplxD(X.re/iowahills::CPLXDMATH_ZERO_TEST, X.im/iowahills::CPLXDMATH_ZERO_TEST); // to avoid / 0
    return CplxD(X.re/B, X.im/B);
   }

//  double / CplxD
  inline CplxD operator / (const double& A, const CplxD& Y)
   {
    double Mag = Y.re * Y.re + Y.im * Y.im;
    if(Mag < iowahills::CPLXDMATH_ZERO_TEST)
        Mag = iowahills::CPLXDMATH_ZERO_TEST; // to avoid / 0
    return( CplxD(A * Y.re / Mag , -A * Y.im / Mag) ) ;
   }


// Compound Addition ---------------------------------------------------------------

// +=  CplxD + CplxD
  inline CplxD & CplxD::operator += (const CplxD& X)
   {
    re += X.re;
    im += X.im;
    return *this;
   }

// +=  CplxD + double
  inline CplxD & CplxD::operator += (const double& A)
   {
    re += A;
    return *this;
   }

// Compound Subtraction ---------------------------------------------------------------

// Negative sign
  inline CplxD operator - (const CplxD& X)
   {
    return CplxD(-X.re, -X.im);
   }

// -=   CplxD - CplxD
  inline CplxD & CplxD::operator -= (const CplxD& X)
   {
    re -= X.re;
    im -= X.im;
    return *this;
   }

// -=  CplxD - double
  inline CplxD & CplxD::operator -= (const double& A)
   {
    re -= A;
    return *this;
   }

// Compound Multiplication ---------------------------------------------------------------

// *=  CplxD * CplxD
  inline CplxD & CplxD::operator *= (const CplxD& X)
   {
    double Temp;
    Temp = X.re * im + X.im * re;
    re = X.re * re - X.im * im;
    im = Temp;
    return *this;
   }

// *=  CplxD * double
  inline CplxD & CplxD::operator *= (const double& A)
   {
    re *= A;
    im *= A;
    return *this;
   }

// Compound Division ---------------------------------------------------------------

//  /=  CplxD / CplxD
  inline CplxD & CplxD::operator /= (const CplxD& X)
   {
    double Temp, Mag = X.re * X.re + X.im * X.im;
    Temp = ((X.re * re) + (X.im * im)) / Mag;  // Use Temp because we still need re.
    im =   ((X.re * im) - (X.im * re)) / Mag;
    re = Temp;
    return *this;
   }


//  /=  CplxD / double
  inline CplxD & CplxD::operator /= (const double& A)
   {
    re /= A;
    im /= A;
    return *this;
   }

// Comparisons -------------------------------------------------------------------

// CplxD == CplxD
  inline bool CplxD::operator == (const CplxD& X)
   {
    return( re == X.re && im == X.im );
   }

// CplxD == double
  inline bool CplxD::operator == (const double& X)
   {
    return( re == X && im == 0.0 );
   }

// CplxD != CplxD
  inline bool CplxD::operator != (const CplxD& X)
   {
    return( re != X.re || im != X.im );
   }

// CplxD != CplxD
  inline bool CplxD::operator != (const double& X)
   {
    return( re != X || im != 0.0 );
   }

#if 0
// Complex Conjugate  ------------------------------------------------------------

  inline CplxD CplxD::conj() const   // as in X.conj()
   {
    return(CplxD(re, -im));
   }

// Absolute Value ------------------------------------------------------------

  inline double CplxD::abs() const // as in X.abs()
   {
    return( std::sqrt(re * re + im * im) ); // std::sqrt in this member function so there isn't ambiguity with sqrt(Cplx)
   }
#endif

namespace iowahills {

// Complex Conjugate  ------------------------------------------------------------

  inline CplxD conj(const CplxD& X) // as in conj(X)
   {
    return(CplxD(X.re, -X.im) );
   }

  // Absolute Value ------------------------------------------------------------

  inline double abs(const CplxD& X) // as in abs(X)
   {
    return( sqrt(X.re * X.re + X.im * X.im) );
   }

// Exponentials  ---------------------------------------------------------------------------------

// e^X   X is a radians arg
  inline CplxD exp(const CplxD& X)
   {
    double Mag = std::exp(X.re);
    return CplxD(Mag * cos(X.im), Mag * sin(X.im) );
   }

// log(X)  = ln(CplxD X)  This returns -CPLXDMATH_ZERO_RESULT for log(0) instead of throwing a math exception.
  inline CplxD log(const CplxD& X)
   {
    double Angle, Mag;
    Mag = X.re * X.re + X.im * X.im;
    if(Mag < CPLXDMATH_ZERO_TEST) // taking log(0)
	 {
	  return CplxD(-CPLXDMATH_ZERO_RESULT, 0.0);
	 }
    Angle = atan2(X.im, X.re);  // atan2(Y,X) returns +/- PI/2 if Y != 0 && X = 0
    return CplxD(std::log(Mag) / 2.0, Angle);
   }

// log10(CplxD X)   This returns -CPLXDMATH_ZERO_RESULT for log(0) instead of throwing a math exception.
  inline CplxD log10(const CplxD& X)
   {
    double Angle, Mag;
    Mag = X.re * X.re + X.im * X.im;
    if(Mag < CPLXDMATH_ZERO_TEST)  // taking log(0)
	 {
	  return CplxD(-CPLXDMATH_ZERO_RESULT, 0.0);
	 }
    Angle = atan2(X.im, X.re);  // atan2(Y,X) returns +/- PI/2 if Y != 0 && X = 0
    return CplxD(std::log10(Mag) / 2.0, Angle / std::log(10.0) );  // The constant M_LN10 doesn't have enough significance.
   }


// X^Y  CplxD ^ CplxD
  inline CplxD pow(const CplxD& X, const CplxD& Y)
   {
    double Angle, Mag, Radius, Theta;
    if(abs(Y) < CPLXDMATH_ZERO_TEST) // anything, inc 0, to the 0th pow = 1
     {
      return CplxD(1.0, 0.0);
     }
    if(abs(X) < CPLXDMATH_ZERO_TEST) // 0 to any power, except 0, = 0
     {
      return CplxD(0.0, 0.0);
     }
    Mag = X.re * X.re + X.im * X.im;
    Angle = atan2(X.im, X.re);        // atan2(Y,X) returns +/- PI/2 if Y != 0 && X = 0
    Radius = std::exp(Y.re / 2.0 * std::log(Mag) - Y.im * Angle);
    Theta = Y.re * Angle + Y.im / 2.0 * std::log(Mag); // Mag = 0 was checked for above
    return CplxD( Radius * cos(Theta), Radius * sin(Theta) );
   }

// X^Y  re ^ CplxD      Same code as X^Y
  inline CplxD pow(const double& A, const CplxD& Y)
   {
    double Angle, Mag, Radius, Theta;
    if(abs(Y) < CPLXDMATH_ZERO_TEST) // anything, inc 0, to the 0th pow = 1
     {
      return CplxD(1.0, 0.0);
     }
    if(fabs(A) < CPLXDMATH_ZERO_TEST)  // 0 to any power, except 0, = 0
     {
      return CplxD(0.0, 0.0);
     }
    Mag = A*A;
    Angle = atan2(0.0, A);        // atan2(Y,X) returns +/- PI/2 if Y != 0 && X = 0
    Radius = std::exp(Y.re / 2.0 * std::log(Mag) - Y.im * Angle);
    Theta = Y.re * Angle + Y.im / 2.0 * std::log(Mag); // Mag = 0 was checked for above
    return CplxD( Radius * cos(Theta), Radius * sin(Theta) );
   }

// X^Y CplxD ^ re
  inline CplxD pow(const CplxD& X, const double& B)
   {
    double Theta, Radius;
    if(fabs(B) < CPLXDMATH_ZERO_TEST) // anything, inc 0, to the 0th pow = 1
     {
      return CplxD(1.0, 0.0);
     }
    if(abs(X) < CPLXDMATH_ZERO_TEST) // 0 to any power, except 0, = 0
     {
      return CplxD(0.0, 0.0);
     }
    Theta = B * atan2(X.im,X.re); // atan2(Y,X) returns +/- PI/2 if X = 0  errs with X=Y=0
    Radius = std::pow( sqrt(X.re * X.re + X.im * X.im), B);
    return CplxD( Radius * cos(Theta), Radius * sin(Theta) );
   }

// Square Root of X
  inline CplxD sqrt(const CplxD& X)
   {
    double Theta, Radius;
    if(abs(X) < CPLXDMATH_ZERO_TEST) // sqrt(0) = 0 but atan2(0,0) throws an exception
     {
      return CplxD(0.0, 0.0);
     }
    Theta = atan2(X.im,X.re) / 2.0;     // atan2(Y,X) returns +/- PI/2 if X = 0
    Radius = std::pow(X.re * X.re + X.im * X.im, 0.25) ;
    return CplxD( Radius * cos(Theta), Radius * sin(Theta) );
   }


// Trig Functions---------------------------------------------------------------------

// cosine
  inline CplxD cos(const CplxD& X)
   {
    return CplxD(std::cos(X.re) * cosh(X.im), -sin(X.re) * sinh(X.im) );
   }

// sine
  inline CplxD sin(const CplxD& X)
   {
    return CplxD(std::sin(X.re) * cosh(X.im), std::cos(X.re) * sinh(X.im) );
   }

// tangent
  inline CplxD tan(const CplxD& X)
   {
    double A, B, C;
    A = std::cos(X.re) * std::sin(X.re);                              // A = 0 at X.re = +/- Pi/2
    B = cosh(X.im) * sinh(X.im);                            // B = +/-5.77 at X.im = +/- Pi/2
    C = cosh(X.im) * cosh(X.im) - std::sin(X.re) * std::sin(X.re);    // C = 0 at X.re = +/- n*Pi/2 and X.im = 0
    if(fabs(C) < CPLXDMATH_ZERO_TEST)
	 {
      return CplxD(std::tan(X.re), 0.0);
	 }
    return CplxD(A/C, B/C);
   }

// arc cosine
  inline CplxD acos(const CplxD& X)
   {
    CplxD Root, Temp;
    Root = sqrt(1.0 - X*X);
    Temp = log( X + CplxD(-Root.im, Root.re) ); // swapping the re and im parts == i * Root
    return CplxD(Temp.im, -Temp.re);
   }

// arc sine
  inline CplxD asin(const CplxD& X)
   {
    CplxD Root, Temp;
    Root = sqrt(1.0 - X*X);
    Temp = log( CplxD(-X.im + Root.re, X.re + Root.im) );
    return CplxD(Temp.im, -Temp.re);
   }

// arc tangent
  inline CplxD atan(const CplxD& X)
   {
    CplxD Temp;
    Temp = CplxD(X.re, 1.0 + X.im) / CplxD(-X.re, 1.0 - X.im);
    Temp = log(Temp) / 2.0;
    return CplxD(-Temp.im, Temp.re);
   }

// Hyperbolics -----------------------------------------------------------------------

  inline CplxD cosh(const CplxD& X)
   {
    return CplxD(std::cosh(X.re) * std::cos(X.im), sinh(X.re) * std::sin(X.im) );
   }

  inline CplxD sinh(const CplxD& X)
   {
    return CplxD(std::sinh(X.re) * std::cos(X.im), std::cosh(X.re) * std::sin(X.im) );
   }

  inline CplxD tanh(const CplxD& X)
   {
    CplxD Temp1, Temp2;
    Temp1 = CplxD( std::tanh(X.re), std::tan(X.im) );
    Temp2 = CplxD(1.0, std::tanh(X.re) * std::tan(X.im));
    Temp1 /= Temp2;  // Temp2 cannot = 0
    return Temp1;
   }


// arc cosh  = +/- j*acos(X)
  inline CplxD acosh(const CplxD& X)
   {
    CplxD Root, Temp;
    Root = sqrt(1.0 - X*X);
    Temp = log( X + CplxD(-Root.im, Root.re) ); // swapping the re and im parts == i * Root
    if(Temp.re < 0.0)Temp = -Temp;              // There are 2 solutions, we return the positive.
    return Temp;
   }

// arc sinh  = j*asin(-j*X)
  inline CplxD asinh(const CplxD& X)
   {
    CplxD Y, Root, Temp;
    Y = CplxD(X.im, -X.re); // -j*X
    Root = sqrt(1.0 - Y*Y);
    Temp = log( CplxD(-Y.im + Root.re, Y.re + Root.im) );
    return Temp;
   }

// arc tanh
  inline CplxD atanh(const CplxD& X)
   {
    CplxD Temp;
    Temp = CplxD(1.0 + X.re, X.im) / CplxD(1.0 - X.re, -X.im);
    Temp = log(Temp) / 2.0;
    return Temp;
   }

//---------------------------------------------------------------------------

  inline double dB(const CplxD& X) // 20 * log10(|x|)
   {
    if(fabs(X.re) > 1.0E10 || fabs(X.im) > 1.0E10)  return(200.0);
    double Mag = X.re * X.re + X.im * X.im;
    if(Mag <= 1.0E-10)  return(-200.0); // Don't want to return a val that will overflow an int.
    Mag = 10.0 * std::log10(Mag);         // This is 10 * because we didn't take sqrt of mag.
    return(Mag);
   }

  inline double DegreeAngle(const CplxD& X) // DegreeAngle returns values 180 to -180
   {
    if(X.re == 0.0 && X.im == 0.0)return(0.0);
    return( atan2(X.im, X.re) * 180.0 / M_PI ); // atan2(y,x) returns +/- PI/2 if x = 0
   }

  inline double RadianAngle(const CplxD& X)// RadianAngle returns radian values Pi to -Pi
   {
    if(X.re == 0.0 && X.im == 0.0)return(0.0);
    return( atan2(X.im, X.re) ); // atan2(y,x) returns +/- PI/2 if x = 0
   }

//---------------------------------------------------------------------------

  // Use this to zero out an insignificant real or imag part.
  // Used to set nearly real roots to pure reals and nearly imag roots to pure imags.
  inline CplxD Zero(const CplxD& X, const double& Epsilon) // as in X = Zero(X, Epsilon);
   {
    CplxD Y = X;
    if(fabs(X.re) < fabs(X.im) * Epsilon)   Y.re = 0.0;
    if(fabs(X.im) < fabs(X.re) * Epsilon)   Y.im = 0.0;
    return Y;
   }

//---------------------------------------------------------------------------

  // Calc the roots of ax^2 + bx + c = P[0]*x^2 + P[1]*x + P[2]
  inline void QuadFormula(double *P, CplxD *Root0, CplxD *Root1)
   {
    Root0->re = Root1->re = Root0->im = Root1->im = 0.0;
    if(P[0] == 0.0 && P[1] == 0.0) // no roots
     {
      return;
     }
    if(P[0] == 0.0) // 1st order eq
     {
      Root0->re = -P[2]/P[1];  // P[1]=0 was tested for above.
      return;
     }

    // Quadratic formula
    if(P[0] != 1.0)
     {
      P[1] /= P[0];
      P[2] /= P[0];
      P[0] = 1.0;
     }

    double D = P[1]*P[1] - 4.0*P[2];
    if(D >= 0.0)  // 2 real roots
     {
      Root0->re = (-P[1] - std::sqrt(D)) * 0.5;
      Root1->re = (-P[1] + std::sqrt(D)) * 0.5;
     }
    else // D < 0.0   2 complex roots
     {
      Root0->re = Root1->re = -P[1] * 0.5;
      Root0->im = std::sqrt(-D) * 0.5;
      Root1->im = -Root0->im;
     }
   }

} // namespace iowahills

//---------------------------------------------------------------------------

#endif

