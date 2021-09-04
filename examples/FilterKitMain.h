
#ifndef FilterKitMainH
#define FilterKitMainH


 #define NUM_SAMPLES 2048

/*
 void ExampleFIRCall(void);
 void ExampleIIRCall(void);
 void ExampleLowPassProtoType(void);
 void ExampleFreqSampling(void);
 void ExampleFreqSampling2(void);
 void AnalogFreqSampling(void);
*/

#endif

//---------------------------------------------------------------------------

/*
All the constants used in the kit are here. We believe the primary ones (M_PI, M_SQRT2, etc)
are supposed to be defined in the math.h file for all c compilers, but if they are not
in your math.h, or float.h, here is how they are defined.

#define M_LN2       0.693147180559945309417 // log(2) = ln(2)
#define LOG_OF_TWO  0.301029995663981184    // log10(2)
#define ARC_DB_HALF 0.316227766016837952
#define DBL_MIN     2.22507385850720E-308  // from float.h
#define DBL_MAX     1.79769313486232E+308
#define FLT_MIN     1.175494E-38
#define FLT_MAX     3.402823E+38

#define ZERO_PLUS   8.88178419700125232E-16   // 2^-50 = 4*DBL_EPSILON
#define ZERO_MINUS -8.88178419700125232E-16

// from float.h  Epsilon is the smallest value that can be added to 1.0 so that 1.0 + Epsilon != 1.0
#define LDBL_EPSILON   1.084202172485504434E-019L // = 2^-63
#define DBL_EPSILON    2.2204460492503131E-16     // = 2^-52
#define FLT_EPSILON    1.19209290E-07F            // = 2^-23

#define M_E         2.71828182845904523536      // natural e
#define M_PI        3.14159265358979323846      // Pi
#define M_2PI	    6.28318530717958647692      // 2*Pi
#define M_PI_2      1.57079632679489661923      // Pi/2
#define M_PI_4      0.785398163397448309616     // Pi/4
#define M_1_PI      0.318309886183790671538     // 0.1 * Pi
#define M_SQRT2     1.41421356237309504880      // sqrt(2)
#define M_SQRT_2    0.707106781186547524401     // sqrt(2)/2
#define M_SQRT3     1.7320508075688772935274463 // sqrt(3) from Wikipedia
#define M_SQRT3_2   0.8660254037844386467637231 // sqrt(3)/2

*/
