/******************************************************************************
 * Filename : datastructures.hpp
 * Date     : 2018 Aug 25
 * Author   : Willy Clarke (willy@clarke.no)
 * Version  : 0.0.1
 * Copyright: W. Clarke
 * License  : MIT
 * Descripti: Data structures for My Raytracing Challenge.
 ******************************************************************************/
#ifndef COMMON_DATASTRUCTURES_HPP
#define COMMON_DATASTRUCTURES_HPP

#include <iomanip>  // for setw().
#include <iostream>
#include <limits>
#include <strstream>
#include <vector>

//------------------------------------------------------------------------------
// NOTE: The assert will write to a null pointer.
#if HANDMADE_SLOW
#define debug(format, ...) fprintf(stderr, format, ##__VA_ARGS__)
#define Assert(Condition, ...)                                          \
   if (!(Condition))                                                    \
   {                                                                    \
      fprintf(stderr, "ASSERT. Function %s. Line %d\n", ##__VA_ARGS__); \
      *(volatile int *)0 = 0;                                           \
   }
#else
#define Assert(Condition, ...)
#endif

namespace ww
{
constexpr float EPSILON = 0.0000001;  // 1E27 * std::numeric_limits<float>::min();

//------------------------------------------------------------------------------
union tup {
   tup() : X{}, Y{}, Z{}, W{} {};
   tup(float IX, float IY, float IZ, float IW) : X{IX}, Y{IY}, Z{IZ}, W{IW} {};
   ~tup() {}
   tup(tup const &Other)
   {
      X = Other.X;
      Y = Other.Y;
      Z = Other.Z;
      W = Other.W;
   }
   struct
   {
      float X{};
      float Y{};
      float Z{};
      float W{};  //!< is 1.0 when tuple is point and 0.0 when tuple is a vector.
   };
   struct
   {
      float R;
      float G;
      float B;
      float I;  //!< Intensity is 1.0 at max and 0.0 at pitch black.
   };
};  // end of union tup.

//------------------------------------------------------------------------------
struct canvas
{
   canvas(int IW = 10, int IH = 10) : W{IW}, H{IH} { vXY.resize(W * H); }
   canvas(canvas const &Other)
   {
      W = Other.W;
      H = Other.H;
      vXY = Other.vXY;
   }
   ~canvas() {}
   int W{};  //<! Width
   int H{};  //<! Height
   std::vector<tup> vXY{};
};

// NOTE: Declarations.
tup Add(tup const &A, tup const &B);
tup Color(float const R, float const G, float const B);
tup Cross(tup const &A, tup const &B);

// ---
// NOTE: For an explanation of what the dot product actually is you can take a look
//     : at the description at
//     : http://betterexplained.com/articles/vector-calculus-understanding-the-dot-product .
// ---
float Dot(tup const &A, tup const &B);

// ---
// NOTE: For a discussion on how to do comparison with floating point number the following web
//     : site can be consulted:
//     https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
// ---
bool Equal(float const A, float const B);
bool Equal(tup const &A, tup const &B);
bool IsVector(tup const &Tup);
bool IsPoint(tup const &Tup);
float MagSquared(tup const &Tup);
float Mag(tup const &Tup);
tup Multiply(float const S, tup const &Tup);

// NOTE: Also called Hadamard or Schur product.
tup Multiply(tup const A, tup const B);

tup Negate(tup const &Tup);
tup Normal(tup const &Tup);
tup Point(float A, float B, float C);
tup Sub(tup const &A, tup const &B);
tup Vector(float A, float B, float C);

// ---
// NOTE: Canvas declarations.
// ---
void WritePixel(canvas &Canvas, int X, int Y, tup const &Color);
tup ReadPixel(canvas &Canvas, int X, int Y);
std::strstream PPMHeader(canvas const &Canvas, int X, int Y);
std::string PPMHeader(canvas const &Canvas);
int WriteToPPM(canvas const &Canvas, std::string const &Filename = "test.ppm");
void WriteToPPMFile(canvas const &Canvas, std::string const &Filename = "test.ppm");
std::shared_ptr<canvas> ReadFromPPM(std::string const &Filename = "test.ppm");
};  // namespace ww

// ---
// NOTE: Declare the operator for inclusion elsewhere.
// ---
std::ostream &operator<<(std::ostream &stream, const ww::tup &T);
ww::tup operator+(ww::tup const &A, ww::tup const &B);
ww::tup operator-(ww::tup const &Tup);
ww::tup operator-(ww::tup const &A, ww::tup const &B);
ww::tup operator*(float const S, ww::tup const &Tup);
ww::tup operator*(ww::tup const &Tup, float const S);
ww::tup operator*(ww::tup const &A, ww::tup const &B);
ww::tup operator/(ww::tup const &Tup, float const S);
#endif
