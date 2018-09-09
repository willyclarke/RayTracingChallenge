/******************************************************************************
 * Filename : datastructures.cpp
 * Date     : 2018 Aug 26
 * Author   : Willy Clarke (willy@clarke.no)
 * Version  : 0.0.1
 * Copyright: W. Clarke
 * License  : MIT
 * Descripti: This file contains the code that should be common to the
 *          : various applications in this project.
 ******************************************************************************/
#include "datastructures.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>  // for std::setprecision
#include <iostream>

// ---
// NOTE: Stream operator
// ---
std::ostream &operator<<(std::ostream &stream, const ww::tup &T)
{
   // ---
   // NOTE: The width need to be big enough to hold a negative sign.
   // ---
   size_t const P{5};
   size_t const W{P + 5};
   stream << ((T.W != 0) ? "Point:" : "Vector:");
   stream << " " << std::fixed << std::setprecision(P) << std::setw(W) << T.X  //<!
          << " " << std::fixed << std::setprecision(P) << std::setw(W) << T.Y  //<!
          << " " << std::fixed << std::setprecision(P) << std::setw(W) << T.Z  //<!
          << " " << std::fixed << std::setprecision(P) << std::setw(W) << T.W;
   return stream;
}

namespace ww
{
//------------------------------------------------------------------------------
tup Add(tup const &A, tup const &B)
{
   tup const Result{A.X + B.X, A.Y + B.Y, A.Z + B.Z, A.W + B.W};
   return (Result);
}

//------------------------------------------------------------------------------
tup Color(float const R, float const G, float const B)
{
   tup Result{R, G, B, 0.f};
   return (Result);
}

//------------------------------------------------------------------------------
tup Cross(tup const &A, tup const &B)
{
   tup const Result = Vector(A.Y * B.Z - A.Z * B.Y, A.Z * B.X - A.X * B.Z, A.X * B.Y - A.Y * B.X);
   return (Result);
}

//------------------------------------------------------------------------------
float Dot(tup const &A, tup const &B)
{
   float const Result = A.X * B.X +  //!<
                        A.Y * B.Y +  //<!
                        A.Y * B.Y +  //<!
                        A.Y * B.Y;   //<!
   return (Result);
}

//------------------------------------------------------------------------------
bool Equal(float const A, float const B)
{
   if (std::fabs(A - B) < EPSILON)
   {
      return true;
   }
   return false;
}

//------------------------------------------------------------------------------
bool Equal(tup const &A, tup const &B)
{
   bool const Result = Equal(A.X, B.X) &&  //<!
                       Equal(A.Y, B.Y) &&  //<!
                       Equal(A.Z, B.Z) &&  //<!
                       Equal(A.W, B.W);
   return (Result);
}

//------------------------------------------------------------------------------
bool IsPoint(tup const &Tup)
{
   bool Result{};
   Result = (Tup.W != 0);
   return (Result);
}

//------------------------------------------------------------------------------
bool IsVector(tup const &Tup)
{
   bool Result{};
   Result = !(Tup.W != 0);
   return (Result);
}

//------------------------------------------------------------------------------
tup Multiply(float const S, tup const &Tup)
{
   tup const Result{S * Tup.X, S * Tup.Y, S * Tup.Z, S * Tup.W};
   return (Result);
}

//------------------------------------------------------------------------------
tup Multiply(tup const A, tup const B)
{
   tup const Result{A.R * B.R, A.G * B.G, A.B * B.B, 0.f};
   return (Result);
}

//------------------------------------------------------------------------------
float MagSquared(tup const &Tup)
{
   float const Result = Tup.X * Tup.X +  //<!
                        Tup.Y * Tup.Y +  //<!
                        Tup.Z * Tup.Z +  //<!
                        Tup.W * Tup.W;   //<!
   return (Result);
}

//------------------------------------------------------------------------------
float Mag(tup const &Tup) { return (std::sqrt(MagSquared(Tup))); }

//------------------------------------------------------------------------------
tup Negate(tup const &Tup)
{
   tup const Result{-Tup.X, -Tup.Y, -Tup.Z, -Tup.W};
   return (Result);
}

//------------------------------------------------------------------------------
tup Normal(tup const &Tup)
{
   tup const Result = Tup / Mag(Tup);
   return (Result);
}

//------------------------------------------------------------------------------
tup Point(float A, float B, float C)
{
   tup Result{A, B, C, 1.f};
   return (Result);
}

//------------------------------------------------------------------------------
tup Sub(tup const &A, tup const &B)
{
   tup const Result = {A.X - B.X, A.Y - B.Y, A.Z - B.Z, A.W - B.W};
   return (Result);
}

//------------------------------------------------------------------------------
tup Vector(float A, float B, float C)
{
   tup Result{A, B, C, 0.f};
   return (Result);
}

// ---
// NOTE: Canvas methods/functions.
// ---
void WritePixel(canvas &Canvas, int X, int Y, tup const &Color)
{
   Assert(Canvas.vXY.size() > (X + Y * Canvas.W), __FILE__, __LINE__);
   Canvas.vXY[X + Y * Canvas.W] = Color;
}
//------------------------------------------------------------------------------
tup ReadPixel(canvas &Canvas, int X, int Y)
{
   Assert(Canvas.vXY.size() > (X + Y * Canvas.W), __FILE__, __LINE__);
   tup const Result = Canvas.vXY[X + Y * Canvas.W];
   return (Result);
}

//------------------------------------------------------------------------------
// ---
// NOTE: The Portable Pix Map header.
//     : Returns a string stream with a magic number, the dimension and the pixel
//     : resolution.
// ---
std::string PPMHeader(canvas const &Canvas)
{
   std::strstream ss;
   ss << "P3\n" << Canvas.W << " " << Canvas.H << "\n255\n";
   return (ss.str());
}

//------------------------------------------------------------------------------
// NOTE: Write the canvas to a Portable Pix Map file.
void WriteToPPM(canvas const &Canvas, std::string const &Filename)
{
   std::ofstream O(Filename, std::ofstream::out);
   O << PPMHeader(Canvas) << std::flush;
   int PixelCount{};
   for (size_t Idx = 0;           //<!
        Idx < Canvas.vXY.size();  //<!
        ++Idx, ++PixelCount)
   {
      // NOTE: As per the standard the count per line should not exceed 70.
      //     : So 33 pixels x 3 byte per pixel should do the trick.
      if ((!(PixelCount % 33) || (!(Idx % Canvas.W))) && (Idx > 0))
      {
         O << "\n" << std::flush;
         PixelCount = 0;
      }

      // NOTE: Truncate the float values between 0.f and 1.f.
      // clang-format off
      O <<        int(255 * std::max<float>(0.f, std::min<float>(1.f, Canvas.vXY[Idx].R)))  //<!
        << " " << int(255 * std::max<float>(0.f, std::min<float>(1.f, Canvas.vXY[Idx].G)))  //<!
        << " " << int(255 * std::max<float>(0.f, std::min<float>(1.f, Canvas.vXY[Idx].B)))  //<!
        << " ";
      // clang-format on
   }
   O << "\n" << std::flush;
   O.close();
}
};  // namespace ww

// ---
// NOTE: The Negate operator.
// ---
ww::tup operator+(ww::tup const &A, ww::tup const &B) { return (ww::Add(A, B)); }
ww::tup operator-(ww::tup const &Tup) { return (ww::Negate(Tup)); }
ww::tup operator-(ww::tup const &A, ww::tup const &B) { return (ww::Sub(A, B)); }
ww::tup operator*(float const S, ww::tup const &Tup) { return (ww::Multiply(S, Tup)); }
ww::tup operator*(ww::tup const &Tup, float const S) { return (ww::Multiply(S, Tup)); }
ww::tup operator*(ww::tup const &A, ww::tup const &B) { return (ww::Multiply(A, B)); }
// ---
// NOTE: Division operator does not check for divide by zero; Who cares?
// ---
ww::tup operator/(ww::tup const &Tup, float const S) { return (ww::Multiply(1.f / S, Tup)); }
