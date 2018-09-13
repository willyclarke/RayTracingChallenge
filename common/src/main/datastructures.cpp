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

#include <algorithm>  // for std::copy
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>  // for std::setprecision
#include <iostream>
#include <iterator>
#include <memory>  // for shared pointer.
#include <sstream>
#include <string>
#include <vector>
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
   // std::cout << "X:" << X << ". Y:" << Y << ". clc:" << (X + Y * Canvas.W) << ". size:" << Canvas.vXY.size()
   //          << std::endl;

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
//     : Returns a string with a magic number, the dimension and the pixel resolution.
// ---
std::string PPMHeader(canvas const &Canvas)
{
   std::string const Result = "P3\n" + std::to_string(Canvas.W) + " " + std::to_string(Canvas.H) + "\n255";
   return (Result);
}

//------------------------------------------------------------------------------
int WriteToPPM(canvas const &Canvas, std::string const &Filename)
{
   std::string Tmp{Filename};
   if (!Tmp.size())
   {
      // NOTE: Create a temporary filename.
      static int FilenameModifier{};
      wchar_t const FilenameM = L"ABCDEFG=#"[FilenameModifier];
      FilenameModifier = (FilenameModifier + 1) % 7;
      std::stringstream ss;
      ss << FilenameM;
      Tmp = Filename + ss.str();
   }

   FILE *fp = std::fopen(Tmp.c_str(), "w");
   std::string const Header = PPMHeader(Canvas);
   for (size_t Idx = 0;       //<!
        Idx < Header.size();  //<!
        ++Idx)
   {
      std::putc(Header[Idx], fp);
   }
   // NOTE: Body need to come on a separate line, so we do a new line here.
   std::putc('\n', fp);

   Assert(Canvas.vXY.size() == Canvas.W * Canvas.H, __FILE__, __LINE__);

   // char const *PPM =
   //    "\n255 0 153 255 0 153 255 0 153 255 0 153 255 0 153 255 0 153 255 0 153 255 0 153 255 0 153 255 0 "
   //    "153 255 0 153 255 0 153 255 0 153 255 0 153 255 0 153\n255 0 153 255 0 153 255 0 153 255 0 153 255 0 153 255 "
   //    "0 153 255 0 153 255 0 153 255 0 153 255 0 153 255 0 153 255 0 153 255 0 153 255 0 153 255 0 153\n";

   // std::cout << ". SizeofBody:" << std::strlen(PPM) << ". SizeofHeader:" << Header.size() << std::endl;
   // std::cout << Header << PPM << std::endl;

   // for (size_t Idx = 0;          ///<!
   //     Idx < std::strlen(PPM);  ///<!
   //     ++Idx)
   //{
   //   std::putc(PPM[Idx], fp);
   //}
#if 1
   int PixelCount{};
   for (size_t Idx = 0;           //<!
        Idx < Canvas.vXY.size();  //<!
        ++Idx, ++PixelCount)
   {
      // NOTE: As per the standard the count per line should not exceed 70.
      //     : So 33 pixels x 3 byte per pixel should do the trick.
      if ((!(PixelCount % 15) || (!(Idx % Canvas.W))) && (Idx > 0))
      {
         std::putc('\n', fp);
         PixelCount = 0;
      }

      // NOTE: Truncate the float values between 0.f and 1.f.
      // clang-format off
      int const cR = int(255 * std::max<float>(0.f, std::min<float>(1.f, Canvas.vXY[Idx].R)));
      int const cG = int(255 * std::max<float>(0.f, std::min<float>(1.f, Canvas.vXY[Idx].G)));
      int const cB = int(255 * std::max<float>(0.f, std::min<float>(1.f, Canvas.vXY[Idx].B)));
      // clang-format on
      std::string const Output = std::to_string(cR) + " " + std::to_string(cG) + " " + std::to_string(cB) + " ";

      // NOTE: So we have a stringstream with the data, lets write to the file.
      for (size_t I = 0;       //<!
           I < Output.size();  //<!
           ++I)
      {
         std::putc(Output[I], fp);
      }
   }
   std::putc('\n', fp);
#endif
   std::fclose(fp);

   std::cout << std::endl;
   std::cout << "\nFinished with:" << Tmp << std::endl;

   // ---
   // NOTE: Reopen the file to get it flushed. Not sure if this works, but its worth a try.
   // ---
   //{
   //   FILE *fp = std::fopen(Tmp.c_str(), "r");
   //   std::fclose(fp);
   //}
   // NOTE: move the file to the name it is supposed to have.
   int const Result = std::rename(Tmp.c_str(), Filename.c_str());
   return (Result);
}

//------------------------------------------------------------------------------
// NOTE: Write the canvas to a Portable Pix Map file.
void WriteToPPMFile(canvas const &Canvas, std::string const &Filename)
{
   std::ofstream O(Filename, std::ofstream::out | std::ofstream::trunc);

   O << PPMHeader(Canvas) << std::endl;

   int PixelCount{};
   for (size_t Idx = 0;           //<!
        Idx < Canvas.vXY.size();  //<!
        ++Idx, ++PixelCount)
   {
      std::string Output{};
      // NOTE: As per the standard the count per line should not exceed 70.
      //     : So 33 pixels x 3 byte per pixel should do the trick.
      if ((!(PixelCount % 15) || (!(Idx % Canvas.W))) && (Idx > 0))
      {
         Output += "\n";
         PixelCount = 0;
      }

      // NOTE: Truncate the float values between 0.f and 1.f.
      // clang-format off
      Output += std::to_string(int(255 * std::max<float>(0.f, std::min<float>(1.f, Canvas.vXY[Idx].R))))  //!<
         + " " + std::to_string(int(255 * std::max<float>(0.f, std::min<float>(1.f, Canvas.vXY[Idx].G))))  //<!
         + " " + std::to_string(int(255 * std::max<float>(0.f, std::min<float>(1.f, Canvas.vXY[Idx].B))))  //<!
         + " ";
      // clang-format on
      O << Output;

      Assert(Output.size() < 12, __FUNCTION__, __LINE__);
   }
   O << std::endl;
   O.close();
}

//------------------------------------------------------------------------------
// NOTE: Read a PPM file from disk.
// Return: A shared pointer of the allocated canvas.
// TODO: Fault checking and robustness.
//------------------------------------------------------------------------------
std::shared_ptr<canvas> ReadFromPPM(std::string const &Filename)
{
   std::shared_ptr<canvas> Result{};
   std::ifstream I(Filename, std::ifstream::in);
   if (I.is_open())
   {
      // std::istream_iterator<std::string> begin(I);
      // std::istream_iterator<std::string> end{};
      std::vector<std::string> vstrings{};  //(begin, end);
      while (!I.eof() && I.good())
      {
         std::string myString{};
         I >> myString;
         if (myString.size())
         {
            vstrings.push_back(myString);
         }
      }
      // std::cout << "Size of vstrings:" << vstrings.size() << ". Line:" << __LINE__ << std::endl;

      // NOTE: Check to see if magic number is correct.
      if (vstrings.size() && (vstrings[0] == "P3"))
      {
         if (vstrings.size() > 3)
         {
            int const W = std::stoi(vstrings[1]);    // Width
            int const H = std::stoi(vstrings[2]);    // Height
            float const L = std::stof(vstrings[3]);  // How many levels, aka resolution

            // NOTE: Allocate on the heap, not stack.
            std::shared_ptr<canvas> ptrCanvas = std::make_shared<canvas>(W, H);
            canvas &Cv = *ptrCanvas;
            Assert(Cv.W == W && Cv.H == H, __FUNCTION__, __LINE__);

            if ((Cv.W == W) && (Cv.H == H) && (L > 0))
            {
               // NOTE: The rest of the vector consist of the individual colors, 3 bytes per pixel.
               int CvX{};
               int CvY{};
               for (size_t Idx = 4;         //<!
                    Idx < vstrings.size();  //<!
                    Idx += 3)
               {
                  float const R = std::stof(vstrings[Idx]) / L;
                  float const G = std::stof(vstrings[Idx + 1]) / L;
                  float const B = std::stof(vstrings[Idx + 2]) / L;

                  WritePixel(Cv, CvX, CvY, Color(R, G, B));

                  // Calculate new position
                  CvX++;
                  CvX %= W;
                  if (CvX == 0) ++CvY;
               }
               // NOTE: Update the return pointer.
               Result = ptrCanvas;
            }
            else
            {
               std::cerr << __FUNCTION__ << ": Loading of " << Filename << " failed. Resolution:"  //<!
                         << L << ". W:" << W << ". H:" << H << "." << std::endl;
            }
         }
      }
   }
   return (Result);
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
