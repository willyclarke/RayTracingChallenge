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

std::ostream &operator<<(std::ostream &stream, const ww::matrix &M)
{
  size_t const P{5};
  size_t const W{P + 5};
  switch (M.Dimension)
  {
    case 2:
      stream << "\n |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[0].X         //<!
             << " |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[0].Y << " |";  //<!
      stream << "\n |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[1].X         //<!
             << " |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[1].Y << " |";  //<!
      break;
    case 3:
      stream << "\n |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[0].X         //<!
             << " |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[0].Y           //<!
             << " |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[0].Z << " |";  //<!
      stream << "\n |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[1].X         //<!
             << " |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[1].Y           //<!
             << " |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[1].Z << " |";  //<!
      stream << "\n |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[2].X         //<!
             << " |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[2].Y           //<!
             << " |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[2].Z << " |";  //<!
      break;
    default:
      stream << "\n |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[0].X         //<!
             << " |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[0].Y           //<!
             << " |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[0].Z           //<!
             << " |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[0].W << " |";  //<!
      stream << "\n |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[1].X         //<!
             << " |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[1].Y           //<!
             << " |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[1].Z           //<!
             << " |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[1].W << " |";  //<!
      stream << "\n |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[2].X         //<!
             << " |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[2].Y           //<!
             << " |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[2].Z           //<!
             << " |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[2].W << " |";  //<!
      stream << "\n |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[3].X         //<!
             << " |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[3].Y           //<!
             << " |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[3].Z           //<!
             << " |" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[3].W << " |";  //<!
      break;
  };
  stream << "\n";
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
tup Mul(float const S, tup const &Tup)
{
  tup const Result{S * Tup.X, S * Tup.Y, S * Tup.Z, S * Tup.W};
  return (Result);
}

//------------------------------------------------------------------------------
tup Mul(tup const A, tup const B)
{
  tup const Result{A.R * B.R, A.G * B.G, A.B * B.B, A.W * B.W};
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
// NOTE: This is actually take two of creating a PPM file and it uses a slightly different
//       approach compared to the file writing functin below.
//       For now I will keep the code around for reference.
//------------------------------------------------------------------------------
int WriteToPPMFile(canvas const &Canvas, std::string const &Filename)
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

    // NOTE: So we have a string with the data, lets write to the file.
    for (size_t I = 0;       //<!
         I < Output.size();  //<!
         ++I)
    {
      std::putc(Output[I], fp);
    }
  }
  std::putc('\n', fp);
  std::fclose(fp);

  // NOTE: move the file to the name it is supposed to have.
  int const Result = std::rename(Tmp.c_str(), Filename.c_str());
  return (Result);
}

//------------------------------------------------------------------------------
// NOTE: Write the canvas to a Portable Pix Map file.
void WriteToPPM(canvas const &Canvas, std::string const &Filename)
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
    if ((!(PixelCount % 33) || (!(Idx % Canvas.W))) && (Idx > 0))
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

// ---
// NOTE: Matrix functions.
// ---
matrix Matrix44()
{
  matrix Result{};
  return (Result);
}

//------------------------------------------------------------------------------
matrix Matrix44(tup const &R0, tup const &R1, tup const &R2, tup const &R3)
{
  matrix Result{R0, R1, R2, R3};
  return (Result);
}

//------------------------------------------------------------------------------
// NOTE: Use the data structure for 4x4 matrix for all types of matrixes.
matrix Matrix33(tup const &R0, tup const &R1, tup const &R2)
{
  matrix M{R0, R1, R2, tup{}};
  M.Dimension = 3;
  return (M);
}

//------------------------------------------------------------------------------
// NOTE: Use the data structure for 4x4 matrix for all types of matrixes.
matrix Matrix22(tup const &R0, tup const &R1)
{
  matrix M{R0, R1, tup{}, tup{}};
  M.Dimension = 2;
  return (M);
}

//------------------------------------------------------------------------------
float Get(matrix const &M, int Row, int Col)
{
  Assert(Row < M.Dimension, __FUNCTION__, __LINE__);
  Assert(Col < M.Dimension, __FUNCTION__, __LINE__);
  return M.R[Row].C[Col];
};

void Set(matrix &M, int Row, int Col, float Value) { M.R[Row].C[Col] = Value; }

//------------------------------------------------------------------------------
matrix Mul(matrix const &A, matrix const &B)
{
  matrix M{};
  for (size_t Row = 0;  ///<!
       Row < 4;         ///<!
       ++Row)
  {
    for (size_t Col = 0;  ///<!
         Col < 4;         ///<!
         ++Col)
    {
      float const Mrc = Get(A, Row, 0) * Get(B, 0, Col) +  //
                        Get(A, Row, 1) * Get(B, 1, Col) +  //
                        Get(A, Row, 2) * Get(B, 2, Col) +  //
                        Get(A, Row, 3) * Get(B, 3, Col);   //
      Set(M, Row, Col, Mrc);
    }
  }
  return (M);
}

//------------------------------------------------------------------------------
tup Mul(matrix const &A, tup const &T)
{
  tup Result{};
  for (size_t Idx = 0;  ///<!
       Idx < 4;         ///<!
       ++Idx)
  {
    Result.C[Idx] = Get(A, Idx, 0) * T.C[0] +  //
                    Get(A, Idx, 1) * T.C[1] +  //
                    Get(A, Idx, 2) * T.C[2] +  //
                    Get(A, Idx, 3) * T.C[3];
  }
  return (Result);
}

matrix I()
{
  matrix Identity{};
  Set(Identity, 0, 0, 1.f);
  Set(Identity, 1, 1, 1.f);
  Set(Identity, 2, 2, 1.f);
  Set(Identity, 3, 3, 1.f);
  return (Identity);
}

//------------------------------------------------------------------------------
bool Equal(matrix const &A, matrix const &B)
{
  return Equal(A.R[0], B.R[0]) &&  //
         Equal(A.R[1], B.R[1]) &&  //
         Equal(A.R[2], B.R[2]) &&  //
         Equal(A.R[3], B.R[3]);
}

//------------------------------------------------------------------------------
matrix Transpose(matrix const &M)
{
  matrix R{};
  for (size_t Row = 0;  ///<!
       Row < 4;         ///<!
       ++Row)
  {
    for (size_t Col = 0;  ///<!
         Col < 4;         ///<!
         ++Col)
    {
      R.R[Row].C[Col] = M.R[Col].C[Row];
    }
  }
  return (R);
}

//------------------------------------------------------------------------------
float Determinant22(matrix const &M)
{
  // NOTE: The determinant is D = a*d - b*c
  float const Result = M.R[0].C[0] * M.R[1].C[1] - M.R[1].C[0] * M.R[0].C[1];
  return (Result);
}

//------------------------------------------------------------------------------
matrix SubMatrix(matrix const &M, int RemoveRow, int RemoveCol)
{
  matrix R{};
  int ShiftR{};

  // NOTE: For each row copy source until we get to the removerow
  //       When we get to the remove row the source shifts by one,
  //       so that we skip the RemoveRow.
  for (int Row = 0;            ///<!
       Row < M.Dimension - 1;  ///<!
       ++Row)
  {
    if (Row == RemoveRow) ShiftR++;
    R.R[Row] = M.R[Row + ShiftR];
  }

  // NOTE: Dimension reduces by 1 when we remove one row/column.
  R.Dimension = M.Dimension - 1;

  // NOTE: Copy the columns.
  for (int Row = 0;        ///<!
       Row < R.Dimension;  ///<! The rows are already in order, so we only need to
       ++Row)              ///<! iterate over the result rows.
  {
    for (int Col = 0;        ///<!
         Col < M.Dimension;  ///<!
         ++Col)
    {
      if (Col < RemoveCol) continue;  // No point in copy of data that is already correctly placed.

      if (Col < M.Dimension - 1)
        R.R[Row].C[Col] = R.R[Row].C[Col + 1];
      else
        R.R[Row].C[Col] = 0.f;
    }
  }

  return (R);
}
};  // namespace ww

// ---
// NOTE: The Negate operator.
// ---
ww::tup operator+(ww::tup const &A, ww::tup const &B) { return (ww::Add(A, B)); }
ww::tup operator-(ww::tup const &Tup) { return (ww::Negate(Tup)); }
ww::tup operator-(ww::tup const &A, ww::tup const &B) { return (ww::Sub(A, B)); }
ww::tup operator*(float const S, ww::tup const &Tup) { return (ww::Mul(S, Tup)); }
ww::tup operator*(ww::tup const &Tup, float const S) { return (ww::Mul(S, Tup)); }
ww::tup operator*(ww::tup const &A, ww::tup const &B) { return (ww::Mul(A, B)); }
ww::matrix operator*(ww::matrix const &A, ww::matrix const &B) { return (ww::Mul(A, B)); }
ww::tup operator*(ww::matrix const &A, ww::tup const &T) { return (ww::Mul(A, T)); }
// ---
// NOTE: Division operator does not check for divide by zero; Who cares?
// ---
ww::tup operator/(ww::tup const &Tup, float const S) { return (ww::Mul(1.f / S, Tup)); }
