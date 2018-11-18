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
      stream << "\n {" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[0].X         //<!
             << " ," << std::fixed << std::setprecision(P) << std::setw(W) << M.R[0].Y           //<!
             << " ," << std::fixed << std::setprecision(P) << std::setw(W) << M.R[0].Z           //<!
             << " ," << std::fixed << std::setprecision(P) << std::setw(W) << M.R[0].W << " }";  //<!
      stream << "\n {" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[1].X         //<!
             << " ," << std::fixed << std::setprecision(P) << std::setw(W) << M.R[1].Y           //<!
             << " ," << std::fixed << std::setprecision(P) << std::setw(W) << M.R[1].Z           //<!
             << " ," << std::fixed << std::setprecision(P) << std::setw(W) << M.R[1].W << " }";  //<!
      stream << "\n {" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[2].X         //<!
             << " ," << std::fixed << std::setprecision(P) << std::setw(W) << M.R[2].Y           //<!
             << " ," << std::fixed << std::setprecision(P) << std::setw(W) << M.R[2].Z           //<!
             << " ," << std::fixed << std::setprecision(P) << std::setw(W) << M.R[2].W << " }";  //<!
      stream << "\n {" << std::fixed << std::setprecision(P) << std::setw(W) << M.R[3].X         //<!
             << " ," << std::fixed << std::setprecision(P) << std::setw(W) << M.R[3].Y           //<!
             << " ," << std::fixed << std::setprecision(P) << std::setw(W) << M.R[3].Z           //<!
             << " ," << std::fixed << std::setprecision(P) << std::setw(W) << M.R[3].W << " }";  //<!
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
                       A.Z * B.Z +  //<!
                       A.W * B.W;   //<!
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
tup Normalize(tup const &Tup)
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

//------------------------------------------------------------------------------
float Radians(float Deg) { return (PI_F * Deg / 180.f); }
// ---
// NOTE: Canvas methods/functions.
// ---
void WritePixel(canvas &Canvas, int X, int Y, tup const &Color)
{
  // std::cout << "X:" << X                       //!<
  //          << ". Y:" << Y                     //!<
  //          << ". clc:" << (X + Y * Canvas.W)  //!<
  //          << ". size:" << Canvas.vXY.size()  //!<
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
is_invertible_return IsInvertible(matrix const &M)
{
  // NOTE: When the determinant has already been calculated we just return that result.
  //       Otherwise the Determinant is calculated and a tuple is returned..
  is_invertible_return Result{M.ID};

  if (!Result.IsInvertible)
  {
    Result.Determinant = Determinant(M);
    Result.IsInvertible = Result.Determinant != 0;
  }
  return (Result);
}

//------------------------------------------------------------------------------
matrix Inverse(matrix const &M)
{
  // NOTE: Inverse of matrix is done by
  // 1. Calculate the determinant. If different than zero ok
  // 2. Calculate the cofactor of the input matrix.
  // 3. Transpose the Resulting matrix.
  // 4. Divide each element by the determinant.

  matrix Result{};
  float &DetM = Result.ID.Determinant;
  DetM = Determinant(M);

  if (Equal(DetM, 0.f)) return (Result);

  // NOTE: Since we did not return above the matrix must be invertible
  Result.ID.IsInvertible = true;

  for (int Row = 0;        ///<!
       Row < M.Dimension;  ///<!
       ++Row)
  {
    for (int Col = 0;        ///<!
         Col < M.Dimension;  ///<!
         ++Col)
    {
      Result.R[Row].C[Col] = Cofactor(M, Row, Col) / DetM;
    }
  }

  // NOTE: The transposed matrix is not updated with the Determinant and the IsInvertible flag.
  Result = Transpose(Result);

  return (Result);
}
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

//------------------------------------------------------------------------------
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
float Determinant33(matrix const &M)
{
  float const CF0 = Cofactor33(M, 0, 0);
  float const CF1 = Cofactor33(M, 0, 1);
  float const CF2 = Cofactor33(M, 0, 2);

  float const Result = M.R0.C[0] * CF0 + M.R0.C[1] * CF1 + M.R0.C[2] * CF2;
  return (Result);
}

//------------------------------------------------------------------------------
float Determinant44(matrix const &M)
{
  matrix const M0 = SubMatrix(M, 0, 0);
  matrix const M1 = SubMatrix(M, 0, 1);
  matrix const M2 = SubMatrix(M, 0, 2);
  matrix const M3 = SubMatrix(M, 0, 3);
  float const DetM0 = Determinant33(M0);
  float const DetM1 = Determinant33(M1);
  float const DetM2 = Determinant33(M2);
  float const DetM3 = Determinant33(M3);
  float const Result = M.R0.C[0] * DetM0 - M.R0.C[1] * DetM1 + M.R0.C[2] * DetM2 - M.R0.C[3] * DetM3;
  return (Result);
}

//------------------------------------------------------------------------------
float Determinant(matrix const &M)
{
  float Result{};
  if (M.Dimension == 4)
    Result = Determinant44(M);
  else if (M.Dimension == 3)
    Result = Determinant33(M);
  else if (M.Dimension == 2)
    Result = Determinant22(M);
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

//------------------------------------------------------------------------------
float Minor(matrix const &M, int RemoveRow, int RemoveCol)
{
  matrix const SM = SubMatrix(M, RemoveRow, RemoveCol);
  float Result = Determinant22(SM);

  return (Result);
}

//------------------------------------------------------------------------------
float Cofactor33(matrix const &M, int RemoveRow, int RemoveCol)
{
  Assert(M.Dimension == 3, __FUNCTION__, __LINE__);
  // NOTE: Change sign for the Cofactor when the sum of Row and Col is an odd number.
  //       So; move to -2 for sign and then add 1.
  int const Sign = -((RemoveRow + RemoveCol) % 2) * 2 + 1;
  float const Result = Sign * Minor(M, RemoveRow, RemoveCol);
  return (Result);
}

//------------------------------------------------------------------------------
float Cofactor44(matrix const &M, int RemoveRow, int RemoveCol)
{
  Assert(M.Dimension == 4, __FUNCTION__, __LINE__);
  matrix const A = SubMatrix(M, RemoveRow, RemoveCol);

  int const Sign = -((RemoveRow + RemoveCol) % 2) * 2 + 1;
  float const Result = Sign * Determinant33(A);
  return (Result);
}

//------------------------------------------------------------------------------
float Cofactor(matrix const &M, int RemoveRow, int RemoveCol)
{
  float Result{};

  if (M.Dimension == 4)
    Result = Cofactor44(M, RemoveRow, RemoveCol);
  else if (M.Dimension == 3)
    Result = Cofactor33(M, RemoveRow, RemoveCol);

  return (Result);
}

matrix Translation(float X, float Y, float Z)
{
  matrix M{I()};
  Set(M, 0, 3, X);
  Set(M, 1, 3, Y);
  Set(M, 2, 3, Z);

  return (M);
}

//------------------------------------------------------------------------------
matrix Scale(float X, float Y, float Z)
{
  matrix M{I()};

  Set(M, 0, 0, X);
  Set(M, 1, 1, Y);
  Set(M, 2, 2, Z);

  return (M);
}

//------------------------------------------------------------------------------
matrix RotateX(float Alfa)
{
  matrix M{I()};
  Set(M, 1, 1, std::cos(Alfa));
  Set(M, 1, 2, -std::sin(Alfa));
  Set(M, 2, 1, std::sin(Alfa));
  Set(M, 2, 2, std::cos(Alfa));

  return (M);
}

//------------------------------------------------------------------------------
matrix RotateY(float Alfa)
{
  matrix M{I()};
  Set(M, 0, 0, std::cos(Alfa));
  Set(M, 0, 2, std::sin(Alfa));
  Set(M, 2, 0, -std::sin(Alfa));
  Set(M, 2, 2, std::cos(Alfa));

  return (M);
}

//------------------------------------------------------------------------------
matrix RotateZ(float Alfa)
{
  matrix M{I()};
  Set(M, 0, 0, std::cos(Alfa));
  Set(M, 0, 1, -std::sin(Alfa));
  Set(M, 1, 0, std::sin(Alfa));
  Set(M, 1, 1, std::cos(Alfa));

  return (M);
}

//------------------------------------------------------------------------------
matrix Shearing(float Xy, float Xz, float Yx, float Yz, float Zx, float Zy)
{
  matrix M{I()};
  Set(M, 0, 1, Xy);
  Set(M, 0, 2, Xz);
  Set(M, 1, 0, Yx);
  Set(M, 1, 2, Yz);
  Set(M, 2, 0, Zx);
  Set(M, 2, 1, Zy);
  return (M);
}

//------------------------------------------------------------------------------
matrix TranslateScaleRotate(                   //!<
    float TransX, float TransY, float TransZ,  //!< Translation is in m(?)
    float ScaleX, float ScaleY, float ScaleZ,  //!< Scale input is unitless.
    float AlfaX, float AlfaY, float AlfaZ      //!< Input rotation in radians.
)
{
  matrix const M = Translation(TransX, TransY, TransZ) *             //!<
                   Scale(ScaleX, ScaleY, ScaleZ) *                   //!<
                   RotateX(AlfaX) * RotateY(AlfaY) * RotateZ(AlfaZ)  //!<
      ;

  return (M);
}

/// ---
/// \fn Sphere releated functions
/// ---
//------------------------------------------------------------------------------
sphere Sphere(tup const &Center, float Radius)
{
  sphere S{};
  S.Center = Center;
  S.R = Radius;
  return (S);
}

//------------------------------------------------------------------------------
intersect Intersect(sphere const &Sphere, ray &Ray)
{
  intersect Result{};
  return (Result);
}

/// ---
/// Ray releated functions.
/// ---
//------------------------------------------------------------------------------
ray Ray(tup const &O, tup const &D)
{
  Assert(IsPoint(O), __FUNCTION__, __LINE__);
  Assert(IsVector(D), __FUNCTION__, __LINE__);

  ray Result{};
  Result.O = O;
  Result.D = D;

  return (Result);
}

//------------------------------------------------------------------------------
int Count(intersect const &I) { return I.Count(); }

//------------------------------------------------------------------------------
tup Position(ray const &R, float t)
{
  tup Result{};

  Result = R.O + R.D * t;
  Assert(IsPoint(Result), __FUNCTION__, __LINE__);

  return (Result);
}

// intersection Intersection(float t, object *Object);
// intersections Intersections(intersection const &I1, intersection const &I2);
// intersect Intersect(object *pObject, ray const &Ray);
//------------------------------------------------------------------------------
intersect IntersectSphere(sphere const &Sphere, ray const &Ray)
{
  intersect Result{};
  // NOTE: See explanation from:
  // https://stackoverflow.com/questions/1073336/circle-line-segment-collision-detection-algorithm#1084899
  //
  // NOTE: The vector from the sphere's center to the ray origin
  //       Remember that the sphere is centered at the world origin
  tup const Sphere2Ray = Ray.O - Point(0.f, 0.f, 0.f);
  // std::cout << "Sphere2Ray : " << Sphere2Ray << std::endl;

  float const A = Dot(Ray.D, Ray.D);
  float const B = 2 * Dot(Ray.D, Sphere2Ray);
  float const C = Dot(Sphere2Ray, Sphere2Ray) - 1.f;
  float const Discriminant = B * B - 4 * A * C;
  // std::cout << "A:" << A << ". B:" << B << ". C:" << C << ". Discriminant:" << Discriminant << std::endl;

  if (Discriminant >= 0)
  {
    float const t1 = (-B - std::sqrt(Discriminant)) / (2 * A);
    float const t2 = (-B + std::sqrt(Discriminant)) / (2 * A);
    intersection I{};
    I.Object = Sphere;
    if (t1 > t2)
    {
      I.t = t2;
      Result.vI.push_back(I);
      I.t = t1;
      Result.vI.push_back(I);
    }
    else
    {
      I.t = t1;
      Result.vI.push_back(I);
      I.t = t2;
      Result.vI.push_back(I);
    }
  }
  return (Result);
}

//------------------------------------------------------------------------------
intersect Intersect(sphere const &Object, ray const &RayIn)
{
  intersect Result{};

  // ---
  // NOTE: The object to which we are trying to calculate the intersect may
  //       kind of not be placed at origin. So use its transform to 'move' the
  //       ray by calculation of the inverse.
  // ---
  ray const Ray = Transform(RayIn, Inverse(Object.T));

  // ---
  // NOTE: See explanation from:
  // https://stackoverflow.com/questions/1073336/circle-line-segment-collision-detection-algorithm#1084899
  //
  // NOTE: The vector from the sphere's center to the ray origin
  //       Remember that the sphere is centered at the world origin
  tup const Object2Ray = Ray.O - Point(0.f, 0.f, 0.f);
  // std::cout << "Sphere2Ray : " << Sphere2Ray << std::endl;

  float const A = Dot(Ray.D, Ray.D);
  float const B = 2 * Dot(Ray.D, Object2Ray);
  float const C = Dot(Object2Ray, Object2Ray) - 1.f;
  float const Discriminant = B * B - 4 * A * C;
  // std::cout << "A:" << A << ". B:" << B << ". C:" << C << ". Discriminant:" << Discriminant << std::endl;

  if (Discriminant >= 0)
  {
    float const t1 = (-B - std::sqrt(Discriminant)) / (2 * A);
    float const t2 = (-B + std::sqrt(Discriminant)) / (2 * A);
    intersection I{};
    I.Object = Object;
    if (t1 > t2)
    {
      I.t = t2;
      Result.vI.push_back(I);
      I.t = t1;
      Result.vI.push_back(I);
    }
    else
    {
      I.t = t1;
      Result.vI.push_back(I);
      I.t = t2;
      Result.vI.push_back(I);
    }
  }
  return (Result);
}

//------------------------------------------------------------------------------
intersection Intersection(float t, object *pObject)
{
  intersection Result{};
  Result.t = t;
  Result.pObject = pObject;

  return (Result);
}

//------------------------------------------------------------------------------
intersections Intersections(intersection const &I1, intersection const &I2)
{
  intersections Result{};
  Result.vI.push_back(I1);
  Result.vI.push_back(I2);

  return (Result);
}

//------------------------------------------------------------------------------
intersections &Intersections(intersections &XS, intersection const &I)
{
  XS.vI.push_back(I);
  return (XS);
}

//------------------------------------------------------------------------------
bool Equal(sphere const &A, sphere const &B) { return (Equal(A.Center, B.Center) && Equal(A.R, B.R)); }

//------------------------------------------------------------------------------
bool Equal(intersection const &A, intersection const &B)
{
  bool const Result = Equal(A.t, B.t);
  // std::cout << "Result:" << Result << ". Equal -> A:" << A.t << " B:" << B.t << std::endl;

  return (Result);
}

//------------------------------------------------------------------------------
intersection Hit1(intersections const &Intersections)
{
  intersection Result{};

  // ---
  // NOTE: Check that there is anything in the vector.
  // ---
  if (Intersections.Count() == 2)
  {
    float const &H0 = Intersections.vI[0].t;
    float const &H1 = Intersections.vI[1].t;
    // Decide on the minimum t. Seed it with the first entry of intersections.
    float Checkt = std::max<float>(H0, H1);

    // NOTE: There is at least one hit since the checked value is not negative.
    if (Checkt >= 0)
    {
      if ((H0 >= 0) && (H0 <= H1))
      {
        Result = Intersections.vI[0];
        // std::cout << "Count2A:Updated result with t:" << Result.t << std::endl;
      }
      else if ((H1 >= 0) && ((H1 <= H0) || (H0 < 0)))
      {
        Result = Intersections.vI[1];
        // std::cout << "Count2B:Updated result with t:" << Result.t << std::endl;
      }
    }
  }
  // NOTE: This part of the code is untested .... FIXME: (Willy Clarke)
  else if (Intersections.Count() == 1)
  {
    // NOTE: assert here until this code is ready for testing.
    Assert(Intersections.Count() == 0, __FUNCTION__, __LINE__);
    float const &H0 = Intersections.vI[0].t;
    if (H0 > 0)
    {
      Result = Intersections.vI[0];
      std::cout << "Count1:Updated result with t:" << Result.t << std::endl;
    }
  }
  return (Result);
}

//------------------------------------------------------------------------------
intersection Hit(intersections const &Intersections)
{
  intersection Result{};

  bool DoFirstUpdate{true};

  for (auto I : Intersections.vI)
  {
    if (I.t > 0)
    {
      if (DoFirstUpdate)
      {
        Result = I;
        DoFirstUpdate = false;
      }
      else if (I.t < Result.t)
      {
        Result = I;
      }
    }
  }
  return (Result);
}

//------------------------------------------------------------------------------
ray Mul(matrix const &M, ray const &R)
{
  ray Result{};
  Result.O = M * R.O;
  Result.D = M * R.D;
  return (Result);
}

//------------------------------------------------------------------------------
ray Transform(ray const &R, matrix const &M) { return M * R; }

//------------------------------------------------------------------------------
tup NormalAt(object const &O, tup const P)
{
  tup Result{};
  return (Result);
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
ww::ray operator*(ww::matrix const &M, ww::ray const &R) { return (ww::Mul(M, R)); }
// bool operator==(ww::sphere const &A, ww::sphere const &B) {return (ww::Equal(A, B));}
bool operator==(ww::sphere const &A, ww::sphere const &B) { return (ww::Equal(A, B)); }
bool operator==(ww::intersection const &A, ww::intersection const &B) { return (ww::Equal(A, B)); }
// ---
// NOTE: Division operator does not check for divide by zero; Who cares?
// ---
ww::tup operator/(ww::tup const &Tup, float const S) { return (ww::Mul(1.f / S, Tup)); }
