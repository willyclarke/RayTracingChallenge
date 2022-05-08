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

#include <algorithm>  // for std::copy, std::sort
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
  stream << ((T.W != 0) ? "Point :" : "Vector:");
  stream << " " << std::fixed << std::setprecision(P) << std::setw(W) << T.X  //<!
         << " " << std::fixed << std::setprecision(P) << std::setw(W) << T.Y  //<!
         << " " << std::fixed << std::setprecision(P) << std::setw(W) << T.Z  //<!
         << " " << std::fixed << std::setprecision(P) << std::setw(W) << T.W;
  return stream;
}

//------------------------------------------------------------------------------
std::ostream &operator<<(std::ostream &stream, const ww::ray &R)
{
  stream << "Origin   : " << R.Origin << "\n"
         << "Direction: " << R.Direction;
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

//------------------------------------------------------------------------------
std::ostream &operator<<(std::ostream &stream, const ww::material &M)
{
  stream << "\nMaterial\nColor:" << M.Color << "\nSpecular:" << M.Specular << " Ambient:" << M.Ambient
         << " Diffuse:" << M.Diffuse << " Shininess:" << M.Shininess << std::endl;
  return stream;
}

//------------------------------------------------------------------------------
std::ostream &operator<<(std::ostream &stream, const ww::sphere &S)
{
  stream << "\nSphere\nCenter:" << S.Center << "\nSphere Material:" << S.Material << "Sphere Radius:" << S.Radius
         << "\nSphere Transform:" << S.Transform << std::endl;
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
  // if (std::fabs(A - B) < 2.f * EPSILON)
  if (std::fabs(A - B) < 1.f * EPSILON)
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
tup PixelAt(canvas const &Canvas, int X, int Y)
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

    // Assert(Output.size() < 12, __FUNCTION__, __LINE__);
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
matrix Scaling(float X, float Y, float Z)
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
                   Scaling(ScaleX, ScaleY, ScaleZ) *                 //!<
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
  S.Radius = Radius;
  return (S);
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
  Result.Origin = O;
  Result.Direction = Normalize(D);

  return (Result);
}

//------------------------------------------------------------------------------
int Count(intersections const &I) { return I.Count(); }

//------------------------------------------------------------------------------
tup PositionAt(ray const &R, float t)
{
  tup Result{};

  Result = R.Origin + R.Direction * t;
  Assert(IsPoint(Result), __FUNCTION__, __LINE__);

  return (Result);
}

// intersection Intersection(float t, object *Object);
// intersections Intersections(intersection const &I1, intersection const &I2);
// intersect Intersect(object *pShape, ray const &Ray);
//------------------------------------------------------------------------------
intersections IntersectSphere(sphere const &Sphere, ray const &Ray)
{
  intersections Result{};
  // NOTE: See explanation from:
  // https://stackoverflow.com/questions/1073336/circle-line-segment-collision-detection-algorithm#1084899
  //
  // NOTE: The vector from the sphere's center to the ray origin
  //       Remember that the sphere is centered at the world origin
  tup const Sphere2Ray = Ray.Origin - Point(0.f, 0.f, 0.f);
  // std::cout << "Sphere2Ray : " << Sphere2Ray << std::endl;

  float const A = Dot(Ray.Direction, Ray.Direction);
  float const B = 2 * Dot(Ray.Direction, Sphere2Ray);
  float const C = Dot(Sphere2Ray, Sphere2Ray) - 1.f;
  float const Discriminant = B * B - 4 * A * C;
  // std::cout << "A:" << A << ". B:" << B << ". C:" << C << ". Discriminant:" << Discriminant << std::endl;

  if (Discriminant >= 0)
  {
    float const t1 = (-B - std::sqrt(Discriminant)) / (2 * A);
    float const t2 = (-B + std::sqrt(Discriminant)) / (2 * A);
    intersection I{};

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

/**
 * Check if a ray have a local intersect with a plane.
 * Param: PtrShape: A ww::plane is needed to do some actual checks for intersections.
 * Return: ww::intersections.
 **/
ww::intersections LocalIntersectPlane(shared_ptr_shape PtrShape, ray const &Ray)
{
  Assert(PtrShape->isA<plane>(), __FUNCTION__, __LINE__);

  ww::intersections Result{};

  if (!PtrShape->isA<ww::plane>()) return Result;

  if (std::abs(Ray.Direction.Y) < EPSILON) return Result;

  float const t = -Ray.Origin.Y / Ray.Direction.Y;

  Result.vI.push_back(Intersection(t, PtrShape));

  return Result;
}

/**
 * Check if a ray have a local intersect with a sphere.
 * Param: PtrShape: A ww::sphere is needed to do some actual checks for intersections.
 * Return: ww::intersections.
 **/
ww::intersections LocalIntersectSphere(shared_ptr_shape PtrShape, ray const &Ray)
{
  Assert(PtrShape->isA<sphere>(), __FUNCTION__, __LINE__);

  ww::intersections Result{};

  if (!PtrShape->isA<ww::sphere>()) return Result;

  // ---
  // NOTE: See explanation from:
  // https://stackoverflow.com/questions/1073336/circle-line-segment-collision-detection-algorithm#1084899
  //
  // NOTE: The vector from the sphere's center to the ray origin
  //       Remember that the sphere is centered at the world origin
  tup const Sphere2Ray = Ray.Origin - Point(0.f, 0.f, 0.f);

  float const A = ww::Dot(Ray.Direction, Ray.Direction);
  float const B = 2 * ww::Dot(Ray.Direction, Sphere2Ray);
  float const C = ww::Dot(Sphere2Ray, Sphere2Ray) - 1.f;
  float const Discriminant = B * B - 4 * A * C;

  if (Discriminant >= 0)
  {
    float const t1 = (-B - std::sqrt(Discriminant)) / (2 * A);
    float const t2 = (-B + std::sqrt(Discriminant)) / (2 * A);
    intersection I{};

    // NOTE: take a copy of the object for future reference.
    I.pShape = PtrShape;
    Assert(I.pShape, __FUNCTION__, __LINE__);

    if (t1 > t2)
    {
      I.t = t2;
      Assert(I.pShape, __FUNCTION__, __LINE__);
      Result.vI.push_back(I);
      I.t = t1;
      Assert(I.pShape, __FUNCTION__, __LINE__);
      Result.vI.push_back(I);
    }
    else
    {
      I.t = t1;
      Assert(I.pShape, __FUNCTION__, __LINE__);
      Result.vI.push_back(I);
      I.t = t2;
      Assert(I.pShape, __FUNCTION__, __LINE__);
      Result.vI.push_back(I);
    }
  }
  return (Result);
}

//------------------------------------------------------------------------------
/// \brief Generic Local Intersect
///
///        Looks up the type of shape and calls the specific
///        LocalIntersect for the shape
//------------------------------------------------------------------------------
intersections LocalIntersect(shared_ptr_shape PtrShape, ray const &RayIn)
{
  intersections Result{};

  if (PtrShape->isA<sphere>())
  {
    // Assert(0 == 1, __FUNCTION__, __LINE__);
    Result = LocalIntersectSphere(PtrShape, RayIn);
  }
  else if (PtrShape->isA<plane>())
  {
    Result = LocalIntersectPlane(PtrShape, RayIn);
  }

  return (Result);
}

//------------------------------------------------------------------------------
intersections Intersect(shared_ptr_shape PtrShape, ray const &Ray, ray *PtrLocalRayOutput)
{
  intersections Result{};

  // ---
  // NOTE: The object to which we are trying to calculate the intersect may
  //       kind of not be placed at origin. So use its transform to 'move' the
  //       ray by calculation of the inverse.
  //
  // Excerpt From The Ray Tracer Challenge Jamis Buck This material may be protected by copyright.
  // Another way to think about transformation matrices is to think of them as
  // converting points between two different coordinate systems. At the scene level,
  // everything is in world space coordinates, relative to the overall world.
  // But at the object level, everything is in object space coordinates, relative to the object itself.
  // Multiplying a point in object space by a transformation matrix converts that point
  // to world space—scaling it, translating, rotating it, or whatever.
  // Multiplying a point in world space by the inverse of the transformation matrix converts
  // that point back to object space.
  // ---
  // In other words: Whatever transformation you want to apply to the sphere, apply the inverse
  // of that transformation to the ray instead.
  ray const LocalRay = Transform(Ray, Inverse(PtrShape->Transform));

  // TODO: (Willy Clarke) Get rid of the extra variable
  PtrShape->SavedRay = LocalRay;

  // ---
  // NOTE: The local ray output is for testing only. It allows us to deliver the local
  //       computation of the transformed ray.
  // ---
  if (PtrLocalRayOutput)
  {
    *PtrLocalRayOutput = LocalRay;
  }

  if (PtrShape->funcPtrLocalIntersect)
  {
    Result = PtrShape->funcPtrLocalIntersect(PtrShape, LocalRay);
  }

  return (Result);
}

//------------------------------------------------------------------------------
intersection Intersection(float t, shared_ptr_shape pShape)
{
  intersection Result{};
  Result.t = t;
  Result.pShape = pShape;

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
bool Equal(material const &A, material const &B)
{
  bool Result{};
  Result = Equal(A.Ambient, B.Ambient) && Equal(A.Diffuse, B.Diffuse) && Equal(A.Shininess, B.Shininess) &&
           Equal(A.Specular, B.Specular) && Equal(A.Color, B.Color);
  return (Result);
}

//------------------------------------------------------------------------------
bool Equal(sphere const &A, sphere const &B)
{
  bool const EqMaterial = A.Material == B.Material;
  bool const EqTransform = A.Transform == B.Transform;
  return (Equal(A.Center, B.Center) && Equal(A.Radius, B.Radius) && EqMaterial && EqTransform);
}

//------------------------------------------------------------------------------
bool Equal(intersection const &A, intersection const &B)
{
  bool const Result = Equal(A.t, B.t);
  // std::cout << "Result:" << Result << ". Equal -> A:" << A.t << " B:" << B.t << std::endl;

  return (Result);
}

//------------------------------------------------------------------------------
bool Equal(light const &A, light const &B)
{
  bool const Result = Equal(A.Intensity, B.Intensity) && Equal(A.Position, B.Position);
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
    Assert(I.pShape, __FUNCTION__, __LINE__);
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
      Assert(Result.pShape, __FUNCTION__, __LINE__);
    }
  }

  return (Result);
}

/**
  Mul
  Multiply a matrix and a ray
  return : ray
*/
ray Mul(matrix const &M, ray const &R)
{
  ray Result{};
  Result.Origin = M * R.Origin;
  Result.Direction = M * R.Direction;
  return (Result);
}

//------------------------------------------------------------------------------
ray Transform(ray const &R, matrix const &M) { return M * R; }

/**
 * Convert the point to a vector.
 */
tup LocalNormalAt(shape const &Shape, tup const &LocalPoint)
{
  tup Result{};
  Result = LocalPoint - Point(0.f, 0.f, 0.f);
  return Result;
}

/**
 * The local normal for a plane is always 0, 1, 0.
 */
tup LocalNormalAt(plane const &Plane, tup const &LocalPoint)
{
  tup const Result{Vector(0.f, 1.f, 0.f)};
  return Result;
}

/**
 \fn NormalAt
 \brief Calculate normal vector at given point. The resulting vector will
        be normalized to a length of 1.f.
*/
tup NormalAt(shape const &Shape, tup const &PointInput)
{
  tup const LocalPoint = Inverse(Shape.Transform) * PointInput;
  tup const LocalNormal = LocalNormalAt(Shape, LocalPoint);
  tup WorldNormal = Transpose(Inverse(Shape.Transform)) * LocalNormal;
  WorldNormal.W = 0.f;

  tup const Result = Normalize(WorldNormal);
  return (Result);
}

//------------------------------------------------------------------------------
// \fn Reflect
// \brief Calculate the reflection vector based on the input and the surface normal.
// \sa https://math.stackexchange.com/questions/13261/how-to-get-a-reflection-vector
// \sa https://en.wikipedia.org/wiki/Vector_projection
//------------------------------------------------------------------------------
tup Reflect(tup const &In, tup const &Normal)
{
  tup Result{};
  Result = In - Normal * 2.f * Dot(In, Normal);
  return (Result);
}

//------------------------------------------------------------------------------
light PointLight(tup const &Position, tup const &Intensity)
{
  light Result{};
  Result.Position = Position;
  Result.Intensity = Intensity;
  return (Result);
}
//------------------------------------------------------------------------------
tup Lighting(material const &Material,  //!<
             light const &Light,        //!<
             tup const &Position,       //!<
             tup const &vEye,           //!<
             tup const &vNormal,        //!<
             bool const InShadow        //!<
)
{
  tup Result{};
  // Combine the surface color with the light''s color/intensity.
  tup const EffectiveColor = Material.Color * Light.Intensity;

  // Find the direction to the light source.
  tup const vLight = Normalize(Light.Position - Position);

  // Compute the ambient contribution
  tup const Ambient = EffectiveColor * Material.Ambient;

  // Light dot Normal represents the cosine of the angle between the
  // light vector and the normal vector. A negative number means that
  // the light is on the other side of the surface.
  float const LightDotNormal = Dot(vLight, vNormal);

  // Assert(LightDotNormal < 0.f, __FUNCTION__, __LINE__);

  tup Diffuse{};
  tup Specular{};

  if (LightDotNormal < 0.f)  // pointing away ...
  {
    Diffuse = Color(0.f, 0.f, 0.f);
    Specular = Color(0.f, 0.f, 0.f);
  }
  else
  {
    // Compute the diffuse contribution
    Diffuse = EffectiveColor * Material.Diffuse * LightDotNormal;

    // Reflect dot eye represents the cosine of the angle between the
    // reflection vector and the eye vector.
    // A negative value means that the light reflects away from the eye.
    tup vReflect = Reflect(-vLight, vNormal);
    float const ReflectDotEye = Dot(vReflect, vEye);
    if (ReflectDotEye <= 0)
    {
      Specular = Color(0.f, 0.f, 0.f);  // Black
    }
    else
    {
      // Compute the specular contribution
      float const Factor = std::pow(ReflectDotEye, Material.Shininess);
      Specular = Light.Intensity * Material.Specular * Factor;
    }
  }

  // Add the three contributions together to the the final shading
  // NOTE: When Inshadow we ignore the specular and diffuse component.
  if (InShadow)
  {
    Result = Ambient;
  }
  else
  {
    Result = Ambient + Diffuse + Specular;
  }

  return (Result);
};

/// ---
/// World functions.
/// ---
// \fn World - Create a default world.
world World()
{
  world W{};

  light Light = ww::PointLight(ww::Point(-10.f, 10.f, -10.f), ww::Color(1.f, 1.f, 1.f));
  // W.vLights.push_back(Light);
  {
    ww::shared_ptr_light PtrLight{};
    PtrLight.reset(new ww::light);
    *PtrLight = Light;
    W.vPtrLights.push_back(PtrLight);
  }

  sphere S1{};
  S1.Material.Color = ww::Color(0.8f, 1.0f, 0.6f);
  S1.Material.Diffuse = 0.7f;
  S1.Material.Specular = 0.2f;
  S1.funcPtrLocalIntersect = &ww::LocalIntersectSphere;

  {
    ww::shared_ptr_sphere PtrSphere{};
    PtrSphere.reset(new ww::sphere);
    *PtrSphere = S1;
    W.vPtrObjects.push_back(PtrSphere);
  }

  {
    ww::shared_ptr_sphere PtrSphere = ww::PtrDefaultSphere();
    PtrSphere->Transform = ww::Scaling(0.5f, 0.5f, 0.5f);
    PtrSphere->Radius = PtrSphere->Transform.R0.X;
    W.vPtrObjects.push_back(PtrSphere);
  }

  return (W);
}

//------------------------------------------------------------------------------
void WorldAddObject(world &W, shared_ptr_shape pShape) { W.vPtrObjects.push_back(pShape); }
//------------------------------------------------------------------------------
void WorldAddLight(world &W, shared_ptr_light pLight) { W.vPtrLights.push_back(pLight); }

//------------------------------------------------------------------------------
intersections IntersectWorld(world const &World, ray const &Ray)
{
  intersections XS{};

  for (auto const &PtrObject : World.vPtrObjects)
  {
    ww::ray LocalRayComputed{};

    intersections const I = Intersect(PtrObject, Ray, &LocalRayComputed);

#if 0
    if (PtrObject->isA<ww::sphere>())
    {
      ww::sphere const *pSphere = dynamic_cast<ww::sphere *>(PtrObject.get());
      std::cout << "\n\n " << __FUNCTION__ << "-> Sphere has radius " << pSphere->Radius << std::endl;
      std::cout << "\tLocalRayComputed: " << LocalRayComputed << std::endl;
    }
#endif

    for (auto const &Element : I.vI)
    {
      XS.vI.push_back(Element);
    }
  }

  // NOTE: Keep the intersections sorted in ascending order.
  //       Use lambda function to extract the t value for each intersection.
  std::sort(XS.vI.begin(), XS.vI.end(), [](intersection const &A, intersection const &B) { return A.t < B.t; });
  return (XS);
}

//------------------------------------------------------------------------------
shared_ptr_plane PtrDefaultPlane()
{
  ww::plane P{};
  std::shared_ptr<ww::plane> pPlane = ww::SharedPtrSh<ww::plane>(P);
  pPlane->funcPtrLocalIntersect = &ww::LocalIntersectPlane;
  return (pPlane);
}

//------------------------------------------------------------------------------
shared_ptr_sphere PtrDefaultSphere()
{
  shared_ptr_sphere pSphere{};
  pSphere.reset(new sphere);
  pSphere->funcPtrLocalIntersect = &ww::LocalIntersectSphere;
  return (pSphere);
}

//------------------------------------------------------------------------------
prepare_computation PrepareComputations(intersection const &I, ray const &R)
{
  Assert(I.pShape, __FUNCTION__, __LINE__);

  prepare_computation Comps{};

  // NOTE: Assign values we want to keep.
  Comps.t = I.t;
  Comps.pShape = I.pShape;

  // NOTE: Compute some useful values.
  Comps.Point = PositionAt(R, Comps.t);
  Comps.Eye = -R.Direction;
  Comps.Normal = NormalAt(*Comps.pShape, Comps.Point);

  // NOTE: Adjust Point for floating point inaccuracy.
  Comps.Point = Comps.Point + Comps.Normal * EPSILON;

  // NOTE: We use the dot product between the Normal and the Eye to figure out if the normal points
  //       away from the Eye. If negative they are (roughly) pointing in opposite directions.
  if (Dot(Comps.Normal, Comps.Eye) < 0.f)
  {
    Comps.Inside = true;  // NOTE: Default for the flag is false, no need to clear it once again.
    Comps.Normal = -Comps.Normal;
  }

  return (Comps);

  tup Test{};
}

//------------------------------------------------------------------------------
tup ShadeHit(world const &W, prepare_computation const &Comps)
{
  tup Color{};

  for (auto pWorldLight : W.vPtrLights)
  {
    light const &WorldLight = *pWorldLight;

    bool const Shadowed = IsShadowed(W, Comps.Point);

    tup C = Lighting(Comps.pShape->Material,  //!<
                     WorldLight,              //!<
                     Comps.Point,             //!<
                     Comps.Eye,               //!<
                     Comps.Normal,            //!<
                     Shadowed                 //!<
    );

    // NOTE: Add the colors from the various lights.
    Color = Color + C;
  }
  return (Color);
}

//------------------------------------------------------------------------------
tup ColorAt(world const &World, ray const &Ray)
{
  tup Result{};
  // 1. Call IntersectWorld() to find out the intersections of the given ray with the world.
  intersections const IS = IntersectWorld(World, Ray);

  // 2. Find the Hit from the resulting intersections.
  intersection const I = Hit(IS);

  // 3. Return the Color black if there is no such intersection.
  if (IS.Count() == 0) return Result;
  if (I.pShape == nullptr) return Result;
  Assert(IS.Count() != 0, __FUNCTION__, __LINE__);

  // 4. Otherwise pre-compute the necessary values with PrepareComputations
  Assert(I.pShape, __FUNCTION__, __LINE__);
  prepare_computation const PC = PrepareComputations(I, Ray);

  // 5. Call shade hit to find the color at the hit.
  Result = ShadeHit(World, PC);

  return (Result);
}

//------------------------------------------------------------------------------
matrix ViewTransform(tup const &From, tup const &To, tup const &Up)
{
  // 1. Compute the Forward vector by subtracting the From from To.
  tup const Forward = Normalize(To - From);

  // 2. Compute the Left vector by taking the cross product of Forward and the
  //    normalized Up vector.
  tup const UpNorm = Normalize(Up);
  tup const Left = Cross(Forward, UpNorm);

  // 3. Compute the TrueUp vector by taking the cross product of Left and Forward.
  //    This allows us to use an original Up vector to be approximately up. Framing
  //    a scene is a lot easier when we dont have to calculate a precise up vector.
  // tup const Left = Normalize(Cross(Forward, UpNorm));
  tup const TrueUp = Cross(Left, Forward);

  // 4. With Left, TrueUp and Forward we can now construct a matrix that represent
  //    an orientation transformation.
  matrix const Orientation = Matrix44(tup{Left.X, Left.Y, Left.Z, 0.f},              //!<
                                      tup{TrueUp.X, TrueUp.Y, TrueUp.Z, 0.f},        //!<
                                      tup{-Forward.X, -Forward.Y, -Forward.Z, 0.f},  //!<
                                      tup{0.f, 0.f, 0.f, 1.f}                        //!<
  );

  // 5. Finally, append a translation to the transformation to move the scene into place
  //    before orienting it. Multiply Orientation by Translation(-From.X, -From.Y, -From.Z)
  //    and we are done.
  matrix const Result = Orientation * Translation(-From.X, -From.Y, -From.Z);

  return (Result);
}

//------------------------------------------------------------------------------
camera Camera(int const HSize, int const VSize, float const FieldOfView)
{
  camera C{};

  C.HSize = HSize;
  C.VSize = VSize;
  C.FieldOfView = FieldOfView;

  // 1. You know the canvas is one unit away, and you know the angle of the field of view.
  //    By cutting the field of view in half, you create a right triangle,
  float HalfView = std::tan(C.FieldOfView / 2.f);

  // 2. The aspect ratio is the ratio of the Horisontal to the Vertical size.
  float const Aspect = float(C.HSize) / float(C.VSize);  // Who is going to create a canvas of zero height?

  // 3. Now if the Horisontal size is greater than or equal to the vertical size (aspect >= 1),
  //    then the HalfView is Half of the Width of the canvas and HalfView/Aspect is half the
  //    canvas' Height.
  //    If the vertial size is greater than the horizontal size (Aspect < 1), then the HalfView
  //    is instead the height of the canvas, and the half of the canvas' Width is HalfView * Aspect.
  if (Aspect >= 1.f)
  {
    C.HalfWidth = HalfView;
    C.HalfHeight = HalfView / Aspect;
  }
  else
  {
    C.HalfWidth = HalfView * Aspect;
    C.HalfHeight = HalfView;
  }

  // 4. Compute the size of a single pixel on the canvas by dividing the full width of the
  //    canvas (half width * 2) by the horizontal size in pixels of the canvas (HSize).
  //    This is the pixel size.
  C.PixelSize = (C.HalfWidth * 2.f) / C.HSize;

  return (C);
}

//------------------------------------------------------------------------------
ray RayForPixel(camera const &Camera, int const Px, int const Py)
{
  // The offset from edge of the canvas to the pixel's center
  float const XOffset = (Px + 0.5f) * Camera.PixelSize;
  float const YOffset = (Py + 0.5f) * Camera.PixelSize;

  // The untransformed coordinates of the pixel in world-space.
  // (remember that the camera looks toward -z, so +x is toward the *left*).
  float const WorldX = Camera.HalfWidth - XOffset;
  float const WorldY = Camera.HalfHeight - YOffset;

  // Using the camera matrix, transform the canvas point and the origin,
  // and compute the ray's direction vector.
  // Remember that the canvas is at z=-1.
  tup const Pixel = Inverse(Camera.Transform) * Point(WorldX, WorldY, -1.f);
  tup const Origin = Inverse(Camera.Transform) * Point(0.f, 0.f, 0.f);
  tup const Direction = Normalize(Pixel - Origin);
  ray const R = Ray(Origin, Direction);

  return (R);
}

//------------------------------------------------------------------------------
canvas Render(camera const &Camera, world const &World)
{
  canvas Image(Camera.HSize, Camera.VSize);

  for (int Y = 0;             ///<!
       Y < Camera.VSize - 1;  ///<!
       ++Y)
  {
    for (int X = 0;             ///<!
         X < Camera.HSize - 1;  ///<!
         ++X)
    {
      ray const R = RayForPixel(Camera, X, Y);
      tup const Color = ColorAt(World, R);
      WritePixel(Image, X, Y, Color);
    }
  }

  return (Image);
}

//------------------------------------------------------------------------------
bool IsShadowed(world const &World, tup const &Point)
{
  // NOTE: Count the number of times the point is in a shadow.
  size_t ShadowCount{};

  for (size_t Idx = 0;                 ///<!
       Idx < World.vPtrLights.size();  ///<!
       ++Idx)
  {
    // 1. Measure the distance from Point to the light source by subtracting
    //    Point from the light posistion, and taking the magnitude of the
    //    resulting vector. Call this distance.
    tup const V = World.vPtrLights[Idx]->Position - Point;
    Assert(IsVector(V), __FUNCTION__, __LINE__);
    float const Distance = Mag(V);
    tup const Direction = Normalize(V);

    // 2. Create a ray from Point toward the light source by normalizing the
    //    vector from step #1.
    ray const R = Ray(Point, Direction);

    // 3. Intersect the world with that ray.
    intersections const Intersections = IntersectWorld(World, R);

    // 4. Check to see if there was a hit, and if so, whether t is less than
    //    distance. If so, the hit lies between the Point and the light source,
    //    and the point is in shadow.
    intersection const H = Hit(Intersections);

    if (H.pShape && H.t < Distance)
    {
      ShadowCount++;
    }
  }

  // NOTE: The Point is in shadow only when it is in shadow from all light sources.
  return (ShadowCount == World.vPtrLights.size());
}

//------------------------------------------------------------------------------
// \brief TestShape
// \detail Function to demonstrate the abstract behavior of the shape class.
//         As shape itself is abstract, the TestShape instansiates and returns
//         a special subclass of shape that we will call test_shape, that implements
//         enough behavior to be concrete.
//------------------------------------------------------------------------------
shape TestShape()
{
  shape S{};
  // NOTE: Test that we can use the function pointer.
  S.funcPtrLocalIntersect = &ww::LocalIntersect;
  return (S);
}

//------------------------------------------------------------------------------
// \fn TestShapePtr
// \brief create an object and return a smart pointer to the object
// \return smart pointer to shape
shared_ptr_shape SharedPtrShape(shape const &Shape)
{
  ww::shared_ptr_shape PtrShape{};
  PtrShape.reset(new ww::shape);
  *PtrShape = Shape;

  return (PtrShape);
}

//------------------------------------------------------------------------------
/**
 * Plane
 */
tup Plane()
{
  tup Result{};
  return Result;
}

/**
 * StripePattern
 */
pattern StripePattern(tup const &C1, tup const &C2)
{
  return pattern{C1, C2};
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
bool operator==(ww::intersection const &A, ww::intersection const &B) { return (ww::Equal(A, B)); }
bool operator==(ww::material const &A, ww::material const &B) { return (ww::Equal(A, B)); }
bool operator==(ww::matrix const &A, ww::matrix const &B) { return (ww::Equal(A, B)); }
bool operator==(ww::light const &A, ww::light const &B) { return (ww::Equal(A, B)); }
bool operator==(ww::sphere const &A, ww::sphere const &B) { return (ww::Equal(A, B)); }
bool operator==(ww::tup const &A, ww::tup const &B) { return (ww::Equal(A, B)); }
// ---
// NOTE: Division operator does not check for divide by zero; Who cares?
// ---
ww::tup operator/(ww::tup const &Tup, float const S) { return (ww::Mul(1.f / S, Tup)); }
