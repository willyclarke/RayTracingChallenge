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
#include <list>    // for std::list
#include <memory>  // for shared pointer.
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "perlinnoise.hpp"

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
std::ostream &operator<<(std::ostream &stream, const ww::shape &S)
{
  stream << "\nShape:"
         << "\nShape Material:" << S.Material << "\nShape Transform:" << S.Transform << std::endl;
  return stream;
}

//------------------------------------------------------------------------------
std::ostream &operator<<(std::ostream &stream, const ww::sphere &S)
{
  stream << "\nSphere\nCenter:" << S.Center << "\nSphere Material:" << S.Material << "Sphere Radius:" << S.Radius
         << "\nSphere Transform:" << S.Transform << std::endl;
  return stream;
}

std::ostream &operator<<(std::ostream &stream, const ww::pattern &P)
{
  stream << "\n----\nPattern contents\nC1:" << P.A << "\nC2:" << P.B << "\nPrint: " << (P.Print ? "YES" : "NO")
         << "\nTransform:\n"
         << P.Transform << "\n-----" << std::endl;
  return stream;
}

std::ostream &operator<<(std::ostream &stream, const ww::intersections &XS)
{
  stream << "\nIntersections\n";
  int Idx{};
  for (auto const &I : XS.vI)
  {
    stream << "XS.vI[" << Idx << "]."
           << "I:" << I.pShape << ". t: " << I.t;
    if (I.pShape->isA<ww::sphere>())
    {
      stream << ". Sphere. Center: " << I.pShape->Center;
    }
    else if (I.pShape->isA<ww::cube>())
      stream << ". Cube";
    else if (I.pShape->isA<ww::plane>())
      stream << ". Plane";
    else if (I.pShape->isA<ww::shape>())
      stream << ". Shape";
    else
      stream << ". Undefinded shape type";

    stream << ". T:" << I.pShape->Transform;
    stream << std::endl;

    ++Idx;
  }
  return stream;
}

std::ostream &operator<<(std::ostream &stream, const ww::prepare_computation &Comps)
{
  stream << "\n---\nPrepareComputation. Address:" << Comps.pShape << std::endl;

  if (Comps.pShape->isA<ww::sphere>())
  {
    ww::sphere *ptrSphere = dynamic_cast<ww::sphere *>(Comps.pShape.get());
    if (ptrSphere)
    {
      stream << "\nSphere\nradius:" << ptrSphere->Radius;
      stream << "\ncenter:" << ptrSphere->Center << std::endl;
    }
  }
  stream << "\nn1:" << Comps.n1 << " (refration factor from where ray is leaving)."  //!<
         << "\nn2:" << Comps.n2 << " (refraction factor to where ray is entering)."  //!<
         << "\nHit at   t=" << Comps.t                                               //!<
         << "\nInside    :" << (Comps.Inside ? "Yes" : "No")                         //!<
         << "\n     Point:" << Comps.Point                                           //!<
         << "\nUnderPoint:" << Comps.UnderPoint                                      //!<
         << "\nOverPoint :" << Comps.OverPoint                                       //!<
         << "\nEye       :" << Comps.vEye                                            //!<
         << "\nNormal    :" << Comps.vNormal                                         //!<
         << "\nReflect   :" << Comps.vReflect << std::endl;
  if (Comps.pShape.get())
  {
    if (Comps.pShape->Material.Pattern.funcPtrPatternAt == &ww::FuncDefaultPatternAt)
      stream << "\nPattern: FuncDefaultPatternAt";
    else
      stream << "\nPattern: Unknown";
  }
  stream << "\n---" << std::endl;

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
/**
 * Return the maximum of tup and variable.
 * Works for both point and vector. Ignores W which is returned unaltered.
 */
tup Max(tup const &A, float const B) { return tup{std::max(A.X, B), std::max(A.Y, B), std::max(A.Z, B), A.W}; }

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
float Mag(tup const &Tup) { return (std::sqrtf(MagSquared(Tup))); }

//------------------------------------------------------------------------------
tup Abs(tup const &Tup)
{
  tup const Result{std::abs(Tup.X), std::abs(Tup.Y), std::abs(Tup.Z), std::abs(Tup.W)};
  return (Result);
}

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
/**
 * Convert from degrees to radians.
 */
float Radians(float Deg) { return (PI_F * Deg / 180.f); }
// ---
// NOTE: Canvas methods/functions.
// ---

//------------------------------------------------------------------------------
/**
 * Write a pixel color to the XY vector of the Canvas.
 * @Param Canvas: Destination for receiving the color.
 * @Param X: Offset in the X direction.
 * @Param Y: Offset in the Y direction.
 * @Param Color: Tuple containing the color.
 */
void WritePixel(canvas &Canvas, int X, int Y, tup const &Color)
{
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
  Result.ID.IsInvertible = true;
  Result.ID.Determinant = DetM;
  Result.ID.IsComputed = true;

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

/**
 * @Param: O - Origin
 * @Param: D - Direction
 * @Return: ray with normalized direction.
 */
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

  if (std::abs(Ray.Direction.Y) < EPSILON)
  {
    return Result;
  }

  float const t = -Ray.Origin.Y / Ray.Direction.Y;

  // ---
  // NOTE: Keep negative and positive t-hits. But avoid the ones where the ray has penetrated the plane.
  //       Ref RTC page 95: Identifying hits.
  // ---
  if (std::abs(t) > EPSILON)
  {
    Result.vI.push_back(Intersection(t, PtrShape));
  }

  return Result;
}

/**
 * Check if a ray has a local intersect with a sphere.
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

/**
 * Check if a ray has a local intersect with a cube.
 * Param: PtrShape: A ww::cube is needed to do some actual checks for intersections.
 * Return: ww::intersections.
 **/
ww::intersections LocalIntersectCube(shared_ptr_shape PtrShape, ray const &Ray)
{
  Assert(PtrShape->isA<cube>(), __FUNCTION__, __LINE__);

  ww::intersections XS{};

  if (!PtrShape->isA<ww::cube>()) return XS;

  auto const CheckAxis = [](float const Origin, float const Direction) -> std::pair<float, float>
  {
    std::pair<float, float> Result{};
    float &TMin = Result.first;
    float &TMax = Result.second;

    float const TMinNumerator = (-1.f - Origin);
    float const TMaxNumerator = (1.f - Origin);

    if (std::abs(Direction) >= EPSILON)
    {
      TMin = TMinNumerator / Direction;
      TMax = TMaxNumerator / Direction;
    }
    else
    {
      TMin = TMinNumerator * INIFINITY;
      TMax = TMaxNumerator * INIFINITY;
    }

    // ---
    // Swap when needed ...
    // ---
    if (TMin > TMax)
    {
      float const Tmp = TMin;
      TMin = TMax;
      TMax = Tmp;
    }

    return Result;
  };

  std::pair<float, float> const XtMinMax = CheckAxis(Ray.Origin.X, Ray.Direction.X);
  std::pair<float, float> const YtMinMax = CheckAxis(Ray.Origin.Y, Ray.Direction.Y);
  std::pair<float, float> const ZtMinMax = CheckAxis(Ray.Origin.Z, Ray.Direction.Z);

  // ---
  // Get the maximum value of minimums.
  // ---
  float TMin = std::max<float>(XtMinMax.first, YtMinMax.first);
  TMin = std::max<float>(TMin, ZtMinMax.first);

  // ---
  // Get the minimum value of the maximums.
  // ---
  float TMax = std::min<float>(XtMinMax.second, YtMinMax.second);
  TMax = std::min<float>(TMax, ZtMinMax.second);

  // ---
  // NOTE: Return no intersections when the ray misses the cube.
  // ---
  if (TMin > TMax) return {};

  intersection I0{TMin, PtrShape};
  intersection I1{TMax, PtrShape};

  XS.vI.push_back(I0);
  XS.vI.push_back(I1);

  return XS;
}

/**
 * Check if a ray has a local intersect with a cone.
 * Param: PtrShape: A cone is needed to do some actual checks for intersections.
 * Return: intersections.
 * NOTE: For reference :
 * stackoverflow:
 *https://gamedev.stackexchange.com/questions/112382/how-do-i-test-for-intersection-between-a-ray-and-a-cone
 *
 * In geometry quadric shapes are any surface that can be defined by an algebraic equation of second degree
 *
 * On the normal form this equation looks like this:
 *
 * Ax^2 + By^2 + Cz^2 +2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz + J = 0
 *
 * For the local intersect of a cone A, B and C above are 1 and the remaining capital letters are 0.
 **/
intersections LocalIntersectCone(shared_ptr_shape PtrShape, ray const &Ray)
{
  Assert(PtrShape->isA<cone>(), __FUNCTION__, __LINE__);

  auto IntersectCaps = [&](cone const *pCone, ray const &Ray, intersections &XS)
  {
    auto CheckCap = [](ray const &Ray, float t, float Radius) -> bool
    {
      float const X = Ray.Origin.X + t * Ray.Direction.X;
      float const Z = Ray.Origin.Z + t * Ray.Direction.Z;
      bool const Result = (X * X + Z * Z) <= std::abs(Radius);
      return Result;
    };

    // ---
    // NOTE: Caps only matter if the cone is closed, and might possibly be intersected
    //       by the ray.
    // ---
    if (!pCone->Closed) return;

    // ---
    // NOTE: The ray must not be horizontal in order to hit the cap.
    // ---
    if (std::abs(Ray.Direction.Y) < EPSILON) return;

    // ---
    // NOTE: Check for an intersection with the lower end cap by intersecting the ray
    //       with the plane a y=cyl.minimum.
    // ---
    {
      float const t = (pCone->Minimum - Ray.Origin.Y) / Ray.Direction.Y;
      if (CheckCap(Ray, t, pCone->Minimum)) XS.vI.push_back({t, PtrShape});
    }

    // ---
    // NOTE: Check for an intersection with the upper cap by intersecting
    //       the ray with the plane at y=cyl.maximum.
    // ---
    {
      float const t = (pCone->Maximum - Ray.Origin.Y) / Ray.Direction.Y;
      if (CheckCap(Ray, t, pCone->Maximum)) XS.vI.push_back({t, PtrShape});
    }
  };

  intersections XS{};
  cone const *pCone = dynamic_cast<cone *>(PtrShape.get());

  float const A = Ray.Direction.X * Ray.Direction.X          //!<
                  - Ray.Direction.Y * Ray.Direction.Y        //!<
                  + Ray.Direction.Z * Ray.Direction.Z;       //!<
  float const B = 2.f * (Ray.Origin.X * Ray.Direction.X      //!<
                         - Ray.Origin.Y * Ray.Direction.Y    //!<
                         + Ray.Origin.Z * Ray.Direction.Z);  //!<
  float const C = Ray.Origin.X * Ray.Origin.X                //!<
                  - Ray.Origin.Y * Ray.Origin.Y +            //!<
                  Ray.Origin.Z * Ray.Origin.Z;

  // ---
  // NOTE: The ray will miss when A and B are both 0. With A beeing 0 and B at somevalue there will be one hit.
  //       When the cone is closed there may be additional hits for the capped ends though.
  // ---
  if (std::abs(A) < EPSILON)
  {
    if (std::abs(B) < EPSILON)
    {
      return XS;
    }

    IntersectCaps(pCone, Ray, XS);
    float const t = -C / (2.f * B);

    XS.vI.push_back({t, PtrShape});
    return XS;
  }

  float Discriminant = B * B - 4.f * A * C;
  // ---
  // NOTE: No solutions when Discriminant i less than zero.
  // FIXME: (Willy Clarke) Check to see if InterectCaps need to be called here before returning.
  // ---
  if (Discriminant < -EPSILON) return {};
  if (Discriminant < 0.f) Discriminant = 0.f;  // FIXME: (Willy Clarke) What is the correct way of handling -0.f???

  float t0 = (-B - std::sqrtf(Discriminant)) / (2.f * A);
  float t1 = (-B + std::sqrtf(Discriminant)) / (2.f * A);

  if (t0 > t1)  // swap
  {
    float const tmp = t0;
    t0 = t1;
    t1 = tmp;
  }

  float const Y0 = Ray.Origin.Y + t0 * Ray.Direction.Y;
  if (pCone->Minimum < Y0 && Y0 < pCone->Maximum)
  {
    XS.vI.push_back({t0, PtrShape});
  }

  float const Y1 = Ray.Origin.Y + t1 * Ray.Direction.Y;
  if (pCone->Minimum < Y1 && Y1 < pCone->Maximum)
  {
    XS.vI.push_back({t1, PtrShape});
  }

  IntersectCaps(pCone, Ray, XS);

  return XS;
}

/**
 * Check if a ray has a local intersect with a cylinder.
 * Param: PtrShape: A cylinder is needed to do some actual checks for intersections.
 * Return: intersections.
 **/
intersections LocalIntersectCylinder(shared_ptr_shape PtrShape, ray const &Ray)
{
  Assert(PtrShape->isA<cylinder>(), __FUNCTION__, __LINE__);

  intersections XS{};
  cylinder const *pCylinder = dynamic_cast<cylinder *>(PtrShape.get());

  auto IntersectCaps = [&](cylinder const *pCylinder, ray const &Ray, intersections &XS)
  {
    auto CheckCap = [](ray const &Ray, float t) -> bool
    {
      float const X = Ray.Origin.X + t * Ray.Direction.X;
      float const Z = Ray.Origin.Z + t * Ray.Direction.Z;
      bool const Result = (X * X + Z * Z) <= 1.f;
      return Result;
    };

    // ---
    // NOTE: Caps only matter if the cylinder is closed, and might possibly be intersected
    //       by the ray.
    // ---
    if (!pCylinder->Closed) return;

    if (std::abs(Ray.Direction.Y) < EPSILON) return;

    // ---
    // NOTE: Check for an intersection with the lower end cap by intersecting the ray
    //       with the plane a y=cyl.minimum.
    // ---
    {
      float const t = (pCylinder->Minimum - Ray.Origin.Y) / Ray.Direction.Y;
      if (CheckCap(Ray, t)) XS.vI.push_back({t, PtrShape});
    }

    // ---
    // NOTE: Check for an intersection with the upper cap by intersecting
    //       the ray with the plane at y=cyl.maximum.
    // ---
    {
      float const t = (pCylinder->Maximum - Ray.Origin.Y) / Ray.Direction.Y;
      if (CheckCap(Ray, t)) XS.vI.push_back({t, PtrShape});
    }
  };

  float const A = Ray.Direction.X * Ray.Direction.X + Ray.Direction.Z * Ray.Direction.Z;

  // ---
  // NOTE: A ray is parallel to the y axis, return empty set.
  // ---
  if (std::abs(A) < EPSILON)
  {
    IntersectCaps(pCylinder, Ray, XS);
    return XS;
  }

  float const B = 2.f * Ray.Origin.X * Ray.Direction.X +  //!<
                  2.f * Ray.Origin.Z * Ray.Direction.Z;
  float const C = Ray.Origin.X * Ray.Origin.X + Ray.Origin.Z * Ray.Origin.Z - 1.f;
  float const Discriminant = B * B - 4 * A * C;

  // ---
  // NOTE: No solutions when Discriminant i less than zero.
  // FIXME: (Willy Clarke) Check to see if InterectCaps need to be called here before returning.
  // ---
  if (Discriminant < 0.f) return {};

  float t0 = (-B - std::sqrtf(Discriminant)) / (2 * A);
  float t1 = (-B + std::sqrtf(Discriminant)) / (2 * A);

  if (t0 > t1)  // swap
  {
    float const tmp = t0;
    t0 = t1;
    t1 = tmp;
  }

  float const Y0 = Ray.Origin.Y + t0 * Ray.Direction.Y;
  if (pCylinder->Minimum < Y0 && Y0 < pCylinder->Maximum)
  {
    XS.vI.push_back({t0, PtrShape});
  }

  float const Y1 = Ray.Origin.Y + t1 * Ray.Direction.Y;
  if (pCylinder->Minimum < Y1 && Y1 < pCylinder->Maximum)
  {
    XS.vI.push_back({t1, PtrShape});
  }

  // ---
  // NOTE: When the ray has less than two intersections it can still possibly intersect
  //       with the caps at either end.
  // ---
  if (XS.Count() < 2)
  {
    IntersectCaps(pCylinder, Ray, XS);
  }

  return XS;
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
  else if (PtrShape->isA<cube>())
  {
    Result = LocalIntersectCube(PtrShape, RayIn);
  }
  else if (PtrShape->isA<cone>())
  {
    Result = LocalIntersectCone(PtrShape, RayIn);
  }
  else if (PtrShape->isA<cylinder>())
  {
    Result = LocalIntersectCylinder(PtrShape, RayIn);
  }

  return (Result);
}

//------------------------------------------------------------------------------
intersections Intersect(shared_ptr_shape PtrShape, ray const &Ray, ray *PtrLocalRayOutput)
{
  intersections XS{};

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
    XS = PtrShape->funcPtrLocalIntersect(PtrShape, LocalRay);
  }

  return (XS);
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
intersections Intersections(intersection const &I)
{
  intersections XS{};
  XS.vI.push_back(I);
  Assert(I.pShape, __FUNCTION__, __LINE__);
  return (XS);
}

//------------------------------------------------------------------------------
intersections Intersections(intersection const &I1, intersection const &I2)
{
  intersections XS{};

  Assert(I1.pShape, __FUNCTION__, __LINE__);
  Assert(I2.pShape, __FUNCTION__, __LINE__);

  XS.vI.push_back(I1);
  XS.vI.push_back(I2);

  return (XS);
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
 * The local normal for a cube is the side that has absolute value of (1).
 * If none of the x,y,z components are 1 the point is not on the cube.
 */
tup LocalNormalAtCube(shape const &Cube, tup const &LocalPoint)
{
  tup Result = Vector(0.f, 0.f, LocalPoint.Z);

  float MaxC = std::max<float>(std::abs(LocalPoint.X), std::abs(LocalPoint.Y));
  MaxC = std::max<float>(MaxC, std::abs(LocalPoint.Z));

  if (MaxC == std::abs(LocalPoint.X))
  {
    Result = Vector(LocalPoint.X, 0.f, 0.f);
  }
  else if (MaxC == std::abs(LocalPoint.Y))
  {
    Result = Vector(0.f, LocalPoint.Y, 0.f);
  }

  return Result;
}

/**
 * The local normal for a cone is the ...
 */
tup LocalNormalAtCone(shape const &Cone, tup const &LocalPoint)
{
  float Y = std::sqrt(LocalPoint.X * LocalPoint.X + LocalPoint.Z * LocalPoint.Z);
  if (LocalPoint.Y > 0.f) Y = -Y;

  tup Result = Vector(LocalPoint.X, Y, LocalPoint.Z);

  cone const *pCone = dynamic_cast<cone const *>(&Cone);

  // ---
  // NOTE: Compute the square of the distance from the y axis.
  // ---
  float const Distance = LocalPoint.X * LocalPoint.X + LocalPoint.Z * LocalPoint.Z;

  if (Distance < 1.f && LocalPoint.Y >= pCone->Maximum - EPSILON)
  {
    Result = Vector(0.f, 1.f, 0.f);
  }
  else if (Distance < 1.f && LocalPoint.Y <= pCone->Minimum + EPSILON)
  {
    Result = Vector(0.f, -1.f, 0.f);
  }
  return Result;
}

/**
 * The local normal for a cylinder is the point on the cylinder with the
 * Y-component removed.
 */
tup LocalNormalAtCylinder(shape const &Cylinder, tup const &LocalPoint)
{
  tup Result = Vector(LocalPoint.X, 0.f, LocalPoint.Z);

  cylinder const *pCylinder = dynamic_cast<cylinder const *>(&Cylinder);

  // ---
  // NOTE: Compute the square of the distance from the y axis.
  // ---
  float const Distance = LocalPoint.X * LocalPoint.X + LocalPoint.Z * LocalPoint.Z;

  if (Distance < 1.f && LocalPoint.Y >= pCylinder->Maximum - EPSILON)
  {
    Result = Vector(0.f, 1.f, 0.f);
  }
  else if (Distance < 1.f && LocalPoint.Y <= pCylinder->Minimum + EPSILON)
  {
    Result = Vector(0.f, -1.f, 0.f);
  }

  return Result;
}

/**
 * The local normal for a plane is always 0, 1, 0.
 */
tup LocalNormalAtPlane(shape const &Plane, tup const &LocalPoint)
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
  // tup const LocalNormal = LocalNormalAt(Shape, LocalPoint);
  tup const LocalNormal = Shape.funcPtrLocalNormalAt(Shape, LocalPoint);
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
             shape const &Object,       //!<
             light const &Light,        //!<
             tup const &Position,       //!<
             tup const &vEye,           //!<
             tup const &vNormal,        //!<
             bool const InShadow        //!<
)
{
  tup Result{};

  // ---
  // NOTE: Initialize to the material color, then check if a pattern has been applied.
  // ---
  tup Color = Material.Color;

  if (Material.Pattern != pattern{})
  {
    Color = PatternAtShape(Material.Pattern, Object, Position);
  }

  // Combine the surface color with the light''s color/intensity.
  tup const EffectiveColor = Color * Light.Intensity;

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
    Diffuse = ww::Color(0.f, 0.f, 0.f);
    Specular = ww::Color(0.f, 0.f, 0.f);
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
      Specular = ww::Color(0.f, 0.f, 0.f);  // Black
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

  {
    ww::shared_ptr_plane Floor = ww::PtrDefaultPlane();
    Floor->Transform = ww::Translation(0.f, -1.f, 0.f);
    Floor->Material.Transparency = 0.5f;
    Floor->Material.RefractiveIndex = 1.5f;  //!< Glass
    W.vPtrObjects.push_back(Floor);
  }

  {
    ww::shared_ptr_sphere Ball = ww::PtrDefaultSphere();
    Ball->Material.Color = ww::Color(1.f, 0.f, 0.f);
    Ball->Material.Ambient = 0.5f;
    Ball->Transform = ww::Translation(0.f, -3.5f, -0.5f);
    W.vPtrObjects.push_back(Ball);
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
    ww::ray &LocalRayComputed = World.LocalRayComputed;

    intersections const XSPerObject = Intersect(PtrObject, Ray, &LocalRayComputed);

#if 0
    if (World.Print)
    {
      if (PtrObject->isA<ww::sphere>())
      {
        ww::sphere const *pSphere = dynamic_cast<ww::sphere *>(PtrObject.get());
        std::cout << "\n\n " << __FUNCTION__ << "-> Sphere has radius " << pSphere->Radius << std::endl;
        std::cout << "\nLocalRayComputed:\n" << LocalRayComputed << std::endl;
      }
      if (PtrObject->isA<ww::plane>())
      {
        ww::plane const *pPlane = dynamic_cast<ww::plane *>(PtrObject.get());
        std::cout << "\n\n " << __FUNCTION__ << "-> SavedRay:\n" << pPlane->SavedRay << std::endl;
        std::cout << "\nLocalRayComputed:\n" << LocalRayComputed << std::endl;
      }
      std::cout << "Material.Pattern.A: " << PtrObject->Material.Pattern.A << std::endl;
      std::cout << "Material.Pattern.B: " << PtrObject->Material.Pattern.B << std::endl;
    }
#endif

    for (auto const &I : XSPerObject.vI)
    {
      if (I.pShape) XS.vI.push_back(I);
    }
  }

  // NOTE: Keep the intersections sorted in ascending order.
  //       Use lambda function to extract the t value for each intersection.
  std::sort(XS.vI.begin(), XS.vI.end(), [](intersection const &A, intersection const &B) { return A.t < B.t; });
  return (XS);
}

//------------------------------------------------------------------------------
shared_ptr_cone PtrDefaultCone()
{
  shared_ptr_cone pCone = SharedPtrSh<cone>(cone{});
  pCone->funcPtrLocalIntersect = &ww::LocalIntersectCone;
  pCone->funcPtrLocalNormalAt = &ww::LocalNormalAtCone;
  return pCone;
}

//------------------------------------------------------------------------------
shared_ptr_cone PtrCappedCone(float Minimum, float Maximum)
{
  shared_ptr_cone pCone = PtrDefaultCone();
  pCone->Closed = true;
  pCone->Minimum = Minimum;
  pCone->Maximum = Maximum;
  return pCone;
}

//------------------------------------------------------------------------------
shared_ptr_cylinder PtrDefaultCylinder()
{
  shared_ptr_cylinder pCylinder = SharedPtrSh<cylinder>(cylinder{});
  pCylinder->funcPtrLocalIntersect = &ww::LocalIntersectCylinder;
  pCylinder->funcPtrLocalNormalAt = &ww::LocalNormalAtCylinder;
  return pCylinder;
}

//------------------------------------------------------------------------------
shared_ptr_cylinder PtrCappedCylinder(float Minimum, float Maximum)
{
  shared_ptr_cylinder pCylinder = PtrDefaultCylinder();
  pCylinder->Closed = true;
  pCylinder->Minimum = Minimum;
  pCylinder->Maximum = Maximum;
  return pCylinder;
}

//------------------------------------------------------------------------------
shared_ptr_plane PtrDefaultPlane()
{
  ww::plane P{};
  std::shared_ptr<ww::plane> pPlane = ww::SharedPtrSh<ww::plane>(P);
  pPlane->funcPtrLocalIntersect = &ww::LocalIntersectPlane;
  pPlane->funcPtrLocalNormalAt = &ww::LocalNormalAtPlane;
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
shared_ptr_sphere PtrGlassSphere()
{
  shared_ptr_sphere pSphere{};
  pSphere.reset(new sphere);
  pSphere->Material.Transparency = 1.f;
  pSphere->Material.RefractiveIndex = 1.5f;
  pSphere->funcPtrLocalIntersect = &ww::LocalIntersectSphere;
  return (pSphere);
}

//------------------------------------------------------------------------------
shared_ptr_cube PtrDefaultCube()
{
  shared_ptr_cube pCube{};
  pCube.reset(new cube);
  pCube->funcPtrLocalIntersect = &ww::LocalIntersectCube;
  pCube->funcPtrLocalNormalAt = &ww::LocalNormalAtCube;
  return (pCube);
}

//------------------------------------------------------------------------------
prepare_computation PrepareComputations(intersection const &Hit, ray const &R, intersections const *ptrXS)
{
  Assert(Hit.pShape, __FUNCTION__, __LINE__);

  prepare_computation Comps{};

  // NOTE: Assign values we want to keep.
  Comps.t = Hit.t;
  Comps.pShape = Hit.pShape;

  // NOTE: Compute some useful values.
  Comps.Point = PositionAt(R, Comps.t);
  Comps.vEye = -R.Direction;
  Comps.vNormal = NormalAt(*Comps.pShape, Comps.Point);
  Comps.vReflect = Reflect(R.Direction, Comps.vNormal);
  Comps.OverPoint = Comps.Point + Comps.vNormal * EPSILON;
  Comps.UnderPoint = Comps.Point - Comps.vNormal * EPSILON;

  // NOTE: We use the dot product between the Normal and the Eye to figure out if the normal points
  //       away from the Eye. If negative they are (roughly) pointing in opposite directions.
  if (Dot(Comps.vNormal, Comps.vEye) < 0.f)
  {
    Comps.Inside = true;  // NOTE: Default for the flag is false, no need to clear it once again.
    Comps.vNormal = -Comps.vNormal;
  }

  // ---
  // NOTE: Determine the refractive index.
  //       The variables n1 and n2 is the refractive index at either side of a ray-object
  //       intersection.
  // ---
  if (ptrXS)
  {
    std::list<shared_ptr_shape> Container{};

    for (size_t Idx = 0;        //!<
         Idx < ptrXS->Count();  //!<
         ++Idx)
    {
      intersection const &I = ptrXS->vI[Idx];
      // ---
      // NOTE: 1. If the intersection is the hit, set the refractive index to the refractive index of the last object in
      // the
      //          containers list. If that list is empty, then there is no containing object and n1 should be set to 1.
      // ---
      if (I == Hit)
      {
        if (Container.empty())
        {
          Comps.n1 = 1.f;
        }
        else
        {
          Comps.n1 = Container.back()->Material.RefractiveIndex;
        }
      }

      // ---
      // NOTE: 2. If the intersection's object is already in the Containers list, then this intersection must be exiting
      //          the object. Remove the object from the Containers list in this case. Otherwise, the intersection is
      //          entering the object, and the object should be added to the end of the list.
      // ---
      //
      // ---
      // NOTE: Look for intersection I in the Container with all the shapes.
      // ---
      Assert(I.pShape, __FUNCTION__, __LINE__);
      auto it = std::find(Container.begin(), Container.end(), I.pShape);
      if (it != Container.end())  // found it so now delete the intersection.
      {
        Container.remove(*it);
      }
      else  // it was not there so add the object with the intersection.
      {
        Container.push_back(I.pShape);
      }

      // ---
      // NOTE: 3. If the intersection is the Hit, set n2 to the refractive index of the last object in the containers
      //          list. If that list is empty, then again, there is no containing object and n2 should be set to 1.
      // ---
      if (I == Hit)
      {
        if (Container.empty())
        {
          Comps.n2 = 1.f;
        }
        else
        {
          Comps.n2 = Container.back()->Material.RefractiveIndex;
        }

        // ---
        // NOTE: Since the Intersection is the Hit; break here.
        // ---
        break;
      }
    }  // end for
  }    // end ptrXS

  return (Comps);
}

//------------------------------------------------------------------------------
tup ShadeHit(world const &World, prepare_computation const &Comps, int const Remaining)
{
  tup Color{};

  for (auto pWorldLight : World.vPtrLights)
  {
    light const &WorldLight = *pWorldLight;

    bool const Shadowed = IsShadowed(World, Comps.OverPoint);

    tup const Surface = Lighting(Comps.pShape->Material,  //!<
                                 *Comps.pShape,           //!<
                                 WorldLight,              //!<
                                 Comps.OverPoint,         //!<
                                 Comps.vEye,              //!<
                                 Comps.vNormal,           //!<
                                 Shadowed                 //!<
    );

    tup const Reflected = ReflectedColor(World, Comps, Remaining);
    tup const Refracted = RefractedColor(World, Comps, Remaining);

    // ---
    // NOTE: Add the colors from the various lights.
    // ---
    material const &Material = Comps.pShape->Material;
    if (Material.Reflective > 0.f && Material.Transparency > 0.f)
    {
      float const Reflectance = Schlick(Comps);
      Color = Color + Surface + Reflected * Reflectance + Refracted * (1.f - Reflectance);
    }
    else
    {
      Color = Color + Surface + Reflected + Refracted;
    }
  }
  return (Color);
}

//------------------------------------------------------------------------------
tup ColorAt(world const &World, ray const &Ray, int const Remaining)
{
  tup Result{};
  // 1. Call IntersectWorld() to find out the intersections of the given ray with the world.
  intersections const XS = IntersectWorld(World, Ray);

  // 2. Find the Hit from the resulting intersections.
  intersection const I = Hit(XS);

  // 3. Return the Color black if there is no such intersection.
  if (XS.Count() == 0)
  {
    return Result;
  }

  if (I.pShape == nullptr)
  {
    return Result;
  }
  Assert(XS.Count() != 0, __FUNCTION__, __LINE__);

  // 4. Otherwise pre-compute the necessary values with PrepareComputations
  Assert(I.pShape, __FUNCTION__, __LINE__);
  prepare_computation const Comps = PrepareComputations(I, Ray, &XS);

  // 5. Call shade hit to find the color at the hit.
  Result = ShadeHit(World, Comps, Remaining);

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
void RenderSingleThread(camera const &Camera, world const &World, canvas &Image)
{
  for (int Y = 0;         ///<!
       Y < Camera.VSize;  ///<!
       ++Y)
  {
    for (int X = 0;         ///<!
         X < Camera.HSize;  ///<!
         ++X)
    {
      ray const R = RayForPixel(Camera, X, Y);
      tup const Color = ColorAt(World, R);
      WritePixel(Image, X, Y, Color);
    }
  }
}

/**
 * Render a block of pixels defined in the struct render_block.
 * This function locks on a mutex defined in the render block.
 */
void RenderBlock(render_block const &RB)
{
  for (int Y = RB.VStart;           ///<!
       Y < RB.VStart + RB.VHeigth;  ///<!
       ++Y)
  {
    for (int X = RB.HStart;           ///<!
         X < RB.HStart + RB.HLength;  ///<!
         ++X)
    {
      ray const R = RayForPixel(*RB.ptrCamera, X, Y);
      tup const Color = ColorAt(*RB.ptrWorld, R);
      WritePixel(*RB.ptrImage, X, Y, Color);
    }
  }
}

/**
 * Use a number of threads to render the image on to the canvas.
 * NOTE: The function is not thread safe but since the canvas is
 *       split into blocks that are not overlapping there should (?)
 *       not be any undefinded behavior.
 */
void RenderMultiThread(camera const &Camera, world const &World, canvas &Image)
{
  int const NumBlocksH = Camera.NumBlocksH;
  int const NumBlocksV = Camera.NumBlocksV;
  int const HSizeBlock = Camera.HSize / NumBlocksH;
  int const VSizeBlock = Camera.VSize / NumBlocksV;

  std::vector<render_block> vRenderBlocks{};

  for (int HIdx = 0;       //!<
       HIdx < NumBlocksH;  //!<
       ++HIdx)
  {
    for (int VIdx = 0;       //!<
         VIdx < NumBlocksV;  //!<
         ++VIdx)
    {
      render_block RB{};
      RB.HLength = HSizeBlock;
      RB.VHeigth = VSizeBlock;
      RB.HStart = HIdx * HSizeBlock;
      RB.VStart = VIdx * VSizeBlock;
      RB.ptrImage = &Image;
      RB.ptrCamera = &Camera;
      RB.ptrWorld = &World;
      vRenderBlocks.push_back(RB);
    }
  }

  struct WorkerThread
  {
    std::thread T;
    render_block _RB{};
    WorkerThread() { std::cout << __PRETTY_FUNCTION__ << " -> Called XXXXXXXXXXXXXXXXXXXXXXXXXX " << std::endl; }
    WorkerThread(render_block const &RB) { _RB = RB; }

    void Start() { T = std::thread(RenderBlock, _RB); }

    WorkerThread(WorkerThread const &Orig) { _RB = Orig._RB; };
    ~WorkerThread()
    {
      if (T.joinable())
      {
        T.join();
      }
    };
  };

  std::vector<WorkerThread> vWorkerThreads{};

  for (auto const &RB : vRenderBlocks)
  {
    vWorkerThreads.push_back(RB);
  }

  for (auto &T : vWorkerThreads)
  {
    T.Start();
  }
}

//------------------------------------------------------------------------------
canvas Render(camera const &Camera, world const &World)
{
  canvas Image(Camera.HSize, Camera.VSize);

  if (Camera.RenderSingleThread)
    RenderSingleThread(Camera, World, Image);
  else
    RenderMultiThread(Camera, World, Image);
  return Image;
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

/**
 * Use the prepare_computation Reflect vector to calculate the color to return.
 * When the Reflect is 0 the color will be black.
 */
tup ReflectedColor(world const &World, prepare_computation const &Comps, int const Remaining)
{
  if ((Remaining <= 0) || (Comps.pShape->Material.Reflective < EPSILON))
  {
    return ww::Color(0.f, 0.f, 0.f);
  }

  ray const ReflectRay = ray{Comps.OverPoint, Comps.vReflect};
  tup Color = ColorAt(World, ReflectRay, Remaining - 1) * Comps.pShape->Material.Reflective;
  return Color;
}

namespace
{

/**
 */
tup RefractedColorVectorBased(world const &World, prepare_computation const &Comps, int const Remaining)
{
  // ---
  // NOTE: Find the ratio of the first index to the second.
  // ---
  float const nRatio = Comps.n1 / Comps.n2;  //!< This is the inverted defintion from Snells law.

  float const CosI = Dot(Comps.vEye, Comps.vNormal);  //!< Cosine(thetaI) : the dot product between the vectors.
  Assert(CosI > 0.f, __FUNCTION__, __LINE__);  //!< CosI must be positive when the normal points to the light source.

  float const Sin2T = nRatio * nRatio * (1.f - CosI * CosI);  //!< Find sin(ThetaT)^2 by trigonometric identity.

  bool const TotalReflection = Sin2T > 1.f;  //!< There is total reflection going on ?
  if (TotalReflection || (0 >= Remaining) || (Comps.pShape->Material.Transparency == 0.f))
  {
    return {};  //!< Return Black when ...
  }

  // ---
  // NOTE: alternative calc - from Wikipedia.
  // ---
  float const c = -Dot(Comps.vNormal, Comps.vEye);
  float const r = Comps.n1 / Comps.n2;
  tup const &l = Comps.vEye;
  tup const &n = Comps.vNormal;
  tup vRefract = r * l + (r * c - std::sqrtf(1 - r * r * (1 - c * c))) * n;
  vRefract.W = 0.f;
  // ---
  //
  if (Comps.PrintDebug)
  {
    std::cout << __FUNCTION__ << "." << __LINE__ << std::endl;
    std::cout << Comps << std::endl;
    std::cout << "\nvRefract (Wikipedia): " << vRefract << std::endl;
  }

  // ---
  // NOTE: Create the refracted ray.
  // ---
  ray const RefractRay = Ray(Comps.Inside ? Comps.OverPoint : Comps.UnderPoint, vRefract);

  if (Comps.PrintDebug)
  {
    std::cout << __FUNCTION__ << "." << __LINE__ << ": RefractedRay:\n" << RefractRay << std::endl;
  }

  // ---
  // NOTE: Find the color of the refracted ray making sure to multiply
  //       with the transparency value to account for any opacity.
  // ---
  tup const Color = ColorAt(World, RefractRay, Remaining - 1) * Comps.pShape->Material.Transparency;
  std::cout << __FUNCTION__ << "." << __LINE__ << ". Color (Wikipedia):" << Color << std::endl;
  std::cout << "xxxxxxxx" << std::endl;
  return Color;
}

/**
 */
tup RefractedColorCosineBased(world const &World, prepare_computation const &Comps, int const Remaining)
{
  // ---
  // NOTE: Find the ratio of the first index to the second.
  // ---
  float const nRatio = Comps.n1 / Comps.n2;  //!< This is the inverted defintion from Snells law.

  // ---
  // NOTE: Snell law ->
  //                      sin(ThetaI)       n2
  //                      -------------  = -----
  //                      sin(ThetaT)       n1

  float const CosI = Dot(Comps.vEye, Comps.vNormal);  //!< Cosine(thetaI) : the dot product between the vectors.
  Assert(CosI > 0.f, __FUNCTION__, __LINE__);  //!< CosI must be positive when the normal points to the light source.

  float const Sin2T = nRatio * nRatio * (1.f - CosI * CosI);  //!< Find sin(ThetaT)^2 by trigonometric identity.

  bool const TotalReflection = Sin2T > 1.f;  //!< There is total reflection going on ?
  if (TotalReflection || (0 >= Remaining) || (Comps.pShape->Material.Transparency == 0.f))
  {
    return {};  //!< Return Black when ...
  }

  float const CosT = std::sqrtf(1.f - Sin2T);  //!< Find cos(ThetaT) via trigonometric identity.

  // ---
  // NOTE: Compute the direction of the refracted ray.
  // ---
  tup const vDirection = Comps.vNormal * (nRatio * CosI - CosT) - Comps.vEye * nRatio;

  // ---
  // NOTE: Create the refracted ray.
  // ---
  // ray const RefractRay = Ray(Comps.Point, vDirection);
  // ray const RefractRay = Ray(Comps.UnderPoint, vDirection);
  ray const RefractRay = Ray(Comps.Inside ? Comps.OverPoint : Comps.UnderPoint, vDirection);

  // ---
  // NOTE: Find the color of the refracted ray making sure to multiply
  //       with the transparency value to account for any opacity.
  // ---
  tup const Color = ColorAt(World, RefractRay, Remaining - 1) * Comps.pShape->Material.Transparency;
  return Color;
}
};  // end of anonymous namespace

/**
 */
tup RefractedColor(world const &World, prepare_computation const &Comps, int const Remaining)
{
  tup Color = RefractedColorCosineBased(World, Comps, Remaining);
  return Color;
  Color = RefractedColorVectorBased(World, Comps, Remaining);
}

/**
 */
float Schlick(prepare_computation const &Comps)
{
  float Result{};

  // ---
  // NOTE: Find the cosine between the Eye and the Normal vector.
  // ---
  float Cos = Dot(Comps.vEye, Comps.vNormal);

  // ---
  // NOTE: Total internal reflection can only occur if n1 > n2.
  // ---

  if (Comps.n1 > Comps.n2)
  {
    float const n = Comps.n1 / Comps.n2;
    float const Sin2T = n * n * (1.f - Cos * Cos);
    if (Sin2T > 1.f) return 1.f;

    // ---
    // NOTE: Compute the cosine of ThetaT by using trig identity.
    // ---
    float const CosT = std::sqrtf(1.f - Sin2T);

    // ---
    // NOTE: When n1 > n2 use CosT as Cos instead.
    // ---
    Cos = CosT;
  }

  float const r = ((Comps.n1 - Comps.n2) / (Comps.n1 + Comps.n2));
  float const r0 = r * r;
  float const OneMinCos = 1.f - Cos;
  float const OneMinCosPow5 = OneMinCos * OneMinCos * OneMinCos * OneMinCos * OneMinCos;
  Result = r0 + (1.f - r0) * OneMinCosPow5;

  return Result;
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
 * Test Pattern
 */
pattern TestPattern()
{
  pattern P{};
  P.A = ww::Color(1.f, 1.f, 1.f);
  P.B = ww::Color(1.f, 1.f, 1.f);
  return P;
}

/**
 * Set the transform into the pattern and return a reference.
 */
pattern &SetPatternTransform(pattern &P, matrix const &T)
{
  P.Transform = T;
  return P;
}

/**
 * Solid Pattern - Return the same color for every point.
 */
pattern SolidPattern(tup const &Color, char const *ptr)
{
#if 0
  std::cout << __PRETTY_FUNCTION__ << "." << __LINE__ << ". Called from " << std::string(ptr) << std::endl;
#endif

  pattern P{Color, Color};
  P.funcPtrPatternAt = SolidPatternAt;

  // ---
  // NOTE: Now create a new pattern based on P1 which will call whatever
  //       funcPtrPatternAt that has been set up for it.
  // ---
  P.ptrNext = std::make_shared<pattern>();
  if (P.ptrNext != nullptr)
  {
    *P.ptrNext = P;

    // ---
    // NOTE: The end of the linked list of patterns is set up to whatever
    // ---   funcPtrPatternAt has been set up for P.
    P.ptrNext->ptrNext = std::make_shared<pattern>();

    if (P.ptrNext->ptrNext != nullptr) *P.ptrNext->ptrNext = P;
  }
  return P;
}

/**
 */
pattern StripePattern(tup const &C1, tup const &C2)
{
  pattern P = SolidPattern(C1, __PRETTY_FUNCTION__);
  P.A = C1;
  P.B = C2;
  P.funcPtrPatternAt = StripeAt;
  return P;
}

/**
 */
pattern CheckersPattern(tup const &C1, tup const &C2)
{
  pattern P{C1, C2};
  P.funcPtrPatternAt = CheckersPatternAt;
  return P;
}

/**
 */
pattern CheckersGradientPattern(tup const &C1, tup const &C2)
{
  pattern P = SolidPattern(C1, __PRETTY_FUNCTION__);
  P.A = C1;
  P.B = C2;
  P.funcPtrPatternAt = CheckersGradientPatternAt;
  return P;
}

/**
 */
pattern RingPattern(tup const &C1, tup const &C2)
{
  pattern P = SolidPattern(C1, __PRETTY_FUNCTION__);
  P.A = C1;
  P.B = C2;
  P.funcPtrPatternAt = RingPatternAt;
  return P;
}

/**
 */
pattern RadialGradientPattern(tup const &C1, tup const &C2)
{
  pattern P = SolidPattern(C1, __PRETTY_FUNCTION__);
  P.A = C1;
  P.B = C2;
  P.funcPtrPatternAt = RadialGradientPatternAt;
  return P;
}

/**
 */
pattern GradientPattern(tup const &C1, tup const &C2)
{
  pattern P = SolidPattern(C1, __PRETTY_FUNCTION__);
  P.A = C1;
  P.B = C2;
  P.funcPtrPatternAt = GradientPatternAt;
  return P;
}

/**
 * Set up function pointers for blending two patterns.
 * @Param: P1 is used for the first link pattern.
 * @Param: P2 is used for the second link pattern.
 * @Return: P which is set up to call BlendedPatternAt.
 * @Note: Input patterns may not use BlendedPattern as funcPtrPatternAt
 *        since that would lead to an infinite recursion.
 */
pattern BlendedPattern(pattern const &P1, pattern const &P2)
{
  pattern P = SolidPattern(P1.A, __PRETTY_FUNCTION__);
  P.A = P1.A;
  P.B = P1.B;

  // ---
  // NOTE: Set up the resulting pattern to call the blended pattern function.
  // ---
  P.funcPtrPatternAt = BlendedPatternAt;
  *P.ptrNext = P1;
  *P.ptrNext->ptrNext = P2;

  // ---
  // NOTE: Check that the function pointers are not pointing to the BlendedPattern for its funcPtrPatternAt.
  // ---
  auto const ptrResultFunc = P.funcPtrPatternAt;
  auto const ptrP1Func = P1.funcPtrPatternAt;
  auto const ptrP2Func = P2.funcPtrPatternAt;
  if (ptrResultFunc == ptrP1Func || ptrResultFunc == ptrP2Func)
  {
    std::cerr << "ERROR: " << __PRETTY_FUNCTION__ << ": Nested <pattern> functions would lead to infinite recursion"
              << std::endl;
    Assert(false, __FUNCTION__, __LINE__);
    return {};
  }

  return P;
}

/**
 * Set up a linked list with three patterns.
 * @Param: PMain - The main pattern. e.g Checkers.
 * @Param: P1 - Pattern used when selector 1 triggers.
 * @Param: P2 - Pattern used when selector 2 triggers.
 * @Note: The nested patterns can not be of the same type as the
 *        main pattern since that would lead to an infinite
 *        recursion.
 * @Return: A pattern that has a linked list to two other patterns.
 */
pattern NestedPattern(pattern const &PMain, pattern const &P1, pattern const &P2)
{
  pattern Result = PMain;

  Result.ptrNext = std::make_shared<pattern>();
  *Result.ptrNext = P1;
  Result.ptrNext->ptrNext = std::make_shared<pattern>();
  *Result.ptrNext->ptrNext = P2;

  if (Result.ptrNext == nullptr || ((Result.ptrNext != nullptr) && (Result.ptrNext->ptrNext == nullptr)))
  {
    std::cerr << __FUNCTION__ << ". ERROR: Could not create shared pointers to <pattern>." << std::endl;
    Assert(false, __FUNCTION__, __LINE__);
    return {};
  }

  // ---
  // NOTE: Check that the function pointers are not pointing to the Main pattern funcPtrPatternAt.
  // ---
  auto const ptrMainFunc = PMain.funcPtrPatternAt;
  auto const ptrP1Func = P1.funcPtrPatternAt;
  auto const ptrP2Func = P2.funcPtrPatternAt;
  if (ptrMainFunc == ptrP1Func || ptrMainFunc == ptrP2Func)
  {
    std::cerr << "ERROR: " << __PRETTY_FUNCTION__ << ": Nested <pattern> functions would lead to infinite recursion"
              << std::endl;
    Assert(false, __FUNCTION__, __LINE__);
    return {};
  }

  return Result;
}

/**
 * Set up a perturbed pattern so that the input pattern is disturbed by noise.
 */
pattern PerturbPattern(pattern const &P1)
{
  pattern P = SolidPattern(P1.A);
  P.funcPtrPatternAt = PerturbPatternAt;
  *P.ptrNext = P1;

  // ---
  // NOTE: Reset here so that it is possible to check for nullpointer
  //       in the actual pattern.
  // ---
  P.ptrNext->ptrNext.reset();

  // ---
  // NOTE: Check that the function pointers are not pointing to the Main pattern funcPtrPatternAt.
  // ---
  if (P1.funcPtrPatternAt == PerturbPatternAt)
  {
    std::cerr << "ERROR: " << __PRETTY_FUNCTION__ << ": Nested <pattern> functions would lead to infinite recursion"
              << std::endl;
    Assert(false, __FUNCTION__, __LINE__);
    return {};
  }

  return P;
}

/**
 * Return the color for the given pattern on  the given object, at the given
 * world-space point. Respects the transformations on both the pattern and the
 * object when doing so.
 * This function works by calling PatternAtShape. The Object referenced need
 * to have its material set up with a proper function to a pattern.
 */
tup StripeAtObject(pattern const &Pattern, shape const Object, tup const &Point)
{
  tup const Color = PatternAtShape(Pattern, Object, Point);
  return Color;
}

/**
 * Compute a Shape local point by calculating the inverse of the transform.
 * Use  the shape local point when calculating the Pattern point.
 * Then use  the function pointer to get the pattern at the specific point.
 * @Return: A tuple representing the color at the point.
 */
tup PatternAtShape(pattern const &Pattern, shape const &Shape, tup const &Point)
{
  Assert(Shape.Material.Pattern.funcPtrPatternAt != nullptr, __FILE__, __LINE__);
  /**
   * Multiply the given world-space point by the inverse of the objects
   * transformation matrix to convert the point to object space.
   */
  tup const ShapePoint = Inverse(Shape.Transform) * Point;

  tup const Color = Pattern.funcPtrPatternAt(Pattern, ShapePoint);
  return Color;
}

/**
 * Default function for pattern, used as initializer in CTOR.
 */
tup FuncDefaultPatternAt(pattern const &Pattern, tup const &ShapePoint)
{
  /**
   * Multiply the shape-space point by the inverse of the pattern's
   * transformation matrix to convert that point to the Pattern space.
   */
  tup const PatternPoint = Inverse(Pattern.Transform) * ShapePoint;

  tup const C = Color(PatternPoint.X, PatternPoint.Y, PatternPoint.Z);
  return C;
}

/**
 */
tup PatternAt(pattern const &Pattern, tup const &Point) { return Pattern.funcPtrPatternAt(Pattern, Point); }

/**
 * Return the color of the Pattern at the given Point.
 */
tup SolidPatternAt(pattern const &Pattern, tup const &ShapePoint)
{
  tup const Color = Pattern.A;
  return Color;
}

/**
 * Apply Perlin noise to the shape point prior the normal pattern computation.
 *
 * @Param Pattern: The actual pattern that will get the noise applied to.
 * @Param ShapePoint: World coordinate
 */
tup PerturbPatternAt(pattern const &Pattern, tup const &ShapePoint)
{
  tup const NoiseX = rtc::perlinnoise::improved_noise::Noise(ShapePoint.X, ShapePoint.X, ShapePoint.X) * .20f *
                     tup{1.f, 0.f, 0.f, 0.f};
  tup const NoiseY = rtc::perlinnoise::improved_noise::Noise(ShapePoint.Y, ShapePoint.Y, ShapePoint.Y) * 0.10f *
                     tup{0.f, 1.f, 0.f, 0.f};
  tup const NoiseZ = rtc::perlinnoise::improved_noise::Noise(ShapePoint.Z, ShapePoint.Z, ShapePoint.Z) * .20f *
                     tup{0.f, 0.f, 1.f, 0.f};
  tup const NoisyShapePoint = ShapePoint + NoiseX + NoiseY + NoiseZ;

  tup const Color = Pattern.ptrNext->funcPtrPatternAt(*Pattern.ptrNext, NoisyShapePoint);
  Assert(Pattern.Continue, __FUNCTION__, __LINE__);
  return Color;
}

/**
 * Return the color of the Pattern at the given Point.
 */
tup StripeAt(pattern const &Pattern, tup const &ShapePoint)
{
  /**
   * Multiply the shape-space point by the inverse of the pattern's
   * transformation matrix to convert that point to the Pattern space.
   */
  tup const PatternPoint = Inverse(Pattern.Transform) * ShapePoint;

  tup Color{Pattern.B};
  if (int(std::floorf(PatternPoint.X)) % 2 == 0) Color = Pattern.A;
  Assert(Pattern.Continue, __FUNCTION__, __LINE__);
  return Color;
}

/**
 * Return a blend of the two patterns at the given Point.
 */
tup BlendedPatternAt(pattern const &Pattern, tup const &ShapePoint)
{
  tup const Color1 = Pattern.ptrNext->funcPtrPatternAt(Pattern, ShapePoint);
  tup const Color2 = Pattern.ptrNext->ptrNext->funcPtrPatternAt(*Pattern.ptrNext, ShapePoint);
  tup const Color = 0.5f * Color1 + 0.5f * Color2;
  Assert(Pattern.Continue, __FUNCTION__, __LINE__);
  return Color;
}

/**
 * A Linear extrapolation pattern, LERP.
 * Calculate a distance beetween the two colors and use the Distance
 * in the X direction to compute a fraction that multiplied with the
 * Distance will produce a color in between the two colors of the
 * gradient.
 */
tup GradientPatternAt(pattern const &Pattern, tup const &ShapePoint)
{
  /**
   * Multiply the shape-space point by the inverse of the pattern's
   * transformation matrix to convert that point to the Pattern space.
   */
  tup const PatternPoint = Inverse(Pattern.Transform) * ShapePoint;

  pattern const &Gradient = Pattern;
  tup const Distance = Gradient.B - Gradient.A;
  float const Fraction = PatternPoint.X - std::floorf(PatternPoint.X);
  tup const Color = Gradient.A + Distance * Fraction;
  return Color;
}

/**
 * A ring pattern depends on the X and Z dimension.
 */
tup RingPatternAt(pattern const &Pattern, tup const &ShapePoint)
{
  /**
   * Multiply the shape-space point by the inverse of the pattern's
   * transformation matrix to convert that point to the Pattern space.
   */
  tup const PatternPoint = Inverse(Pattern.Transform) * ShapePoint;

  float const Hyp = std::sqrtf(PatternPoint.X * PatternPoint.X + PatternPoint.Z * PatternPoint.Z);
  int const Floor = int(std::floorf(Hyp)) % 2;
  tup const Color = Floor == 0 ? Pattern.A : Pattern.B;
  return Color;
}

/**
 * A radial/ring pattern with a gradient. Depends on the X and Z dimension.
 */
tup RadialGradientPatternAt(pattern const &Pattern, tup const &ShapePoint)
{
  /**
   * Multiply the shape-space point by the inverse of the pattern's
   * transformation matrix to convert that point to the Pattern space.
   */
  tup const PatternPoint = Inverse(Pattern.Transform) * ShapePoint;

  float const Hyp = std::sqrtf(PatternPoint.X * PatternPoint.X + PatternPoint.Z * PatternPoint.Z);
  int const Floor = int(std::floorf(Hyp)) % 2;
  tup const Color1 = Floor == 0 ? RingPatternAt(Pattern, PatternPoint) : GradientPatternAt(Pattern, PatternPoint);
  tup const Color2 =
      Pattern.ptrNext ? Pattern.ptrNext->funcPtrPatternAt(Pattern, PatternPoint) : ww::Color(0.f, 0.f, 0.f);
  tup const Color = Color1 + Color2;
  return Color;
}

/**
 * Get a pattern of alternating cubes by taking the sum of all directions mod 2.
 * Note: Point - In local coordinates.
 */
tup CheckersPatternAt(pattern const &Pattern, tup const &ShapePoint)
{
  /**
   * Multiply the shape-space point by the inverse of the pattern's
   * transformation matrix to convert that point to the Pattern space.
   */
  tup const PatternPoint = Inverse(Pattern.Transform) * ShapePoint;

  int const Floor =
      int(std::floorf(PatternPoint.X)) + int(std::floorf(PatternPoint.Y)) + int(std::floorf(PatternPoint.Z));

  // tup const Color = (0 == Floor % 2) ? Pattern.A : Pattern.B;

  tup const Color =
      (0 == Floor % 2)
          ? (Pattern.ptrNext ? Pattern.ptrNext->funcPtrPatternAt(*Pattern.ptrNext, PatternPoint) : Pattern.A)
          : (Pattern.ptrNext && Pattern.ptrNext->ptrNext
                 ? Pattern.ptrNext->ptrNext->funcPtrPatternAt(*Pattern.ptrNext->ptrNext, PatternPoint)
                 : Pattern.B);
  return Color;
}

/**
 * Get a pattern with gradient of alternating cubes by taking the sum of all directions mod 2.
 */
tup CheckersGradientPatternAt(pattern const &Pattern, tup const &ShapePoint)
{
  /**
   * Multiply the shape-space point by the inverse of the pattern's
   * transformation matrix to convert that point to the Pattern space.
   */
  tup const PatternPoint = Inverse(Pattern.Transform) * ShapePoint;

  int const Floor =
      int(std::floorf(PatternPoint.X)) + int(std::floorf(PatternPoint.Y)) + int(std::floorf(PatternPoint.Z));
  tup const Color =
      (0 == Floor % 2) ? CheckersPatternAt(Pattern, PatternPoint) : GradientPatternAt(Pattern, PatternPoint);
  return Color;
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
bool operator==(ww::pattern const &A, ww::pattern const &B) { return (ww::Equal(A.A, B.A) && ww::Equal(A.B, B.B)); }
bool operator!=(ww::pattern const &A, ww::pattern const &B) { return !(ww::Equal(A.A, B.A) && ww::Equal(A.B, B.B)); }
bool operator==(ww::sphere const &A, ww::sphere const &B) { return (ww::Equal(A, B)); }
bool operator==(ww::tup const &A, ww::tup const &B) { return (ww::Equal(A, B)); }
// ---
// NOTE: Division operator does not check for divide by zero; Who cares?
// ---
ww::tup operator/(ww::tup const &Tup, float const S) { return (ww::Mul(1.f / S, Tup)); }
