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
#define Assert(Condition, ...)                                        \
  if (!(Condition))                                                   \
  {                                                                   \
    fprintf(stderr, "ASSERT. Function %s. Line %d\n", ##__VA_ARGS__); \
    *(volatile int *)0 = 0;                                           \
  }
#else
#define Assert(Condition, ...)
#endif

namespace ww
{
constexpr float EPSILON = 0.0000100;  // 1E27 * std::numeric_limits<float>::min();
constexpr double PI = 3.141592653589793238463;
constexpr float PI_F = 3.14159265358979f;

//------------------------------------------------------------------------------
union tup {
  tup() : X{}, Y{}, Z{}, W{} {};
  tup(float IX, float IY, float IZ, float IW) : X{IX}, Y{IY}, Z{IZ}, W{IW} {};
  tup(float IX, float IY) : X{IX}, Y{IY}, Z{}, W{} {};
  tup(float IX, float IY, float IZ) : X{IX}, Y{IY}, Z{IZ}, W{} {};
  ~tup() {}
  tup(tup const &Other)
  {
    X = Other.X;
    Y = Other.Y;
    Z = Other.Z;
    W = Other.W;
  }
  struct  //!< A tuple is initially a vector with four elements or a 3D point.
  {
    float X{};
    float Y{};
    float Z{};
    float W{};  //!< is 1.0 when tuple is point and 0.0 when tuple is a vector.
  };
  struct  //!< A tuple can also represent colors with Intensity
  {
    float R;
    float G;
    float B;
    float I;  //!< Intensity is 1.0 at max and 0.0 at pitch black.
  };
  struct  //!< A tuple is a vector with four columns.
  {
    float C[4];
  };
};  // end of union tup.

/// ---
/// \struct base struct for the raytracing objects
/// ---
struct object
{
  tup Center{};
  object() {}
  virtual ~object() {}
  template <typename T>
  bool isA()
  {
    return (dynamic_cast<T *>(this) != NULL);
  }
};

/// ---
/// \struct sphere
/// \brief The sphere is defined by its center and the radii.
/// \detailed By making this struct a subclass of object we are able
///           to set up pointers to objeects of different types.
/// ---
struct sphere : public object
{
  float R{1.f};  //!< Radius.
};

/// ---
/// \struct cube
/// \brief Placeholder/stub for an upcoming cube. For now it is used for testing
///        of object pointers.
/// ---
struct cube : public object
{
  float L{1.f};
};

/// ---
/// \struct intersections
/// \brief Connect the time t value with the object for intersection.
struct intersection
{
  float t{};
  object *pObject{};  //!< The pointer need to be cast to a valid object type.
  object Object{};
};

/// ---
/// \struct intersection
/// \brief Collection of t values where a ray intersects an object.
/// ---
struct intersect
{
  std::vector<intersection> vI{};  //!< The points of intersection.
  int Count() const { return (int)vI.size(); }
};

/// ---
/// \struct intersections
/// \brief A collection of intersect's as defined above.
/// ---
struct intersections
{
  std::vector<intersection> vI{};
  int Count() const { return (int)vI.size(); }
};
/// ---
/// \struct ray
/// \brief A ray consist of an origint point and a vector for the direction.
/// ---
struct ray
{
  tup O{};  //!< The origin.
  tup D{};  //!< The direction.
};

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

// NOTE: Use a struct to return multiple values.
struct is_invertible_return
{
  bool IsInvertible{};
  bool IsComputed{};
  float Determinant{};
};

//------------------------------------------------------------------------------
union matrix {
  matrix() : R0{}, R1{}, R2{}, R3{}, Dimension{4} {}
  matrix(tup const &Cr0, tup const &Cr1, tup const &Cr2, tup const &Cr3)
      : R0{Cr0}, R1{Cr1}, R2{Cr2}, R3{Cr3}, Dimension{4}
  {
  }
  ~matrix() {}
  matrix(matrix const &Other)
  {
    R0 = Other.R0;
    R1 = Other.R1;
    R2 = Other.R2;
    R3 = Other.R3;
    Dimension = Other.Dimension;
    ID = Other.ID;
  }

  struct  //!< A matrix can be four tuple rows
  {
    tup R0{};
    tup R1{};
    tup R2{};
    tup R3{};
    int Dimension{4};
    // NOTE: Storeage of invertible and determinant
    is_invertible_return ID{};
  };
  struct  //!< or it can be an array of 4 tuple rows.
  {
    tup R[4];
  };
};

//------------------------------------------------------------------------------
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
tup Mul(float const S, tup const &Tup);
float Radians(float Deg);

// ---
// NOTE: Multiply is also called Hadamard or Schur product.
// ---
tup Mul(tup const A, tup const B);
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
void WriteToPPM(canvas const &Canvas, std::string const &Filename = "test.ppm");
int WriteToPPMFile(canvas const &Canvas, std::string const &Filename = "test.ppm");
std::shared_ptr<canvas> ReadFromPPM(std::string const &Filename = "test.ppm");

// ---
// NOTE: Matrix functions.
// ---

// ---
// \fn Cofactor A minor that possibly have had its sign changed.
// \return
// ---
float Cofactor(matrix const &M, int RemoveRow, int RemoveCol);
float Cofactor33(matrix const &M, int RemoveRow, int RemoveCol);
float Cofactor44(matrix const &M, int RemoveRow, int RemoveCol);

float Determinant22(matrix const &M);
float Determinant33(matrix const &M);
float Determinant44(matrix const &M);
float Determinant(matrix const &M);
bool Equal(matrix const &A, matrix const &B);
float Get(matrix const &M, int Row, int Col);

is_invertible_return IsInvertible(matrix const &M);

/// \fn Inverse Calculate the inverse of the matrix M.
///
/// \brief The inverse is not always possible to calculate. When inversion is
///        not possible the Zero matrix will be returned.
/// \return Inverse when possible, Zero matrix otherwise.
matrix Inverse(matrix const &M);

/// ---
/// \fn Identity matrix
/// \return Returs a 4x4 identity matrix.
/// ---
matrix I();

matrix Matrix44();
matrix Matrix44(tup const &R0, tup const &R1, tup const &R2, tup const &R3);
matrix Matrix33(tup const &R0, tup const &R1, tup const &R2);
matrix Matrix22(tup const &R0, tup const &R1);

/// ---
/// \fn Minor Calculate the determinant of a 2x2 submatrix.
/// \ brief Remove one row and column to be able to calculate the determinant
///         of a 2x2 matrix.
/// \param matrix &
/// \param RemoveRow - The row to be removed before computing the determinant.
/// \param RemoveCol - The column to be removed before computing the determinant.
/// \return Determinant of a submatrix of which a Row,Column have been removed.
float Minor(matrix const &M, int RemoveRow, int RemoveCol);

matrix Mul(matrix const &A, matrix const &B);
matrix Mul(matrix const &A, matrix const &B);
tup Mul(matrix const &A, tup const &T);
void Set(matrix &M, int Row, int Col, float Value);
matrix Transpose(matrix const &M);
matrix SubMatrix(matrix const &M, int RemoveRow, int RemoveCol);

matrix RotateX(float Alfa);
matrix RotateY(float Alfa);
matrix RotateZ(float Alfa);
matrix Scale(float X, float Y, float Z);
matrix Shearing(float Xy, float Xz, float Yx, float Yz, float Zx, float Zy);
matrix Translation(float X, float Y, float Z);
// Combine translation, rotate and scale in one single function.
matrix TranslateScaleRotate(                   //!<
    float TransX, float TransY, float TransZ,  //!< Translation is in m(?)
    float ScaleX, float ScaleY, float ScaleZ,  //!< Scale input is unitless.
    float AlfaX, float AlfaY, float AlfaZ      //!< Input rotation in radians.
);

/// ---
/// \fn Sphere releated functions
/// ---
intersect Intersect(sphere const &Sphere, ray const &Ray);

/// ---
/// Ray releated functions.
/// ---
ray Ray(tup const &P, tup const &V);
tup Position(ray const &R, float t);
intersection Intersection(float t, object *Object);
intersections Intersections(intersection const &I1, intersection const &I2);
intersect Intersect(sphere const &Object, ray const &Ray);
bool Equal(sphere const &A, sphere const &B);
};  // namespace ww

// ---
// NOTE: Declare the operator for inclusion elsewhere.
// ---
std::ostream &operator<<(std::ostream &stream, const ww::tup &T);
std::ostream &operator<<(std::ostream &stream, const ww::matrix &M);
ww::tup operator+(ww::tup const &A, ww::tup const &B);
ww::tup operator-(ww::tup const &Tup);
ww::tup operator-(ww::tup const &A, ww::tup const &B);
ww::tup operator*(float const S, ww::tup const &Tup);
ww::tup operator*(ww::tup const &Tup, float const S);
ww::tup operator*(ww::tup const &A, ww::tup const &B);
ww::tup operator/(ww::tup const &Tup, float const S);
ww::matrix operator*(ww::matrix const &A, ww::matrix const &B);
ww::tup operator*(ww::matrix const &A, ww::tup const &T);
bool operator==(ww::sphere const &A, ww::sphere const &B);
#endif
