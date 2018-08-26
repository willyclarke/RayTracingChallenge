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

#include <limits>

namespace ww
{
struct tup
{
   float X{};
   float Y{};
   float Z{};
   float W{};  //!< is 1.0 when tuple is point and 0.0 when tuple is a vector.
};

constexpr float EPSILON = 0.0000001;//1E27 * std::numeric_limits<float>::min();

// NOTE: Declarations.
bool Equal(float const A, float const B);
bool IsVector(tup const &Tup);
bool IsPoint(tup const &Tup);
tup Vector(float A, float B, float C);
tup Point(float A, float B, float C);

};  // namespace ww
#endif
