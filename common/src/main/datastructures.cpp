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
#include <iostream>

namespace ww
{
bool Equal(float const A, float const B)
{
    if (std::fabs(A -B) < EPSILON)
    {
        return true;
    }
    return false;
}

bool IsPoint(tup const &Tup)
{
   bool Result{};
   Result = (Tup.W != 0);
   return (Result);
}

bool IsVector(tup const &Tup)
{
   bool Result{};
   Result = !(Tup.W != 0);
   return (Result);
}

tup Point(float A, float B, float C)
{
   tup Result{A, B, C, 1.f};
   return (Result);
}

tup Vector(float A, float B, float C)
{
   tup Result{A, B, C, 0.f};
   return (Result);
}
};  // namespace ww
