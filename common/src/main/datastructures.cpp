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

namespace ww
{
    bool IsPoint(tup const &Tup)
    {
        bool Result{};
        Result = !(Tup.W != 0);
        return (Result);
    }

    bool IsVector(tup const &Tup)
    {
        bool Result{};
        Result = (Tup.W != 0);
        return (Result);
    }
};
