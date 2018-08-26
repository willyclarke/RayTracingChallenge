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

namespace ww
{
    struct tup
    {
        float X{};
        float Y{};
        float Z{};
        float W{}; //!< is 1.0 when tuple is point and 0.0 when tuple is a vector.
    };

    // NOTE: Declarations.
    bool IsVector(tup const &Tup);
    bool IsPoint(tup const &Tup);
};
#endif
