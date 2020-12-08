/******************************************************************************
 * Filename : testsmatrix.hpp
 * Date     : 2018 Sep 15
 * Author   : Willy Clarke (willy@clarke.no)
 * Version  : 0.0.1
 * Copyright: W. Clarke
 * License  : MIT
 * Descripti: Google tests for chapter 3 of the Ray Tracing Challenge.
 ******************************************************************************/
#ifndef SRC_RAYTRACE_SRC_MAIN_TESTSMATRIX_HPP
#define SRC_RAYTRACE_SRC_MAIN_TESTSMATRIX_HPP

#include <datastructures.hpp>

#include <cmath>
#include <iostream>
#include <memory>  // for shared pointer.

#include "gtest/gtest.h"

namespace rtcch3
{
//#include "testmatrixch1.hpp"
//#include "testmatrixch6.hpp"
//#include "testmatrixch7.hpp"
#include "testmatrixch8.hpp"
#include "testmatrixch9.hpp"

//------------------------------------------------------------------------------
void RunMatrixTest(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  int Count{};
  for (size_t Idx = 0;  //<!
       Idx < 1;         //<!
       ++Idx)
  {
    Count++;
    int const Result = RUN_ALL_TESTS();
    std::cout << "Testing result was " << (Result ? "Failure." : "Success.")  //<!
              << ". Test count " << Count << std::endl;
    if (Result != 0) break;
  }
}
};  // namespace rtcch3
#endif
