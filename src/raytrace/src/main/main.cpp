/******************************************************************************
 * Filename : main.cpp
 * Date     : 2018 Aug 25
 * Author   : Willy Clarke (willy@clarke.no)
 * Version  : 0.0.1
 * Copyright: W. Clarke
 * License  : MIT
 * Descripti: Test program for my Raytracing challenge.
 ******************************************************************************/
#include <datastructures.hpp>

#include <iostream>

#include "featuretest.hpp"
#include "gtest/gtest.h"

//TEST(ZeroValue, Construction)
//{
//    EXPECT_EQ(0, int(0));
//}

// ---
// NOTE: Main function.
// ---
auto main(int argc, char *argv[]) -> int
{
   bool DoTheTest{};

   if (argc > 1)
   {
      std::string const Argv1{argv[1]};
      if ("-t" == Argv1)
      {
         std::cout << "Do the test" << std::endl;
         DoTheTest = true;
      }
   }

   if (DoTheTest)
   {
      RunTupleTest(argc, argv);
   }

   return 0;
}
