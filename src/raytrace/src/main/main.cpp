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

// ---
//NOTE: The flag GTEST_ENABLED is normally set by the CMakeLists.txt file.
//      Google Test is disabled in order to be able to build with the build.linux file.
//      This functionality may be added later, but it may require Google Test installed
//      on this PC.
// ---
#ifdef GTEST_ENABLED
#include "featuretest.hpp"
#include "projectile.hpp"
#else
// NOTE: Stubs for the tests.
auto RunTupleTest(int argc, char *argv[]) -> int
{
   std::cout << argv[0];
   return argc;
}
namespace rtcch1
{
auto RunProjectileTest(int argc, char *argv[]) -> int
{
   std::cout << argv[0];
   return argc;
}
};  // namespace rtcch1
#endif

// ---
// NOTE: Main function.
// ---
auto main(int argc, char *argv[]) -> int
{
   bool DoTheTest{};
   bool DoTheProjectile{};

   if (argc > 1)
   {
      std::string const Argv1{argv[1]};
      if ("-t" == Argv1)
      {
         std::cout << "Do the test" << std::endl;
         DoTheTest = true;
      }
      else if ("-p" == Argv1)
      {
         DoTheProjectile = true;
      }
   }

   if (DoTheTest)
   {
      RunTupleTest(argc, argv);
   }
   else if (DoTheProjectile)
   {
      rtcch1::RunProjectileTest(argc, argv);
   }

   return 0;
}
