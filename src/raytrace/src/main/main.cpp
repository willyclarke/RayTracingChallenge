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
#include "projectile.hpp"

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
    // rtcch1::RunProjectileTest(argc, argv);
    rtcch2::RunProjectileTest();
  }

  return 0;
}
