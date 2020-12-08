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
#include "testsmatrix.hpp"

namespace
{
// ---
// Simple help message
// ---
void PrintHelp()
{
  // NOTE: For Anis escape sequences:
  // https://stackoverflow.com/questions/4842424/list-of-ansi-color-escape-sequences
  std::cout << "\nCommand line switches:"                                     //!<
               "\n--help        : \033[32;1mShow help\033[0m"                 //!<
               "\n--test        : \033[32;1mRun test\033[0m"                  //!<
               "\n--projectile  : \033[32;1mRun the projectile tests\033[0m"  //!<
               "\n--test-matrix : \033[32;1mRun the Matrix test\033[0m"       //!<
            << std::endl;
}
};  // namespace
// ---
// NOTE: Main function.
// ---
auto main(int argc, char *argv[]) -> int
{
  if (argc > 1)
  {
    std::string const Argv1{argv[1]};
    if ("--help" == Argv1)
    {
      PrintHelp();
    }
    else if ("--test" == Argv1)
    {
      RunTupleTest(argc, argv);
    }
    else if ("--projectile" == Argv1)
    {
      // rtcch1::RunProjectileTest(argc, argv);
      rtcch2::RunProjectileTest();
    }
    else if ("--test-matrix" == Argv1)
    {
      rtcch3::RunMatrixTest(argc, argv);
      //rtcch3::Ch9Planes_AssigningATransformation_Test();
    }
  }
  else
  {
    PrintHelp();
  }

  return 0;
}
