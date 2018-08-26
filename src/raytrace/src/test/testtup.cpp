/******************************************************************************
 * Filename : testtup.cpp
 * Date     : 2018 Aug 25
 * Author   : Willy Clarke (willy@clarke.no)
 * Version  : 0.0.1
 * Copyright: W. Clarke
 * License  : MIT
 * Descripti: Test file to test the datastructures in raytracer
 ******************************************************************************/
#include <datastructures.hpp>

#include <cstdlib>  //<! for size_t definition.
#include <iostream>

bool TestAssign(ww::tup Dest, ww::tup const &Src)
{
   Dest = Src;
   bool const Result = (Dest.X == Src.X);
   return (Result);
}

bool TestIsVector(ww::tup const& A)
{
    bool const Result = ww::IsVector(A);
    return (Result);
}

void RunTests()
{
   ww::tup A{};
   ww::tup B{};
   if (TestAssign(A, B))
   {
      std::cout << "TestAssign Success" << std::endl;
   }
   else
   {
       std::cerr << "TestAssign Failed!" << std::endl;
   }
   if (TestIsVector(A))
   {
       std::cout << "TestIsVector Success" << std::endl;
   }
   else
   {
       std::cerr << "TestIsVector Failed!" << std::endl;
   }
}

auto main() -> int
{
    return 0;
}
