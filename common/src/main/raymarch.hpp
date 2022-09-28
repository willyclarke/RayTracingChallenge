#ifndef COMMON_SRC_MAIN_RAYMARCH_HPP
#define COMMON_SRC_MAIN_RAYMARCH_HPP

#include "datastructures.hpp"

namespace ww
{

/**
 * rm - raymarch
 */
namespace rm
{

/**
 * RayMarch - move along ray to find hit.
 * @param: R - Ray
 * @param: PtrShape - Shape to compute distance against
 */
float RayMarch(ray const &R, shared_ptr_shape PtrShape);
float RayMarch(ray const &R, world const &World, shared_ptr_shape PtrShape);

/**
 */
tup MainImage(camera const &Camera, world const &World, int X, int Y, shared_ptr_shape PtrShape);

/**
 */
canvas Render(camera const &Camera, world const &World);

/**
 */
float GetDistance(tup const &P, shared_ptr_shape PtrShape);
float GetDistance(tup const &P, world const &World);

/**
 */
tup GetNormal(tup const &P, shared_ptr_shape PtrShape);

};      // end namespace rm
};      // end namespace ww
#endif  // COMMON_SRC_MAIN_RAYMARCH_HPP
