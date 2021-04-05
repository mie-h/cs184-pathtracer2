#include "camera.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include "CGL/misc.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::ifstream;
using std::ofstream;

namespace CGL {

using Collada::CameraInfo;

Ray Camera::generate_ray_for_thin_lens(double x, double y, double rndR, double rndTheta) const {

  // TODO Project 3-2: Part 4
  // compute position and direction of ray from the input sensor sample coordinate.
  // Note: use rndR and rndTheta to uniformly sample a unit disk.


    //Look up your code from Project 3-1, Part 1 to figure out the generated ray direction (red segment in the figure).
    double radian_conv = 3.1415926534 / 180;
    double right_corner_width = tan(hFov / 2 * radian_conv);
    double right_corner_height = tan(vFov / 2 * radian_conv);
    double width = 2 * right_corner_width;
    double height = 2 * right_corner_height;

    double x_camera = x * width - right_corner_width;
    double y_camera = y * height - right_corner_height;
    double z_camera = -1;

    Vector3D red = Vector3D(x_camera, y_camera, z_camera);
    
    //Uniformly sample the disk representing the thin lens at ....
    //Calculate pFocus by intersecting the plane of focus with the red segment.
    auto pLens = Vector3D((this->lensRadius) * sqrt(rndR) * cos(rndTheta), (this->lensRadius) * sqrt(rndR) * sin(rndTheta), 0);
    auto pFocus =  red* this->focalDistance - pLens;

    //Calculate the ray that originates from pLens, and set its direction towards pFocus (blue segment in the figure).
    //Normalize the direction of the ray, perform the camera - to - world conversion for both its originand 
    //direction, add pos to the ray's origin, and set the near and far clips to be the same as in Project 3-1, Part 1.
    auto direction = this->c2w * pFocus;
    direction.normalize();

    auto ray = Ray(this->c2w * pLens + this->pos, direction, this->fClip);
    ray.min_t = this->nClip;
    return ray;

}


} // namespace CGL
