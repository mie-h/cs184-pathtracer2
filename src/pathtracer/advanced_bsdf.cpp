#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO Project 3-2: Part 1
  // Implement MirrorBSDF
    *pdf = 1.0;
    reflect(wo, wi);
    return reflectance / abs_cos_theta(*wi);
  return Vector3D();
}

void MirrorBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Mirror BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    ImGui::TreePop();
  }
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
  return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D h) {
  // TODO Project 3-2: Part 2
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.
  return 1.0;
}

Vector3D MicrofacetBSDF::F(const Vector3D wi) {
  // TODO Project 3-2: Part 2
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Vector3D.

  return Vector3D();
}

Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi) {
  // TODO Project 3-2: Part 2
  // Implement microfacet model here.

  return Vector3D();
}

Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Project 3-2: Part 2
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.

  *wi = cosineHemisphereSampler.get_sample(pdf);
  return MicrofacetBSDF::f(wo, *wi);
}

void MicrofacetBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Micofacet BSDF"))
  {
    DragDouble3("eta", &eta[0], 0.005);
    DragDouble3("K", &k[0], 0.005);
    DragDouble("alpha", &alpha, 0.005);
    ImGui::TreePop();
  }
}

// Refraction BSDF //

Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Project 3-2: Part 1
  // Implement RefractionBSDF
   
  // task 3
    *pdf = 1;
    if (!refract(wo, wi, ior)) {
        return Vector3D();
    }
    float eta;
    if (wo.z > 0) {
        eta = ior;
    }
    else {
        eta = 1.0 / ior;
    }
    return transmittance / abs_cos_theta(*wi) / ( eta * eta);





  return Vector3D();
}

void RefractionBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

// Glass BSDF //

Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO Project 3-2: Part 1
  // Compute Fresnel coefficient and either reflect or refract based on it.

  // compute Fresnel coefficient and use it as the probability of reflection
  // - Fundamentals of Computer Graphics page 305

    // total internal reflection
  if (!refract(wo, wi, ior)) {
      // assign reflection of wo to wi
      reflect(wo, wi);
      // set pdf to 1
      *pdf = 1;
      return reflectance / abs_cos_theta(*wi);
  }
  else {
      // compute schlick's R
      float n1 = 1.0;
      float n2 = ior;
      float R0 = pow(((n1 - n2) / (n1 + n2)), 2);
      float schlick_r = R0 + (1 - R0) * pow((1 - abs_cos_theta(wo)), 5);
      if (coin_flip(schlick_r)) {
          // assign ref wo to wi
          reflect(wo, wi);
          *pdf = 1;
          return schlick_r * reflectance / abs_cos_theta(*wi);
     }
      else {
          refract(wo, wi, ior);
          *pdf = 1 - schlick_r;

          float eta;
          if (wo.z > 0) {
              eta = ior;
          }
          else {
              eta = 1.0 / ior;
          }
          return (1.0 - schlick_r) * transmittance / abs_cos_theta(*wi) / (eta * eta);

          
      }
  }
  
  return Vector3D();
}

void GlassBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

void BSDF::reflect(const Vector3D wo, Vector3D* wi) {

  // TODO Project 3-2: Part 1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
    *wi = Vector3D(-wo[0], -wo[1], wo[2]);

}

bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior) {

  // TODO Project 3-2: Part 1
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
    float n;
    if (wo.z > 0) {
        n = 1 / ior;
    }
    else {
        n = ior;
    }

    float inner = 1 - pow(n, 2) * (1 - pow(wo.z, 2));
    if (inner < 0) {
        return false;
    }
    float new_x = -1 * n * wo.x;
    float new_y = -1 * n * wo.y;
    float new_z;
    if (wo.z > 0) {
        new_z = sqrt(inner) * -1;
    }
    else {
        new_z = sqrt(inner);
    }
    *wi = Vector3D(new_x, new_y, new_z);



  return true;

}

} // namespace CGL
