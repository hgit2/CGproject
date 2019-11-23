#ifndef RAY_H_  
#define RAY_H_
#include "Vector.h"
#include "Sphere.h"

class Ray{
	public:
	// =====================================================================
	//                               Constructors
	// =====================================================================
	Ray(Vector v, Vector d); // creates a ray of origin v and direction d
	Ray(Vector C, int i, int j, int H, int W, double fov); // creates a ray originating from a camera of origin C (with a field of view fov in radian) and that goes through the pixel (i,j) in an image of height H and width W
	Ray(); // default constructor (creates a ray of origin and direction (0,0,0) )

	// =========================================================================
	//                                Destructor
	// =========================================================================
	virtual ~Ray();

	// =========================================================================
	//                                	Getters
	// =========================================================================
	inline  Vector v() const; // origin
	inline Vector d() const; // direction

	// ===========================================================================
	//                           Public Function members
	// ===========================================================================
	bool inters_sphere(Sphere S); // returns whether an intersection point between the ray and the sphere S exists
	Vector inters_sphere_point(Sphere S); // calculates the first intersection point between the ray and the sphere S if it exists
	Ray reflect(Vector P, Vector n, double eps); // returns the ray reflected at the intersection point P, with n the normal to the surface 
	Ray refract(Vector P, Vector n, double Nair, double Nsphere, double eps); // returns the refracted ray at the intersection point P, with n the normal to the surface, Nair and Nsphere the refractive indices of the air and the intersected sphere
	
	protected:
	// =====================================================================
	//                              Data members 
	// =====================================================================
	Vector v_; // the origin
	Vector d_; // the direction


};

// ===========================================================================
//                            Getters definitions
// ===========================================================================
inline Vector Ray::v() const {
  return v_;
}


inline Vector Ray::d() const {
  return d_;
}


#endif //RAY_H_
