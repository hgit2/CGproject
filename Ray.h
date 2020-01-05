#ifndef RAY_H_  
#define RAY_H_
#include "Vector.h"
#include "Sphere.h"

class Ray{
	public:
	// =====================================================================
	//                               Constructors
	// =====================================================================

	// creates a ray of origin v and direction d
	Ray(Vector v, Vector d); 

	/** 
	 creates a ray coming from a camera of origin C (with a field of view fov in radian)
	 and that goes through the pixel (i,j) in an image of height H and width W 
	*/
	Ray(Vector C, int i, int j, int H, int W, double fov);

	// default constructor (creates a ray of origin (0,0,0) and direction (1/3,1/3,1/3) )
	Ray(); 

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
	
	/**
	 returns true if an intersection point between the ray and the sphere S exists.
	 The norm of the ray's direction must be equal to 1 
	*/
	bool inters_sphere(Sphere S); 

	/**
	 calculates the first intersection point between the ray and the sphere S if it exists.
	 The norm of the ray's direction must be equal to 1 
	*/
	Vector inters_sphere_point(Sphere S); 

	/**
	 returns the ray reflected at the intersection point P, with n the normal to the surface.
	 The norm of the ray's direction must be equal to 1 
	*/
	Ray reflect(Vector P, Vector n, double eps); 

	/**
	 returns the refracted ray at the intersection point P, with n the normal to the surface,
	 Nair and Nsphere the refractive indices of the air and the intersected sphere.
	 The norm of the ray's direction must be equal to 1
	*/
	Ray refract(Vector P, Vector n, double Nair, double Nsphere, double eps); 


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
