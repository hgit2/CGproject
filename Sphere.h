#ifndef SPHERE_H_  
#define SPHERE_H_
#include "Vector.h"
#include "Material.h"

class Sphere{
	public:
	// =====================================================================
	//                               Constructors
	// =====================================================================
	Sphere(Vector v, double r); // creates a sphere of origin v and radius r
	Sphere(Vector v, double r, Material m); // creates a sphere of origin v and radius r and color defined by the material m
	Sphere(); // default constructor

	// =========================================================================
	//                                Destructor
	// =========================================================================
	virtual ~Sphere();

	// =========================================================================
	//                                Getters
	// =========================================================================
	inline  Vector v() const;
	inline  double r() const;
	inline  Material m() const;


	protected:
	// =====================================================================
	//                              Data members 
	// =====================================================================
	Vector v_; // the origin
	double r_;	// the direction
	Material m_; // the material (color)


};

// ===========================================================================
//                            Getters definitions
// ===========================================================================
inline Vector Sphere::v() const {
  return v_;
}

inline double Sphere::r() const {
  return r_;
}

inline Material Sphere::m() const {
  return m_;
}

#endif //SPHERE_H_

