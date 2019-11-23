#ifndef MATERIAL_H_  
#define MATERIAL_H_
#include "Vector.h"

class Material{
	public:
	// =====================================================================
	//                               Constructors
	// =====================================================================
	Material(Vector c); // creates a Material that corresponds to the color defined by the vector c and that is diffuse
	Material(Vector c, int t); // creates a Material corresponding to the color defined by the vector c and the surface type t (diffuse=0, specular=1, transparent=2)
	Material(Vector c, int t, double Nsphere); // creates a Material corresponding to the color defined by the vector c, the surface type t (diffuse=0, specular=1, transparent=2) and the refractive index Nsphere
	Material(); // default constructor (black diffuse surface with a refractive index corresponding to air)

	// =========================================================================
	//                                Destructor
	// =========================================================================
	virtual ~Material();

	// =========================================================================
	//                                Getters
	// =========================================================================
	inline  Vector c() const;
	inline  int t() const;
	inline  double Nsphere() const;

	// =========================================================================
	//                               Setters
	// =========================================================================
	inline  void set_c(Vector c); 


	protected:
	// =====================================================================
	//                              Data members 
	// =====================================================================
	Vector c_; // the color vector
	int t_; // type of surface (diffuse (t_=0) or specular (t_=1) or transparent (t_=2))
	double Nsphere_;
	
};

// ===========================================================================
//                            Getters definitions
// ===========================================================================
inline Vector Material::c() const {
	return c_;
}

inline int Material::t() const {
	return t_;
}

inline double Material::Nsphere() const {
	return Nsphere_;
}

// ===========================================================================
//                            Setters definitions
// ===========================================================================
inline void Material::set_c(Vector c){
	// to ensure that none of the coefficients are negative or greater than 255
	Vector vect( std::max(0., std::min(c.x(), 255.)), std::max(0., std::min(c.y(), 255.)), std::max(0., std::min(c.z(), 255.)) );
	this->c_ = vect;
}

#endif //MATERIAL_H_

