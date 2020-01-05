#ifndef MATERIAL_H_  
#define MATERIAL_H_
#include "Vector.h"

class Material{
	public:
	// =====================================================================
	//                               Constructors
	// =====================================================================

	// creates a Material that corresponds to the diffuse color defined by the vector c and that is diffuse
	Material(Vector c); 

	/**
	 creates a Material corresponding to the color defined by the vector c 
	 and the surface type t (diffuse=0, specular=1, transparent=2)
	*/
	Material(Vector c, int t); 

	/**
	 creates a Material corresponding to the color defined by the vector c, 
	 the surface type t (diffuse=0, specular=1, transparent=2) and the refractive index Nsphere
	*/
	Material(Vector c, int t, double Nsphere);

 	// default constructor (black diffuse surface with a refractive index corresponding to air)
	Material(); 


	// =========================================================================
	//                                Destructor
	// =========================================================================
	virtual ~Material();

	// =========================================================================
	//                                Getters
	// =========================================================================
	inline  Vector c() const;
	inline  Vector intensity() const;
	inline  int t() const;
	inline  double Nsphere() const;

	// =========================================================================
	//                               Setters
	// =========================================================================
	inline  void set_c(Vector c); 
	inline  void set_intensity(Vector c); 

	
	// =========================================================================
	//                                 Operators
	// =========================================================================
	
	// add two Material objects (i.e., add their colors)
	Material operator+(Material const& m);


	protected:
	// =====================================================================
	//                              Data members 
	// =====================================================================
	Vector c_; // the diffuse color vector
	Vector intensity_; // intensity of the material (meaningful only when considering a single pixel)
	int t_; // type of surface (diffuse (t_=0) or specular (t_=1) or transparent (t_=2))
	double Nsphere_; // refractive index
	
};

// ===========================================================================
//                            Getters definitions
// ===========================================================================
inline Vector Material::c() const {
	return c_;
}

inline Vector Material::intensity() const {
	return intensity_;
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
	// to ensure that none of the diffuse color coefficients are negative or greater than 1
	Vector vect( std::max(0., std::min(c.x(), 1.)), std::max(0., std::min(c.y(), 1.)), std::max(0., std::min(c.z(), 1.)) );
	this->c_ = vect;
}

inline void Material::set_intensity(Vector c){
	// to ensure that none of the intensity coefficients are negative or greater than 1
	Vector vect( std::max(0., std::min(c.x(), 1.)), std::max(0., std::min(c.y(), 1.)), std::max(0., std::min(c.z(), 1.)) );	
	this->intensity_ = vect;
}

#endif //MATERIAL_H_

