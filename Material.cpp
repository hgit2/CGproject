#include <iostream>
#include "Material.h"

// =====================================================================
//                               Constructors
// =====================================================================

Material::Material(Vector c, int t, double Nsphere): c_(c), t_(t), Nsphere_(Nsphere){
  }

Material::Material(Vector c){
	c_=c;
	t_=0;
	Nsphere_=1.;
	Vector vect(0.,0.,0.);
	intensity_=vect;
}

Material::Material(Vector c, int t){
	c_=c;
	t_=t;
	Nsphere_=1.;
	Vector vect(0.,0.,0.);
	intensity_=vect;
  }

Material::Material(){
	Vector vect(0.,0.,0.);
	c_=vect;
	intensity_=vect;
	t_=0;
	Nsphere_=1.;
}

// =========================================================================
//                                 Operators
// =========================================================================

Material Material::operator+(Material const& m){
	Material sum;
	sum.set_c(c_ + m.c());
	sum.set_intensity(intensity_ + m.intensity());	
	return sum;
}


// =====================================================================
//                                 Destructor
// =====================================================================
Material::~Material() = default;
