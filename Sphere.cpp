#include <iostream>
#include "Sphere.h"
#include <math.h> 

// =====================================================================
//                               Constructors
// =====================================================================

Sphere::Sphere(Vector v, double r){
	v_ = v;
	r_ = r;
	Vector vect(0,0,0);
	Material m1(vect);
	m_ = m1;
}

Sphere::Sphere(Vector v, double r, Material m): v_(v), r_(r), m_(m){
}

Sphere::Sphere(){
	Vector vect(0,0,0);
	v_=vect;
	r_=1.;
	Material m1(vect);
	m_=m1;
}


// =====================================================================
//                                 Destructor
// =====================================================================
Sphere::~Sphere() = default;
