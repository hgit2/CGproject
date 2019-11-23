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
	Nsphere_=1;
}

Material::Material(Vector c, int t){
	c_=c;
	t_=t;
	Nsphere_=1;
  }

Material::Material(){
	Vector vect(0,0,0);
	c_=vect;
	t_=0;
	Nsphere_=1;
}


// =====================================================================
//                                 Destructor
// =====================================================================
Material::~Material() = default;
