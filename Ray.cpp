#include <iostream>
#include <math.h> 
#include "Ray.h"

using std::cout;
using std::endl;


// =====================================================================
//                               Constructors
// =====================================================================

Ray::Ray(Vector v, Vector d): v_(v), d_(d){
}

Ray::Ray(Vector C, int i, int j, int H, int W, double fov){
	v_ = C; // origin (= the camera)
	Vector d(j - W/2 + .5, i - H/2 + .5, -H/(2*tan(fov/2.)) ); // direction
	d_ = d*(1/d.norm()); // normalize
}

Ray::Ray(){
	Vector vect(0, 0, 0);
	v_ = vect; // origin 
	Vector vectd(1/3, 1/3, 1/3);
	d_ = vectd; // direction
}

// ===========================================================================
//                           Public Function members
// ===========================================================================

bool Ray::inters_sphere(Sphere S){
	Vector D = v_ - S.v(); // C - O with C the origin of the ray and O the origin of sphere S
	double b = 2*(d_.scalar_prod(D)); // 2 times the scalar product between the direction of the ray and D
	double delta = pow(b,2) - 4*(D.sqnorm() - pow(S.r(), 2));
	if(delta < 0){
		// there is no intersection
		return false;  
	}
	else{
		double t1;
		double t2;
		t1 = (-b-sqrt(delta))/2;
		t2 = (-b+sqrt(delta))/2;
		if(delta = 0 or t1>0 or t2>0){
			// there is an intersection 
			return true;	
		}
		else{
			// there is no intersection
			return false; 
		}
	}
}

Vector Ray::inters_sphere_point(Sphere S){
	Vector D = v_ - S.v(); // C - O with C the origin of the ray and O the origin of sphere S
	double b = 2*(d_.scalar_prod(D)); // 2 times the scalar product between the direction of the ray and D
	double delta = pow(b,2) - 4*(D.sqnorm() - pow(S.r(), 2));
	if(delta >= 0){
		double t1;
		double t2;
		Vector P;
		t1 = (-b-sqrt(delta))/2;
		t2 = (-b+sqrt(delta))/2;		
		if(delta = 0 or t1>0){
			// there is an intersection
			P = v_ + d_*t1;
			return P;		
		}
		else if(t2>0){
			// there is an intersection
			P = v_ + d_*t2;
			return P;
		}
		else{
			// there is no intersection
			// a point behind the origin of the ray is returned as a convention (assuming the camera is pointing toward the negative z)
			Vector vect(0, 0, 10);
			return v_+vect;
	
		}
	}
	else{		
		// there is no intersection. 
		// a point behind the origin of the ray is returned as a convention (assuming the camera is pointing toward the negative z)
		Vector vect(0, 0, 10);
		return v_+vect;
	}

}

/*Vector Ray::inters_sphere_point(Sphere S){
	bool inters = Ray::inters_sphere(S);
	if(not inters){		
		// there is no intersection. We return a point behind the origin of the ray as a convention (assuming the camera is pointing toward the negative z)
		Vector vect(1, 1, 1);
		cout << "not inters" << endl;
		return v_+vect;
	}
	else{
		Vector D = v_+(-S.v()); // C - O with C the origin of the ray and O the origin of sphere S
		double b = 2*(d_.scalar_prod(D));
		double delta = pow(b,2) - 4*(D.sqnorm() - pow(S.r(), 2));
		double t1;
		double t2;
		Vector P;
		t1 = (-b-sqrt(delta))/2;
		t2 = (-b+sqrt(delta))/2;		
		if(delta = 0 or t1>=0){
			P = v_ + d_*t1;
			return P;		
		}
		else if(t2>=0){
			P = v_ + d_*t2;
			return P;
		}
	}
}*/

Ray Ray::reflect(Vector P, Vector n, double eps){
	Vector r = d_ - n*d_.scalar_prod(n)*2; // direction of the reflected ray
	Vector Peps = P+n*eps; // P is shifted by eps to avoid noise
	Ray result(Peps, r); // reflected ray (origin = the intersection point P of the incident ray with the surface, shifted by eps)
	return result;
}


Ray Ray::refract(Vector P, Vector n, double Nair, double Nsphere, double eps){
	if(d_.scalar_prod(n) >= 0){
		n = -n; // inversion of the normal 
		double nair = Nair; // save the value of Nair
	
		// switch the roles of Nair and Nsphere
		Nair = Nsphere;
		Nsphere = nair;
	}
	double D = 1 - pow((Nair/Nsphere), 2)*(1 - pow(d_.scalar_prod(n), 2));
	if(D>=0){
		Vector r = d_*(Nair/Nsphere) - n*( d_.scalar_prod(n)*(Nair/Nsphere) + sqrt(D) ); // direction of the refracted ray
		Vector Peps = P-n*eps; // P is shifted by eps to avoid noise
		Ray result(Peps, r); // refracted ray (origin = the intersection point P of the incident ray with the surface, shifted by eps)
		return result;
	}
	else{
	// the ray is reflected
	Ray result = this->reflect(P, n, eps);
	return result;
	}
}


// =====================================================================
//                                 Destructor
// =====================================================================
Ray::~Ray() = default;



