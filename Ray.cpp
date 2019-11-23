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
	d_ = d*(1/d.norm());
}

Ray::Ray(){
	Vector vect(0, 0, 0);
	v_ = vect; // origin 
	d_ = vect; // direction
}

// ===========================================================================
//                           Public Function members
// ===========================================================================
bool Ray::inters_sphere(Sphere S){
	//cout << "inters sphere" << endl;
	Vector D = v_+(-S.v()); // C - O with C the origin of the ray and O the origin of sphere S
	double b = 2*(d_.scalar_prod(D)); // scalar product between the direction of the ray and D
	double delta = pow(b,2) - 4*(D.sqnorm() - pow(S.r(), 2));
	//std::cout << "delta" << delta << std::endl;
	if(delta < 0){
		// there is no intersection
		return false;  
	}
	else{
		double t1;
		double t2;
		t1 = (-b-sqrt(delta))/2;
		t2 = (-b+sqrt(delta))/2;
		Vector P1 = v_ + d_*t1;
		Vector P2 = v_ + d_*t2;
		//cout << "P1 " << P1.x() << " " << P1.y() << " " << P1.z() << endl;
		//cout << "v_ " << v_.x() << " " << v_.y() << " " << v_.z() << endl;
		if(delta = 0 or t1>0){
		
//if(P1!=v_ and (delta = 0 or t1>0)){
			// the intersection point 
			//cout << "delta " << delta << ", t1 = " << t1 << endl; 
			return true;	
		}
		else if(t2>0){
		//else if(P2!=v_ and t2>0){
			//cout << "t2 " << t2 << endl;
			return true;
		}
		else{
			// there is no intersection
			//cout << "t1 " << t1 << " t2 " << t2 << endl;		
			return false; 
		}
	}
}

Vector Ray::inters_sphere_point(Sphere S){
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
}

Ray Ray::reflect(Vector P, Vector n, double eps){
	Vector r = d_ - n*d_.scalar_prod(n)*2; // direction of the reflected ray
	Vector Peps = P+n*eps; // P is shifted by eps to avoid noise
	Ray result(Peps, r); // reflected ray (origin = the intersection point P of the incident ray with the surface, shifted by eps)
	return result;
}


Ray Ray::refract(Vector P, Vector n, double Nair, double Nsphere, double eps){
	//cout << "n " << n.x() << " " << n.y() << " " << n.z() << " Nair " << Nair << " Nsphere " << Nsphere << endl;
	if(d_.scalar_prod(n) >= 0){
		//cout << "neg" << endl;
		n = -n; // inversion of the normal 
		double nair = Nair; // save the value of Nair
		// switch the roles of Nair and Nsphere
		Nair = Nsphere;
		Nsphere = nair;
	}
	double D = 1 - pow((Nair/Nsphere), 2)*(1 - pow(d_.scalar_prod(n), 2));
	//cout << "D " << D << endl; 
	if(D>=0){
		//cout << "refract" << endl;
		Vector r = d_*(Nair/Nsphere) - n*( d_.scalar_prod(n)*(Nair/Nsphere) + sqrt(D) ); // direction of the refracted ray
		Vector Peps = P-n*eps; // P is shifted by eps to avoid noise
		Ray result(Peps, r); // refracted ray (origin = the intersection point P of the incident ray with the surface, shifted by eps)
	//cout << "r " << r.x() << " " << r.y() << " " << r.z() << endl;
//cout << "result.d() " << result.d().x() << " " << result.d().y() << " " << result.d().z() << endl;

//cout << "P " << P.x() << " " << P.y() << " " << P.z() << endl;
//cout << "Peps " << Peps.x() << " " << Peps.y() << " " << Peps.z() << endl;
		return result;

	}
	else{
	// the ray is reflected
	//cout << "reflect" << endl;
	Ray result = this->reflect(P, n, eps);
	return result;
	}
}


// =====================================================================
//                                 Destructor
// =====================================================================
Ray::~Ray() = default;



