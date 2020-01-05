#include <iostream>
#include <vector>
#include <math.h>
#include "Scene.h"
#include <random>


std::default_random_engine engine;
std::uniform_real_distribution <double> distrib(0,1);

using std::cout;
using std::endl;


// =====================================================================
//                               Constructors
// =====================================================================

Scene::Scene(std::vector<Sphere> Sc, int N): sc_(Sc), N_(N){
}

  
// ===========================================================================
//                           Public Function members
// ===========================================================================
void Scene::add_sphere(Sphere S){
	sc_.push_back(S);
	N_ ++; 
}


bool Scene::inters_scene(Ray R){
	for(std::vector<Sphere>::iterator it = sc_.begin() ; it != sc_.end(); ++it){
		bool inters = R.inters_sphere(*it);
		if(inters){
			// the ray intersects at least one sphere in the scene
			return true;
		}
	}
	return false;
}


std::vector<Vector> Scene::inters_scene_point(Ray R){

	//// initialization
	Vector cmin; // parameter used to define the color of the sphere that is intersected
	Vector typemin; // type of the sphere that is intersected (diffuse or specular or transparent)
	Vector Nspmin; // refractive index of the sphere that is intersected
	Vector Pmin=R.v(); // the closest intersection point (initialized at the origin of the ray R)
	double Distmin=0; // distance between Pmin and the origin of the ray R
	Vector nmin; // normal vector to the sphere S at the point Pmin

	//// going through the different spheres
	for(int i = 0 ; i<N_ ; ++i){
		Sphere S = sc_[i];
		Vector D = R.v()-S.v(); // C - O with C the origin of the ray and O the origin of sphere S
		double b = 2*(R.d().scalar_prod(D));
		double delta = pow(b,2) - 4*(D.sqnorm() - pow(S.r(), 2));			
		if(delta >= 0){
			// then there is an intersection between the ray and the sphere
			
			double t1 = (-b-sqrt(delta))/2;
			double t2 = (-b+sqrt(delta))/2;
			if(t1>0){
				double Dist = (R.d()*t1).norm(); // distance between the intersection point R.v() + R.d()*t1 and the origin of the ray R.v()
				if(Distmin==0 or Dist<Distmin){
					// the new intersection point is closer than the previous closest one

					Pmin = R.v() + R.d()*t1;
					Distmin = Dist;
					cmin = S.m().c(); // getting the color of the sphere
					Vector type(S.m().t(), S.m().t(), S.m().t()); // getting the sphere's type of material (diffuse or specular or transparent) 
					typemin = type;
					Vector Nsphere(S.m().Nsphere(), S.m().Nsphere(), S.m().Nsphere());
					Nspmin = Nsphere;
					nmin = Pmin-(S.v()); // normal vector to the sphere S at the point Pmin			
					nmin = nmin.normalize(); // normalizing nmin
				}			
			}	
			else if(t2>0){
				double Dist = (R.d()*t2).norm(); // distance between the intersection point R.v() + R.d()*t2 and the origin R.v()
				if((Distmin==0 or Dist<Distmin) and Dist != 0){
					// the new intersection point is closer than the previous closest one

					Pmin = R.v() + R.d()*t2; // updating Pmin
					Distmin = Dist; // updating Distmin
					cmin = S.m().c(); // getting the color of the sphere
					Vector type(S.m().t(), S.m().t(), S.m().t()); // getting the sphere's type of material (diffuse or specular or transparent) 
					typemin = type;
					Vector Nsphere(S.m().Nsphere(), S.m().Nsphere(), S.m().Nsphere()); // getting the sphere's refractive index
					Nspmin = Nsphere;
					nmin = Pmin+(-(S.v())); // updating nmin		
					nmin = nmin.normalize(); // normalizing nmin
				}			
			}
		}
	}

	//// check whether there is at least one intersection
	if(Distmin != 0){
		// there is at least one intersection between the ray R and the scene

		//vector containing the results
		std::vector<Vector> result;
		result.push_back(Pmin); // intersection point
		result.push_back(nmin); // normal to the surface
		result.push_back(cmin); // color of the surface
		result.push_back(typemin); // type of the surface (diffuse or specular or transparent)
		result.push_back(Nspmin); // refractive index of the surface
		return result;	
	}
	else{
		// then there is no intersection.
		// a point behind the origin of the ray is returned as a convention (assuming the camera is pointing toward the negative z)
		Vector vect(0, 0, 10);
		std::vector<Vector> result; // vector containing the result
		result.push_back(R.v()+vect);
		return result;
	}
}

Material Scene::get_color(Ray R, int nbrebound, Vector L, double I, double eps, double Nair){
	// Initialization
	Vector vect0(0,0,0);	
	Material result(vect0); // black color

	// calculating where the ray R intersects the scene
	std::vector<Vector> inters = this->inters_scene_point(R);
	// naming the intersection point between the ray R and the scene
	Vector P = inters[0]; // intersection point
	Vector n = inters[1]; // normal
	Vector type = inters[3]; // type of the intersected sphere
		
	Vector vect(0,0,10);
	if(P == (R.v()+vect)){
		//then there is no intersections
		return result;
	}
	else if(type.x() == 1 and nbrebound>0){
		// then the intersected surface is specular (inters[3] is the type of the intersected surface)
		result = result + get_color(R.reflect(P, n, eps), nbrebound-1, L, I, eps, Nair); // inters[1] is the normal 
	}
	else if(type.x() == 2 and nbrebound>0){
		// then the intersected surface is transparent (inters[3] is the type of the intersected surface)
		double Nsphere = inters[4].x(); // refractive index of the intersected sphere
		Ray refract = R.refract(P, n, Nair, Nsphere, eps);
		result = result + get_color(R.refract(P, n, Nair, Nsphere, eps), nbrebound-1, L, I, eps, Nair); 
	}
	else{
		// the intersected surface is diffuse or nbrebound=0

		Vector c = inters[2]; // color of the sphere that is intersected
		Vector Peps = P + n*eps; // point from which we will compute the shadows
		Vector D = L-Peps; // vector going from Peps to the light source L
		double d = D.sqnorm(); // distance between Peps and L squared
		Vector l = D.normalize(); // unitary vector going from Peps to L
		double intensity; // intensity of the pixel
		const double pi = 3.14159265358979323846; // define pi

		//creating shadows
		Ray rayP(Peps, l); // creation of a ray of origin the intersection point (shifted by eps) Peps and direction l
		std::vector<Vector> resbis = this->inters_scene_point(rayP);
		Vector Pbis = resbis[0]; // intersection point between rayP and the scene
		if(Pbis==(Peps+vect)){
			// rayP does not intersect the scene, so there is no shadow
			intensity =  std::max(0.,l.scalar_prod(n)*I/d); 
		}
		else{
			// rayP intersects the scene
			Vector PPbis = Pbis-Peps; // vector going from Peps to Pbis
			double dbis = PPbis.sqnorm(); // distance between Peps and Pbis squared
			if(dbis < d ){
				// there is a shadow
				intensity = 0.;					
			}
			else{
				// there is no shadow
				intensity =  std::max(0.,l.scalar_prod(n)*I/d);
			}
		}
		Vector color(intensity*c.x()/pi, intensity*c.y()/pi, intensity*c.z()/pi);
		result.set_intensity(color); // updating the resulting color
	}

	return result;
}


Material Scene::get_color_brdf_param(Ray R, int nbrebound, Vector L, double I, double eps, double Nair, int Kmax){
	// Initialization
	Vector vect0(0,0,0);	
	Material result(vect0); // black color
	const double pi = 3.14159265358979323846; // define pi
	const double invpi = 1./pi; // define 1/pi

	// calculating where the ray R intersects the scene
	std::vector<Vector> inters = this->inters_scene_point(R);

	// naming the intersection point between the ray R and the scene
	Vector P = inters[0]; // intersection point
	Vector n = inters[1]; // normal
	Vector type = inters[3]; // type of the intersected sphere
		
	Vector vect(0,0,10);
	if(P == (R.v()+vect)){
		//then there is no intersections
		return result;
	}
	else if(type.x() == 1 and nbrebound>0){
		// then the intersected surface is specular 
		result = get_color_brdf_param(R.reflect(P, n, eps), nbrebound-1, L, I, eps, Nair, Kmax); // recursive call 
	}
	else if(type.x() == 2 and nbrebound>0){
		// then the intersected surface is transparent 
		double Nsphere = inters[4].x(); // refractive index of the intersected sphere
		Ray refract = R.refract(P, n, Nair, Nsphere, eps);
		result = get_color_brdf_param(R.refract(P, n, Nair, Nsphere, eps), nbrebound-1, L, I, eps, Nair, Kmax); // recursive call
	}
	else{
		// the intersected surface is diffuse or nbrebound=0

		//// direct lighting

		Vector c = inters[2]; // color of the sphere that is intersected
		Vector Peps = P + n*eps; // point from which we will compute the shadows
		Vector D = L-P; // vector going from Peps to the light source L
		double d = D.sqnorm(); // distance between Peps and L squared
		Vector l = D.normalize(); // unitary vector going from Peps to L
		double intensity; // intensity of the pixel

		//creating shadows
		Ray rayP(Peps, l); // creation of a ray of origin the intersection point (shifted by eps) Peps and direction l
		std::vector<Vector> resbis = this->inters_scene_point(rayP);
		Vector Pbis = resbis[0]; // intersection point between rayP and the scene
		if(Pbis==(Peps+vect)){
			// rayP does not intersect the scene, so there is no shadow
			intensity = std::max(0.,l.scalar_prod(n)*I/d); 
		}
		else{
			// rayP intersects the scene
			Vector PPbis = Pbis-Peps; // vector going from Peps to Pbis
			double dbis = PPbis.sqnorm(); // distance between Peps and Pbis squared
			if(dbis < d ){
				// there is a shadow
				intensity = 0.;		
			}
			else{
				// there is no shadow
				intensity = std::max(0., l.scalar_prod(n)*I/d);
			}
		}
		result.set_intensity(c*intensity*invpi); // updating the resulting intensity 

		if(nbrebound>0 and Kmax>1){
			//// indirect lighting		

	
			// sampling a random ray

			double r1 = distrib(engine);
			double r2 = distrib(engine);
			// random direction X, Y, Z coordinates
			double Xrand = cos(2*pi*r1)*sqrt(1-r2);
			double Yrand = sin(2*pi*r1)*sqrt(1-r2);
			double Zrand = sqrt(r2); 
			// changing coordinates system
			Vector A(distrib(engine)-.5, distrib(engine)-.5, distrib(engine)-.5); // random vector
			Vector T1 = (n.vector_prod(A)).normalize(); // vector perpendicular to the normal n
			Vector T2 = (T1.vector_prod(n)).normalize(); // vector perpendicular to both the normal and T1
			Vector Drand = T1*Xrand + T2*Yrand + n*Zrand; // random direction
			Ray Rrand(Peps, Drand); // random ray corresponding to the random direction Drand (origin = intersection point shifted by eps)

			Material resRand = get_color_brdf_param(Rrand, nbrebound-1, L, I, eps, Nair, Kmax); // recursive call


			// direct + indirect lighting (we assume no surface is emissive)

			result.set_intensity(result.intensity() + resRand.intensity()*c); // update the resulting intensity
		}
	}
	return result;
}


// =====================================================================
//                                 Destructor
// =====================================================================
Scene::~Scene() = default;
