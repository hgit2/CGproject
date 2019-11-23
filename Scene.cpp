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


/*std::vector<Vector> Scene::inters_scene_point(Ray R){
	bool inters = this->inters_scene(R);
	if(not inters){
		// there is no intersection. We return a point behind the origin of the ray as a convention (assuming the camera is pointing toward the negative z)
		Vector vect(1, 1, 1);
		std::vector<Vector> result;
		result.push_back(R.v()+vect);
		result.push_back(R.v()+vect);
		return result;
	}
	else{
		// initialization
		Vector P; // intersection point
		Vector n; // normal to the sphere at P
		Sphere S; // sphere
		Vector c = S.m().c(); // parameter used to define the color of the sphere
		double t = 1000000000000; // min parameter corresponding to the closest intersection point
		// going through the different spheres
		for(int i = 0 ; i<N_ ; ++i){
			//std::cout << t << std::endl;
			S = sc_[i];
			Vector D = R.v()+(-S.v()); // C - O with C the origin of the ray and O the origin of sphere S
			double b = 2*(R.d().scalar_prod(D));
			double delta = pow(b,2) - 4*(D.sqnorm() - pow(S.r(), 2));			
			if(delta >= 0){
				double t1 = (-b-sqrt(delta))/2;
				double t2 = (-b+sqrt(delta))/2;
				//std::cout << "t1 " << t1 << std::endl;
				//std::cout << "t2 " << t2 << std::endl;	
				if(0<=t1 and t1<t){
					t = t1;
					P = R.v() + R.d()*t;					
					n = P+(-(S.v())); // normal vector to the sphere S at the point P				
					n = n.normalize(); // normalizing N
					c = S.m().c(); // getting the color of the sphere
					std::cout << "0<=t1 and t1<t : S " << S.v().x()<< " "  << S.v().y() << " " << S.v().z() << std::endl;
					cout << "P " << P.x() << " " << P.y() << " " << P.z() << endl;
					cout << "v_ " << R.v().x() << " " << R.v().y() << " " << R.v().z() << endl;
				}	
			
				else if(0 <= t2 and t2<t){
					t = t2;
					P = R.v() + R.d()*t;
					n = P+(-(S.v())); // normal vector to the sphere S at the point P				
					n = n.normalize(); // normalizing N
					c = S.m().c(); // getting the color of the sphere
					std::cout << "0 <= t2 and t2<t : S " << S.v().x() << S.v().y() << S.v().z() << std::endl;
					cout << "P " << P.x() << " " << P.y() << " " << P.z() << endl;
					cout << "v_ " << R.v().x() << " " << R.v().y() << " " << R.v().z() << endl;
				
				}
			}		
		}
		if(t==1000000000000){
			t=0;
		}		
		//vector containing the results
		std::vector<Vector> result;
		result.push_back(P);
		result.push_back(n);
		result.push_back(c);
		return result;
	}
}*/

/*std::vector<Vector> Scene::inters_scene_point(Ray R){
	bool inters = this->inters_scene(R);
	if(not inters){
		// there is no intersection. We return a point behind the origin of the ray as a convention (assuming the camera is pointing toward the negative z)
		Vector vect(1, 1, 1);
		std::vector<Vector> result;
		result.push_back(R.v()+vect);
		result.push_back(R.v()+vect);
		return result;
	}
	else{
		// initialization
		Vector P; // intersection point
		Vector n; // normal to the sphere at P
		Sphere S; // sphere
		Vector c = S.m().c(); // parameter used to define the color of the sphere
		double t = 1000000000000; // min parameter corresponding to the closest intersection point
		// going through the different spheres
		for(int i = 0 ; i<N_ ; ++i){
			S = sc_[i];
			bool inters = R.inters_sphere(S);
			if(inters){
				Vector D = R.v()+(-S.v()); // C - O with C the origin of the ray and O the origin of sphere S
				double b = 2*(R.d().scalar_prod(D));
				double delta = pow(b,2) - 4*(D.sqnorm() - pow(S.r(), 2));			
				if(delta >= 0){
					double t1 = (-b-sqrt(delta))/2;
					double t2 = (-b+sqrt(delta))/2;
					//std::cout << "t1 " << t1 << std::endl;
					//std::cout << "t2 " << t2 << std::endl;	
					if(0<=t1 and t1<t){
						t = t1;
						P = R.v() + R.d()*t;					
						n = P+(-(S.v())); // normal vector to the sphere S at the point P				
						n = n.normalize(); // normalizing N
						c = S.m().c(); // getting the color of the sphere
						//std::cout << "0<=t1 and t1<t : S " << S.v().x()<< " "  << S.v().y() << " " << S.v().z() << std::endl;
						//cout << "P " << P.x() << " " << P.y() << " " << P.z() << endl;
						//cout << "v_ " << R.v().x() << " " << R.v().y() << " " << R.v().z() << endl;
					}	
				
					else if(0 <= t2 and t2<t){
						t = t2;
						P = R.v() + R.d()*t;
						n = P+(-(S.v())); // normal vector to the sphere S at the point P				
						n = n.normalize(); // normalizing N
						c = S.m().c(); // getting the color of the sphere
						//std::cout << "0 <= t2 and t2<t : S " << S.v().x() << S.v().y() << S.v().z() << std::endl;
						//cout << "P " << P.x() << " " << P.y() << " " << P.z() << endl;
						//cout << "v_ " << R.v().x() << " " << R.v().y() << " " << R.v().z() << endl;
					
					}
				}
			}	
		}
		if(t==1000000000000){
			t=0;
		}		
		//vector containing the results
		std::vector<Vector> result;
		result.push_back(P);
		result.push_back(n);
		result.push_back(c);
		return result;
	}
}*/

/*std::vector<Vector> Scene::inters_scene_point(Ray R){
	bool inters = this->inters_scene(R);
	if(not inters){
		// there is no intersection. We return a point behind the origin of the ray as a convention (assuming the camera is pointing toward the negative z)
		Vector vect(1, 1, 1);
		std::vector<Vector> result;
		result.push_back(R.v()+vect);
		return result;
	}
	else{
		// there is at least one intersection between the ray and the scene
		//cout << "inters" << endl;
		// initialization
		Vector cmin; // parameter used to define the color of the sphere that is intersected
		Vector tmin; // type of the sphere that is intersected (diffuse or specular or transparent)
		Vector Nspmin; // refractive index of the sphere that is intersected
		Vector Pmin=R.v(); // the closest intersection point (initialized at the origin of the ray R)
		double Distmin=0; // distance between Pmin and the origin of the ray R
		Vector nmin; // normal vector to the sphere S at the point Pmin

		// going through the different spheres
		for(int i = 0 ; i<N_ ; ++i){
			Sphere S = sc_[i];
			Vector D = R.v()-S.v(); // C - O with C the origin of the ray and O the origin of sphere S
			double b = 2*(R.d().scalar_prod(D));
			double delta = pow(b,2) - 4*(D.sqnorm() - pow(S.r(), 2));			
			if(delta >= 0){
				if(delta==0){cout << "delta "<< delta<<endl;}
				// then there is an intersection between the ray and the sphere
				double t1 = (-b-sqrt(delta))/2;
				double t2 = (-b+sqrt(delta))/2;
				//std::cout << "t1 " << t1 << std::endl;
				//std::cout << "t2 " << t2 << std::endl;	
				if(0<t1){
					double Dist = (R.d()*t1).norm(); // distance between the intersection point R.v() + R.d()*t1 and the origin R.v()
					//cout << "dist t1: " << Dist << endl;
					if(Distmin==0 or Dist<Distmin){
						// the new intersection point is closer than the previous closest one
						Pmin = R.v() + R.d()*t1;
						Distmin = Dist;
						cmin = S.m().c(); // getting the color of the sphere
						Vector type(S.m().t(), S.m().t(), S.m().t()); // getting the sphere's type of material (diffuse or specular or transparent) 
						tmin = type;
						Vector Nsphere(S.m().Nsphere(), S.m().Nsphere(), S.m().Nsphere());
						Nspmin = Nsphere;
						nmin = Pmin+(-(S.v())); // normal vector to the sphere S at the point Pmin			
						nmin = nmin.normalize(); // normalizing N
					}			
				}	
				else if(0 < t2){
					double Dist = (R.d()*t2).norm(); // distance between the intersection point R.v() + R.d()*t2 and the origin R.v()
					//cout << "dist t2: " << Dist << "R.d() " << R.d().x() << " " << R.d().y() << " " << R.d().z() << endl;
					if((Distmin==0 or Dist<Distmin) and Dist != 0){
						// the new intersection point is closer than the previous closest one
						Pmin = R.v() + R.d()*t2; // updating Pmin
						Distmin = Dist; // updating Distmin
						cmin = S.m().c(); // getting the color of the sphere
						Vector type(S.m().t(), S.m().t(), S.m().t()); // getting the sphere's type of material (diffuse or specular or transparent) 
						tmin = type;
						//cout << "tmin " << tmin.x() << endl;
						Vector Nsphere(S.m().Nsphere(), S.m().Nsphere(), S.m().Nsphere()); // getting the sphere's refractive index
						Nspmin = Nsphere;
						nmin = Pmin+(-(S.v())); // updating nmin		
						nmin = nmin.normalize(); // normalizing nmin
					}			
				}
			//cout << "i " << i << "distmin: " << Distmin << endl;
			}	
		}
		//vector containing the results
		std::vector<Vector> result;
		//cout << "P " << Pmin.x() << " " << Pmin.y() << " " << Pmin.z() << endl;
		//cout << "v_ " << R.v().x() << " " << R.v().y() << " " << R.v().z() << endl;
		result.push_back(Pmin); // intersection point
		result.push_back(nmin); // normal to the surface
		result.push_back(cmin); // color of the surface
		result.push_back(tmin); // type of the surface (diffuse or specular or transparent)
		result.push_back(Nspmin); // refractive index of the surface
		return result;
	}
}*/

std::vector<Vector> Scene::inters_scene_point(Ray R){

	//// initialization
	Vector cmin; // parameter used to define the color of the sphere that is intersected
	Vector tmin; // type of the sphere that is intersected (diffuse or specular or transparent)
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
			if(0<t1){
				double Dist = (R.d()*t1).norm(); // distance between the intersection point R.v() + R.d()*t1 and the origin of the ray R.v()
				//cout << "dist t1: " << Dist << endl;
				if(Distmin==0 or Dist<Distmin){
					// the new intersection point is closer than the previous closest one
					Pmin = R.v() + R.d()*t1;
					Distmin = Dist;
					cmin = S.m().c(); // getting the color of the sphere
					Vector type(S.m().t(), S.m().t(), S.m().t()); // getting the sphere's type of material (diffuse or specular or transparent) 
					tmin = type;
					Vector Nsphere(S.m().Nsphere(), S.m().Nsphere(), S.m().Nsphere());
					Nspmin = Nsphere;
					nmin = Pmin+(-(S.v())); // normal vector to the sphere S at the point Pmin			
					nmin = nmin.normalize(); // normalizing N
				}			
			}	
			else if(0 < t2){
				double Dist = (R.d()*t2).norm(); // distance between the intersection point R.v() + R.d()*t2 and the origin R.v()
				//cout << "dist t2: " << Dist << "R.d() " << R.d().x() << " " << R.d().y() << " " << R.d().z() << endl;
				if((Distmin==0 or Dist<Distmin) and Dist != 0){
					// the new intersection point is closer than the previous closest one
					Pmin = R.v() + R.d()*t2; // updating Pmin
					Distmin = Dist; // updating Distmin
					cmin = S.m().c(); // getting the color of the sphere
					Vector type(S.m().t(), S.m().t(), S.m().t()); // getting the sphere's type of material (diffuse or specular or transparent) 
					tmin = type;
					//cout << "tmin " << tmin.x() << endl;
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
		//cout << "P " << Pmin.x() << " " << Pmin.y() << " " << Pmin.z() << endl;
		//cout << "v_ " << R.v().x() << " " << R.v().y() << " " << R.v().z() << endl;
		result.push_back(Pmin); // intersection point
		result.push_back(nmin); // normal to the surface
		result.push_back(cmin); // color of the surface
		result.push_back(tmin); // type of the surface (diffuse or specular or transparent)
		result.push_back(Nspmin); // refractive index of the surface
		return result;	
	}
	else{
		// then there is no intersection. We return a point behind the origin of the ray as a convention (assuming the camera is pointing toward the negative z)
		Vector vect(1, 1, 1);
		std::vector<Vector> result; // vector containing the result
		result.push_back(R.v()+vect);
		return result;
	}
}


Material Scene::get_color(Ray R, int nbrebound, Vector L, double I, double eps, double Nair){
	/* returns the color of the ray R leaving from the camera towards the scene */
	//cout << nbrebound << endl;
	// Initialization
	Vector vect0(0,0,0);	
	Material result(vect0); // black color

	// checks whether the ray R intersects the scene
	std::vector<Vector> inters = this->inters_scene_point(R);
	//cout << "n inters " << inters[1].x() << " " << inters[1].y() << " " << inters[1].z() << endl;
	//cout << "type inters " << inters[3].x() << endl;
	Vector vect(1,1,1);
	if(inters[0]==(R.v()+vect)){
		//then there is no intersections
		return result;
	}
	else if(inters[3].x() == 1 and nbrebound>0){
		// then the intersected surface is specular (inters[3] is the type of the intersected surface)
		//cout << "specular" << endl;
		result = get_color(R.reflect(inters[0], inters[1], eps), nbrebound-1, L, I, eps, Nair); // inters[0] is the intersection point and inters[1] is the normal 
	}
	else if(inters[3].x() == 2 and nbrebound>0){
		// then the intersected surface is transparent (inters[3] is the type of the intersected surface)
		//cout << "transparent " << "nbrebound " << nbrebound << endl;
		double Nsphere = inters[4].x();
		Ray refract = R.refract(inters[0], inters[1], Nair, Nsphere, eps);
		//cout << "refracted " << refract.d().x() <<  " " << refract.d().y() << " " << refract.d().z() << endl;
		result = get_color(R.refract(inters[0], inters[1], Nair, Nsphere, eps), nbrebound-1, L, I, eps, Nair); // inters[0] is the intersection point and inters[1] is the normal 
	}
	else{
		result.set_c(inters[2]); // color of the intersected surface
		// the intersected surface is diffuse

		// naming the intersection point between the ray R and the scene
		/*Vector P = inters[0]; // intersection point
		Vector n = inters[1]; // normal
		Vector c = inters[2]; // color of the sphere that is intersected
		Vector Peps = P + n*eps; // point from which we will compute the shadows
		Vector D = L+(-Peps); // vector going from Peps to the light source L
		double d = D.norm(); // distance between Peps and L
		Vector l = D.normalize(); // unitary vector going from Peps to L

		// creating shadows
		Ray rayP(Peps, l); // creation of a ray of origin the intersection point (shifted by epsilon) Peps and direction l
		bool intersbis = this->inters_scene(rayP); // intersection between rayP and the scene?
		if(intersbis){
			// rayP intersects the scene
			std::vector<Vector> resbis = this->inters_scene_point(rayP);
			Vector Pbis = resbis[0]; // intersection point between rayP and the scene
			Vector PPbis = Pbis-Peps; // vector going from Peps to Pbis
			double dbis = PPbis.norm(); // distance between Peps and Pbis
			if(dbis < d ){
				// there is a shadow
				Vector color(0., 0., 0.); // black
				result.set_c(result.c()+color); // updating the resulting color

			}
			else{
				// there is no shadow
				double intensity = std::max(0., std::min(l.scalar_prod(n)*I/pow( d, 2), 255.) )/3.14; // intensity of the pixel i,j
				Vector color(intensity*c.x(), intensity*c.y(), intensity*c.z());
				result.set_c(result.c()+color); // updating the resulting color
				//cout << "nbrebound " << nbrebound << "result " << result.c().x() << " " << result.c().y() << " " << result.c().z() << endl;
			}
		}
		else{
			// rayP does not intersect the scene, so there is no shadow
			double intensity = std::max(0., std::min(l.scalar_prod(n)*I/pow( d, 2), 255.) )/3.14; // intensity of the pixel i,j
			Vector color(intensity*c.x(), intensity*c.y(), intensity*c.z());
			result.set_c(result.c()+color); // updating the resulting color
		}
*/
	}
	return result;
}

Material Scene::get_color_brdf(Ray R, int nbrebound, Vector L, double I, double eps, double Nair){
	/* returns the color of the ray R leaving from the camera towards the scene */
	//cout << nbrebound << endl;
	// Initialization
	Vector vect0(0,0,0);	
	Material result(vect0); // black color

	// checks whether the ray R intersects the scene
	std::vector<Vector> inters = this->inters_scene_point(R);
	//cout << "n inters " << inters[1].x() << " " << inters[1].y() << " " << inters[1].z() << endl;
	//cout << "type inters " << inters[3].x() << endl;
	Vector vect(1,1,1);
	if(inters[0]==(R.v()+vect)){
		//then there is no intersections
		//cout << "no inters" << endl;
		return result;
	}
	else if(inters[3].x() == 1 and nbrebound>0){
		// then the intersected surface is specular (inters[3] is the type of the intersected surface)
		//cout << "specular" << endl;
		result = get_color_brdf(R.reflect(inters[0], inters[1], eps), nbrebound-1, L, I, eps, Nair); // inters[0] is the intersection point and inters[1] is the normal 
	}
	else if(inters[3].x() == 2 and nbrebound>0){
		// then the intersected surface is transparent (inters[3] is the type of the intersected surface)
		//cout << "transparent " << "nbrebound " << nbrebound << endl;
		double Nsphere = inters[4].x();
		Ray refract = R.refract(inters[0], inters[1], Nair, Nsphere, eps);
		//cout << "refracted " << refract.d().x() <<  " " << refract.d().y() << " " << refract.d().z() << endl;
		result = get_color_brdf(refract, nbrebound-1, L, I, eps, Nair); // inters[0] is the intersection point and inters[1] is the normal 
	}
	
	else{
		// the intersected surface is diffuse

//// direct lighting

		// naming the intersection point between the ray R and the scene
		/*Vector P = inters[0]; // intersection point
		Vector n = inters[1]; // normal
		Vector c = inters[2]; // color of the sphere that is intersected
		Vector Peps = P + n*eps; // point from which we will compute the shadows
		Vector D = L+(-Peps); // vector going from Peps to the light source L
		double d = D.norm(); // distance between Peps and L
		Vector l = D.normalize(); // unitary vector going from Peps to L

		// creating shadows
		Ray rayP(Peps, l); // creation of a ray of origin the intersection point (shifted by epsilon) Peps and direction l
		bool intersbis = this->inters_scene(rayP); // intersection between rayP and the scene?
		if(intersbis){
			// rayP intersects the scene
			std::vector<Vector> resbis = this->inters_scene_point(rayP);
			Vector Pbis = resbis[0]; // intersection point between rayP and the scene
			Vector PPbis = Pbis-Peps; // vector going from Peps to Pbis
			double dbis = PPbis.norm(); // distance between Peps and Pbis
			if(dbis < d ){
				// there is a shadow
				Vector color(0., 0., 0.); // black
				result.set_c(result.c()+color); // updating the resulting color
			}
			else{
				// there is no shadow
				double intensity = std::max(0., std::min(l.scalar_prod(n)*I/pow( d, 2), 255.) )/3.14; // intensity of the pixel i,j
				Vector color(intensity*c.x(), intensity*c.y(), intensity*c.z());
				result.set_c(result.c()+color); // updating the resulting color
				//cout << "nbrebound " << nbrebound << "result " << result.c().x() << " " << result.c().y() << " " << result.c().z() << endl;
			}
		}
		else{
			// rayP does not intersect the scene, so there is no shadow
			double intensity = std::max(0., std::min(l.scalar_prod(n)*I/pow( d, 2), 255.) )/3.14; // intensity of the pixel i,j
			Vector color(intensity*c.x(), intensity*c.y(), intensity*c.z());
			result.set_c(result.c()+color); // updating the resulting color
		}*/


		if(nbrebound>0){
			//// indirect lighting
		
			// sampling a random ray

			double r1 = distrib(engine);
			double r2 = distrib(engine);
			// random direction X, Y, Z coordinates
			double Xrand = cos(2*3.14*r1)*sqrt(1-r2);
			double Yrand = sin(2*3.14*r1)*sqrt(1-r2);
			double Zrand = sqrt(r2); 
			// changing coordinates system
			Vector A(distrib(engine), distrib(engine), distrib(engine)); // random vector
			Vector T1 = (inters[1].vector_prod(A)).normalize(); // vector perpendicular to the normal inters[1]
			Vector T2 = T1.vector_prod(inters[1]); // vector perpendicular to both the normal and T1
			Vector Drand = T1*Xrand + T2*Yrand + inters[1]*Zrand; // random direction
			Ray Rrand(inters[0]+inters[1]*eps, Drand); // random ray corresponding to the random direction Drand
	
			// computing the intensity (from the indirect lighting)
			result = get_color_brdf(Rrand,  nbrebound-1, L, I, eps, Nair);
//cout << "result resrqnd " << resRand.c().x() << " " << resRand.c().y() << " " << resRand.c().z() << endl;
		//cout << "result interm" << result.c().x() << " " << result.c().y() << " " << result.c().z() << endl;
	
	//		result.set_c(result.c() + resRand.c());
	//	cout << "result interm bis " << result.c().x() << " " << result.c().y() << " " << result.c().z() << endl;
		
		}
else{
	result.set_c(inters[2]); // color of the intersected surface


}	

		
	}


	//cout << "result " << result.c().x() << " " << result.c().y() << " " << result.c().z() << endl;
	return result;
}


// =====================================================================
//                                 Destructor
// =====================================================================
Scene::~Scene() = default;
