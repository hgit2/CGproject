#include <iostream>
#include <vector>
#include "Vector.h"
#include "Ray.h"
#include "Sphere.h"
#include "Scene.h"
#include "Material.h"
#include "template_for_saving_bmp_files.cpp"
#include <math.h>
#include <omp.h> // OpenMP

using std::cout;
using std::endl;

std::vector<unsigned char> create_bin_img(int W, int H, double fov, Vector C, Sphere S){
	// creates the image of a white disk
	std::vector<unsigned char> img(W*H * 3, 0);
	for(int i = 0; i < H; i++){
		for(int j = 0; j < W; j++){
			Ray rayij(C, i, j, H, W, fov);
			bool inters = rayij.inters_sphere(S);	
			if(inters){
				cout << inters << endl;
				// pixel at (i,j) is white
				img[(i * W + j)*3] = 255;  
				img[(i * W + j)*3+1] = 255;
				img[(i * W + j)*3+2] = 255;
			}				
			else{
				img[(i * W + j)*3] = 0;  
				img[(i * W + j)*3+1] = 0;
				img[(i * W + j)*3+2] = 0;
			}
		}
	}
	return img;
}

std::vector<unsigned char> create_vol_img(int W, int H, double fov, Vector C, Sphere S, Vector L, double I){
	/* Creates an image from the point of view of a camera located at C (with a field of view fov in radian)  
	of the sphere S illuminated by the light source located at L with an intensity I*/
	std::vector<unsigned char> img(W*H * 3, 0);
	for(int i = 0; i < H; i++){
		for(int j = 0; j < W; j++){
			Ray rayij(C, i, j, H, W, fov);
			bool inters = rayij.inters_sphere(S);
			if(inters){
				//cout << inters << endl;
				Vector P = rayij.inters_sphere_point(S); // intersection point between the ray rayij and the sphere S
				Vector D = L+(-P); // vector going from P to the light source L
				double d = D.norm(); // distance between P and L
				Vector l = D.normalize(); // unitary vector going from P to L				
				Vector N = P+(-S.v()); // normal vector to the sphere S at the point P				
				Vector n = N.normalize(); // normalizing N
				double intensity = std::max(0., std::min(l.scalar_prod(n)*I/pow( d, 2), 255.) ); // intensity of the pixel i,j
				cout << l.scalar_prod(n)*I/pow( d, 2) << endl;
				img[((H-i-1) * W + j)*3] =  intensity;
				img[((H-i-1) * W + j)*3+1] = intensity;
				img[((H-i-1) * W + j)*3+2] = intensity;
				//cout << l.scalar_prod(n)*I/pow( d, 2) << endl;
			}
			else{
				img[((H-i-1) * W + j)*3] = 0;  
				img[((H-i-1) * W + j)*3+1] = 0;
				img[((H-i-1) * W + j)*3+2] = 0;
			}
		}
	}
	return img;
}

std::vector<unsigned char> create_scene_img(int W, int H, double fov, Vector C, Scene Sc, Vector L, double I){
	/* Creates an image from the point of view of a camera located at C (with a field of view fov in radian)  
	of a Scene Sc (consisting of Spheres) illuminated by the light source located at L with an intensity I*/
	std::vector<unsigned char> img(W*H * 3, 0);
	for(int i = 0; i < H; i++){
		for(int j = 0; j < W; j++){
			Ray rayij(C, i, j, H, W, fov);
			bool inters = Sc.inters_scene(rayij);
			if(inters){
				//cout << inters << endl;
				std::vector<Vector> res = Sc.inters_scene_point(rayij); // intersection point between the ray rayij and the scene Sc and the corresponding normal
				Vector P = res[0]; // intersection point
				//cout << "P " << P.x() << endl;;
				Vector n = res[1]; // normal
				Vector c = res[2]; // color of the sphere
				Vector D = L+(-P); // vector going from P to the light source L
				double d = D.norm(); // distance between P and L
				Vector l = D.normalize(); // unitary vector going from P to L
				double intensity = std::max(0., std::min(l.scalar_prod(n)*I/pow( d, 2), 255.) ); // intensity of the pixel i,j
				img[((H-i-1) * W + j)*3] =  intensity*c.x();
				img[((H-i-1) * W + j)*3+1] = intensity*c.y();
				img[((H-i-1) * W + j)*3+2] = intensity*c.z();
				//cout << "l " << l.x() << endl;
				//cout << l.scalar_prod(n)*I/pow( d, 2) << endl;
			}
			else{
				img[((H-i-1) * W + j)*3] = 0;  
				img[((H-i-1) * W + j)*3+1] = 0;
				img[((H-i-1) * W + j)*3+2] = 0;
			}
		}
	}
	return img;
}

std::vector<unsigned char> create_scene_shadow_img(int W, int H, double fov, Vector C, Scene Sc, Vector L, double I, double eps){
	/* Creates an image, with shadows, from the point of view of a camera located at C (with a field of view fov in radian)  
		of a Scene Sc (consisting of Spheres) illuminated by the light source located at L with an intensity I.
		To calculate shadows, we will shift by eps the intersection point between the ray comming from the camera and a sphere.*/
	std::vector<unsigned char> img(W*H * 3, 0);
	for(int i = 0; i < H; i++){
		for(int j = 0; j < W; j++){
			Ray rayij(C, i, j, H, W, fov);
			bool inters = Sc.inters_scene(rayij);
			if(inters){
				//calculating the intersection point between rayij and the scene
				std::vector<Vector> res = Sc.inters_scene_point(rayij); // intersection point between rayij and Sc, the corresponding normal and "color"
				Vector P = res[0]; // intersection point
				Vector n = res[1]; // normal
				Vector c = res[2]; // color of the sphere that is intersected
				Vector Peps = P + n*eps; // point from which we will compute the shadows
				Vector D = L+(-Peps); // vector going from Peps to the light source L
				double d = D.norm(); // distance between Peps and L
				Vector l = D.normalize(); // unitary vector going from Peps to L

				//creating shadows
				Ray rayP(Peps, l); // creation of a ray of origin the intersection point (shifted by epsilon) Peps and direction l
				bool intersbis = Sc.inters_scene(rayP); // intersection between rayP and the scene?
				if(intersbis){
					// rayP intersects the scene
					std::vector<Vector> resbis = Sc.inters_scene_point(rayP);
					Vector Pbis = resbis[0]; // intersection point between rayP and the scene
					//Vector nbis = resbis[1]; // normal
					Vector PPbis = Pbis-Peps; // vector going from Peps to Pbis
					double dbis = PPbis.norm(); // distance between Peps and Pbis
					//cout << "d " << d << "dbis " << dbis << endl;
					//if(c.x()==0 and c.y()==1 and c.z()==0){cout<< "green" << endl;}
					//cout << resbis[2].x() << " " << resbis[2].y() << " " << resbis[2].z() << endl;

					if(dbis < d ){
						// there is a shadow
						//cout << c.x() << " " << c.y() << " " << c.z() << endl;
						//cout << "shadow" << endl;
						img[((H-i-1) * W + j)*3] = 0;  
						img[((H-i-1) * W + j)*3+1] = 0;
						img[((H-i-1) * W + j)*3+2] = 0;
					}
					else{
						// there is no shadow
						//cout << "no shadow" << endl;
						//if(c.x()==0 and c.y()==1 and c.z()==0){cout<< "green" << endl;}
						double intensity = std::max(0., std::min(l.scalar_prod(n)*I/pow( d, 2), 255.) ); // intensity of the pixel i,j
						img[((H-i-1) * W + j)*3] =  intensity*c.x();
						img[((H-i-1) * W + j)*3+1] = intensity*c.y();
						img[((H-i-1) * W + j)*3+2] = intensity*c.z();
					}
				}
				else{
					// rayP does not intersect the scene, so there is no shadow
					//cout << "no shadow" << endl;
					//if(c.x()==0 and c.y()==1 and c.z()==0){cout<< "green" << endl;}
					double intensity = std::max(0., std::min(l.scalar_prod(n)*I/pow( d, 2), 255.) ); // intensity of the pixel i,j
					img[((H-i-1) * W + j)*3] =  intensity*c.x();
					img[((H-i-1) * W + j)*3+1] = intensity*c.y();
					img[((H-i-1) * W + j)*3+2] = intensity*c.z();
				}
			}
			else{
				// no intersections between rayij and the scene
				img[((H-i-1) * W + j)*3] = 0;  
				img[((H-i-1) * W + j)*3+1] = 0;
				img[((H-i-1) * W + j)*3+2] = 0;
			}
		}
	}
	return img;
}

//std::vector<unsigned char> create_scene_shadow_specular_img(int W, int H, double fov, Vector C, Scene Sc, Vector L, double I, double eps, double Nair, double Kmax){
	/* Creates an image, with shadows, from the point of view of a camera located at C (with a field of view fov in radian)  
		of a Scene Sc (consisting of Spheres) illuminated by the light source located at L with an intensity I.
		To calculate shadows, we will shift by eps the intersection point between the ray comming from the camera and a sphere.*/
	/*std::vector<unsigned char> img(W*H * 3, 0);
	for(int i = 0; i < H; i++){
		for(int j = 0; j < W; j++){
			Ray rayij(C, i, j, H, W, fov);
			bool inters = Sc.inters_scene(rayij);
			if(inters){
				//calculating the intersection point between rayij and the scene
				std::vector<Vector> res = Sc.inters_scene_point(rayij); // intersection point between rayij and Sc, the corresponding normal and "color"
				Vector P = res[0]; // intersection point
				Vector n = res[1]; // normal
				Vector Peps = P + n*eps; // point from which we will compute the shadows
				Vector D = L+(-Peps); // vector going from Peps to the light source L
				double d = D.norm(); // distance between Peps and L
				Vector l = D.normalize(); // unitary vector going from Peps to L

				//creating shadows
				Ray rayP(Peps, l); // creation of a ray of origin the intersection point (shifted by epsilon) Peps and direction l
				bool intersbis = Sc.inters_scene(rayP); // intersection between rayP and the scene?
				if(intersbis){
					// rayP intersects the scene
					std::vector<Vector> resbis = Sc.inters_scene_point(rayP);
					Vector Pbis = resbis[0]; // intersection point between rayP and the scene
					//Vector nbis = resbis[1]; // normal
					Vector PPbis = Pbis-Peps; // vector going from Peps to Pbis
					double dbis = PPbis.norm(); // distance between Peps and Pbis
					//cout << "d " << d << "dbis " << dbis << endl;
					//if(c.x()==0 and c.y()==1 and c.z()==0){cout<< "green" << endl;}
					//cout << resbis[2].x() << " " << resbis[2].y() << " " << resbis[2].z() << endl;

					if(dbis < d ){
						// there is a shadow
						//cout << c.x() << " " << c.y() << " " << c.z() << endl;
						//cout << "shadow" << endl;
						img[((H-i-1) * W + j)*3] = 0;  
						img[((H-i-1) * W + j)*3+1] = 0;
						img[((H-i-1) * W + j)*3+2] = 0;
					}
					else{
						// there is no shadow
						//cout << "no shadow" << endl;
						//if(c.x()==0 and c.y()==1 and c.z()==0){cout<< "green" << endl;}
						double intensity = pow(std::max(0., std::min(l.scalar_prod(n)*I/pow( d, 2), 255.) ), 1.); // intensity of the pixel i,j
						Material color = Sc.get_color_brdf(rayij, 5, L, I, eps, Nair);
						//cout << color.c().x() << " " << color.c().y() << " " << color.c().z() << endl;
						img[((H-i-1) * W + j)*3] =  intensity*color.c().x();
						img[((H-i-1) * W + j)*3+1] = intensity*color.c().y();
						img[((H-i-1) * W + j)*3+2] = intensity*color.c().z();
					}
				}
				else{
					// rayP does not intersect the scene, so there is no shadow
					//cout << "no shadow" << endl;
					//if(c.x()==0 and c.y()==1 and c.z()==0){cout<< "green" << endl;}
					double intensity = pow(std::max(0., std::min(l.scalar_prod(n)*I/pow( d, 2), 255.) ), 1.); // intensity of the pixel i,j
					Material color = Sc.get_color_brdf(rayij, 5, L, I, eps, Nair);
					cout << color.c().x() << " " << color.c().y() << " " << color.c().z() << endl;					
					img[((H-i-1) * W + j)*3] =  intensity*color.c().x();
					img[((H-i-1) * W + j)*3+1] = intensity*color.c().y();
					img[((H-i-1) * W + j)*3+2] = intensity*color.c().z();

				}
			}
			else{
				// no intersections between rayij and the scene
				img[((H-i-1) * W + j)*3] = 0;  
				img[((H-i-1) * W + j)*3+1] = 0;
				img[((H-i-1) * W + j)*3+2] = 0;
			}
		}
	}
	return img;
}*/

std::vector<unsigned char> create_scene_shadow_specular_img(int W, int H, double fov, Vector C, Scene Sc, Vector L, double I, double eps, double Nair, double Kmax){
	/* Creates an image, with shadows, from the point of view of a camera located at C (with a field of view fov in radian)  
		of a Scene Sc (consisting of Spheres) illuminated by the light source located at L with an intensity I.
		To calculate shadows, we will shift by eps the intersection point between the ray comming from the camera and a sphere.*/
	std::vector<unsigned char> img(W*H * 3, 0);
	#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < H; i++){
		cout << i << endl;
		for(int j = 0; j < W; j++){
			Ray rayij(C, i, j, H, W, fov);
			bool inters = Sc.inters_scene(rayij);
			Material color;
			if(inters){
				//calculating the intersection point between rayij and the scene
				std::vector<Vector> res = Sc.inters_scene_point(rayij); // intersection point between rayij and Sc, the corresponding normal and "color"
				Vector P = res[0]; // intersection point
				Vector n = res[1]; // normal
				Vector Peps = P + n*eps; // point from which we will compute the shadows
				Vector D = L+(-Peps); // vector going from Peps to the light source L
				double d = D.norm(); // distance between Peps and L
				Vector l = D.normalize(); // unitary vector going from Peps to L

				//creating shadows
				Ray rayP(Peps, l); // creation of a ray of origin the intersection point (shifted by epsilon) Peps and direction l
				bool intersbis = Sc.inters_scene(rayP); // intersection between rayP and the scene?
				if(intersbis){
					// rayP intersects the scene
					std::vector<Vector> resbis = Sc.inters_scene_point(rayP);
					Vector Pbis = resbis[0]; // intersection point between rayP and the scene
					//Vector nbis = resbis[1]; // normal
					Vector PPbis = Pbis-Peps; // vector going from Peps to Pbis
					double dbis = PPbis.norm(); // distance between Peps and Pbis
					//cout << "d " << d << "dbis " << dbis << endl;
					//if(c.x()==0 and c.y()==1 and c.z()==0){cout<< "green" << endl;}
					//cout << resbis[2].x() << " " << resbis[2].y() << " " << resbis[2].z() << endl;

					if(dbis < d ){
						// there is a shadow
						//cout << c.x() << " " << c.y() << " " << c.z() << endl;
						//cout << "shadow" << endl;
						img[((H-i-1) * W + j)*3] = 0;  
						img[((H-i-1) * W + j)*3+1] = 0;
						img[((H-i-1) * W + j)*3+2] = 0;
					}
					if(dbis >= d){
						// there is no shadow
						//cout << "no shadow" << endl;
						//if(c.x()==0 and c.y()==1 and c.z()==0){cout<< "green" << endl;}
						double intensity = pow(std::max(0., std::min(l.scalar_prod(n)*I/pow( d, 2), 255.) ), 1.); // intensity of the pixel i,j
						
						for(int k=0; k<Kmax; k++){
			
							Material kcolor = Sc.get_color_brdf(rayij, 5, L, I, eps, Nair);
							color.set_c(color.c()+kcolor.c());
						}
						color.set_c(color.c()*(1./Kmax));
Material Ecolor = Sc.get_color(rayij, 5, L, I, eps, Nair);
			color.set_c(color.c()+Ecolor.c());			// ATTENTION INTENSITE et OMBRES

						//cout << color.c().x() << " " << color.c().y() << " " << color.c().z() << endl;
						img[((H-i-1) * W + j)*3] =  intensity*color.c().x();
						img[((H-i-1) * W + j)*3+1] = intensity*color.c().y();
						img[((H-i-1) * W + j)*3+2] = intensity*color.c().z();
					}
				}
				else{
					// rayP does not intersect the scene, so there is no shadow
					//cout << "no shadow" << endl;
					//if(c.x()==0 and c.y()==1 and c.z()==0){cout<< "green" << endl;}
					double intensity = pow(std::max(0., std::min(l.scalar_prod(n)*I/pow( d, 2), 255.) ), 1.); // intensity of the pixel i,j
					for(int k=0; k<Kmax; k++){
			
						Material kcolor = Sc.get_color_brdf(rayij, 5, L, I, eps, Nair);
						color.set_c(color.c()+kcolor.c());
					}
					color.set_c(color.c()*(1./Kmax));
Material Ecolor = Sc.get_color(rayij, 5, L, I, eps, Nair);
			color.set_c(color.c()+Ecolor.c());			

					//cout << color.c().x() << " " << color.c().y() << " " << color.c().z() << endl;					
					img[((H-i-1) * W + j)*3] =  intensity*color.c().x();
					img[((H-i-1) * W + j)*3+1] = intensity*color.c().y();
					img[((H-i-1) * W + j)*3+2] = intensity*color.c().z();

				}
			}
			else{
				// no intersections between rayij and the scene
				img[((H-i-1) * W + j)*3] = 0;  
				img[((H-i-1) * W + j)*3+1] = 0;
				img[((H-i-1) * W + j)*3+2] = 0;
			}
		}
	}
	return img;
}

//std::vector<unsigned char> create_scene_shadow_specular_img(int W, int H, double fov, Vector C, Scene Sc, Vector L, double I, double eps, double Nair, int Kmax){
	/* Creates an image, with shadows, from the point of view of a camera located at C (with a field of view fov in radian)  
		of a Scene Sc (consisting of Spheres) illuminated by the light source located at L with an intensity I.
		To calculate shadows, we will shift by eps the intersection point between the ray comming from the camera and a sphere.
		The surfaces might be specular.*/
/*	std::vector<unsigned char> img(W*H * 3, 0);
	//#pragma omp parallel for schedule(dynamic,1)
	for(int i = 0; i < H; i++){
		cout << i << endl;
		for(int j = 0; j < W; j++){
			Ray rayij(C, i, j, H, W, fov);
			Material color;
			for(int k=0; k<Kmax; k++){
				//cout << "k " << k << endl;
				Material kcolor = Sc.get_color_brdf(rayij, 5, L, I, eps, Nair);
				color.set_c(color.c()+kcolor.c());
			//cout << (kcolor.c()).x() << " " << (kcolor.c()).y() << " " << (kcolor.c()).z() << endl;
			
			}
			color.set_c(color.c()*(1./Kmax));
//cout << (color.c()).x() << " " << (color.c()).y() << " " << (color.c()).z() << endl;
			//Material Ecolor = Sc.get_color(rayij, 5, L, I, eps, Nair);
			//color.set_c(color.c()+Ecolor.c());			

			img[((H-i-1) * W + j)*3] = pow(color.c().x(), (1.));  
			img[((H-i-1) * W + j)*3+1] = pow(color.c().y(), (1.));
			img[((H-i-1) * W + j)*3+2] = pow(color.c().z(), (1.));
		}
	}
	return img;
}*/

// ============================================================================
//                                   Main 
// ============================================================================
int main(int argc, char* argv[]){

	// parallelization
	int core= 4; // number of cores
	omp_set_num_threads(core);


	//// First algorithm ////
	// create the center of the camera
	Vector C(0,0,35);
	// create the "color" vectors
	Vector vectm1(1,1,1);
	Material m1(vectm1); // color white, diffuse
	// create a sphere
	Vector vect1(0,0,0);
	Sphere sphere1(vect1, 10, m1);
	cout << sphere1.v().x() << ", " << sphere1.r() << endl;
	//create an image of a white circle and a black background
	int W = 500;
	int H = 500;
	//std::v4ctor<unsigned char> img = create_bin_img(W, H, 1.0472, C, sphere1); // fov = 60 degree = 1.0472 rad
	//save_image("WhiteDisk.bmp", &img[0], W, H);	

	//// Second algorithm ////
	// create a light source
	Vector L(-10, 20, 40);
	double I = 300000;
	//create an image of a white sphere and a black background
	//std::vector<unsigned char> img2 = create_vol_img(W, H, 1.0472, C, sphere1, L, I); // fov = 60 degree = 1.0472 rad
	//save_image("WhiteSphere_I300000.bmp", &img2[0], W, H);	

	/// Creation of a colored scene ///
	// create the "color" vectors
	Vector vectm2(0,0,1); // blue
	Material m2(vectm2); 
	Vector vectm3(0,1,0); // green
	Material m3(vectm3); 
	Vector vectm4(1,0,0); // red
	Material m4(vectm4); 
	Vector vectm5(1,1,1); // white
	Material m5(vectm5); 
	Vector vectm6(1,0,1); // magenta
	Material m6(vectm6);
	Vector vectm7(0,1,1); // cyan
	Material m7(vectm7);
	// transparent glass
	double Nsphere = 1.3;
	Material m1t(vectm1, 2, Nsphere); // white diffuse glass
	// create the walls (big spheres)
	Vector vect2(0,1000,0); // sphere2 (ceiling and red)
	Sphere sphere2(vect2, 940, m4);
	Vector vect3(0,0,-1000); // sphere3 (behind the main sphere sphere1 and green)
	Sphere sphere3(vect3, 940, m3);
	Vector vect4(0,-1000, 0); // sphere4 (floor and blue)	
	Sphere sphere4(vect4, 990, m2);
	Vector vect5(0,0,1000); // sphere5 (behind the camera and white)
	Sphere sphere5(vect5, 940, m1);
	Vector vect6(1000, 0,0); // sphere6 (right-hand side and magenta)
	Sphere sphere6(vect6, 940, m6);
	Vector vect7(-1000, 0,0); // sphere7 (left-hand side and cyan)
	Sphere sphere7(vect7, 940, m7);
	// create an additional white specular sphere
	Vector vect8(30,0,0);
	Sphere sphere8(vect8, 5, m1);
	// create an additional white transparent glass sphere
	Vector vect9(-30,0,0);
	Sphere sphere9(vect9, 10, m1t);
	// create a scene
	std::vector<Sphere> sc1;
	sc1.push_back(sphere1);
	Scene scene1(sc1 , 1);
	scene1.add_sphere(sphere2);
	scene1.add_sphere(sphere3);
	scene1.add_sphere(sphere4);
	scene1.add_sphere(sphere5);
	scene1.add_sphere(sphere6);
	scene1.add_sphere(sphere7);
	//scene1.add_sphere(sphere8);
	//scene1.add_sphere(sphere9);
	// create an image of the scene	
	//std::vector<unsigned char> img3 = create_scene_img(W, H, 1.0472, C, scene1, L, I); // fov = 60° = 1.0472 rad
	//save_image("TestColorScene_I300000.bmp", &img3[0], W, H);

	/// Creation of a colored scene with shadows ///
	double eps = .01;
	/*std::vector<unsigned char> img4 = create_scene_shadow_img(W, H, 1.22173, C, scene1, L, I, eps); // fov = 60° = 1.0472 rad
	save_image("ShadowColorScene_I300000_eps0pt01.bmp", &img4[0], W, H);*/

	/// Creation of a colored scene with shadows and specular surfaces ///
	double Nair = 1.;
	//std::vector<unsigned char> img5 = create_scene_shadow_specular_img(W, H, 1.22173, C, scene1, L, I, eps, Nair); // fov = 60° = 1.0472 rad
	//save_image("SpecularScene_2spheres_I1000000_eps0pt01.bmp", &img5[0], W, H);
	
	/// Creation of a colored scene with shadows and transparent surfaces ///
	//std::vector<unsigned char> img5 = create_scene_shadow_specular_img(W, H, 1.22173, C, scene1, L, I, eps, Nair); // fov = 60° = 1.0472 rad
	//save_image("TransparentScene_I1000000_eps0pt01_Nsphere1pt3.bmp", &img5[0], W, H);

	/// Creation of a colored scene with shadows and diffuse surfaces with indirect lighting ///
	int Kmax = 100;
	std::vector<unsigned char> img5 = create_scene_shadow_specular_img(W, H, 1.22173, C, scene1, L, I, eps, Nair, Kmax); // fov = 60° = 1.0472 rad
	save_image("Test_IndirectScene_I300000_eps0pt01_Nsphere1pt3_Kmax100.bmp", &img5[0], W, H);

	
	return 0;
}
