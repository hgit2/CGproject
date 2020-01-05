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
#include <time.h> // for the timer

using std::cout;
using std::endl;

// creates the image of a white disk
std::vector<unsigned char> create_bin_img(int W, int H, double fov, Vector C, Sphere S){
	std::vector<unsigned char> img(W*H * 3, 0);
	for(int i = 0; i < H; i++){
		for(int j = 0; j < W; j++){
			Ray rayij(C, i, j, H, W, fov); // ray going toward pixel (j,i)
			bool inters = rayij.inters_sphere(S);	
			if(inters){
				// pixel at (j,i) is white
				img[(i * W + j)*3] = 255;  
				img[(i * W + j)*3+1] = 255;
				img[(i * W + j)*3+2] = 255;
			}				
			else{
				// pixel at (j,i) is black
				img[(i * W + j)*3] = 0;  
				img[(i * W + j)*3+1] = 0;
				img[(i * W + j)*3+2] = 0;
			}
		}
	}
	return img;
}


/**
 creates an image from the point of view of a camera located at C (with a field of view fov in radian)  
 of the sphere S illuminated by the light source located at L with an intensity I
*/	
std::vector<unsigned char> create_vol_img(int W, int H, double fov, Vector C, Sphere S, Vector L, double I){
	std::vector<unsigned char> img(W*H * 3, 0);
	Vector vect(0,0,10);
	double intensity;

	for(int i = 0; i < H; i++){
		for(int j = 0; j < W; j++){
			Ray rayij(C, i, j, H, W, fov); 
			Vector P = rayij.inters_sphere_point(S); // intersection point between the ray rayij and the sphere S
			if(P==(C+vect)){
				// no intersection
				intensity = 0.;
			}
			else{
				// there is an intersection
				Vector D = L-P; // vector going from P to the light source L
				double d = D.sqnorm(); // distance between P and L squared
				Vector l = D.normalize(); // unitary vector going from P to L				
				Vector N = P-S.v(); // normal vector to the sphere S at the point P				
				Vector n = N.normalize(); // normalizing N
				intensity = std::max(0., std::min(l.scalar_prod(n)*I/d, 1.) ); // intensity of the pixel i,j
			}
			// output pixel
			img[((H-i-1) * W + j)*3] =  std::max(0., std::min(intensity*255., 255.)); // to make sure the pixel intensity is correctly bounded
			img[((H-i-1) * W + j)*3+1] = std::max(0., std::min(intensity*255., 255.));
			img[((H-i-1) * W + j)*3+2] = std::max(0., std::min(intensity*255., 255.));	
		}
	}
	return img;
}

/**
 creates an image from the point of view of a camera located at C (with a field of view fov in radian)  
 of a Scene Sc (consisting of Spheres) illuminated by the light source located at L with an intensity I
*/
std::vector<unsigned char> create_scene_img(int W, int H, double fov, Vector C, Scene Sc, Vector L, double I){
	std::vector<unsigned char> img(W*H * 3, 0);

	Vector vect(0,0,10);
	for(int i = 0; i < H; i++){
		for(int j = 0; j < W; j++){
			Ray rayij(C, i, j, H, W, fov);
			std::vector<Vector> res = Sc.inters_scene_point(rayij); // intersection between the ray rayij and the scene Sc
			Vector P = res[0]; // intersection point
			if(P==(C+vect)){
				// no intersection
				img[((H-i-1) * W + j)*3] =  0.;
				img[((H-i-1) * W + j)*3+1] = 0.;
				img[((H-i-1) * W + j)*3+2] = 0.;	
				
			}
			else{
				// there is an intersection
				Vector n = res[1]; // normal
				Vector c = res[2]; // color of the intersected sphere
				Vector D = L+(-P); // vector going from P to the light source L
				double d = D.norm(); // distance between P and L
				Vector l = D.normalize(); // unitary vector going from P to L
				double intensity = std::max(0., std::min(l.scalar_prod(n)*I/pow( d, 2), 1.) ); // intensity of the pixel i,j
				img[((H-i-1) * W + j)*3] =  std::max(0., std::min(intensity*c.x()*255., 255.));
				img[((H-i-1) * W + j)*3+1] = std::max(0., std::min(intensity*c.y()*255., 255.));
				img[((H-i-1) * W + j)*3+2] = std::max(0., std::min(intensity*c.z()*255., 255.));
			}
		}
	}
	return img;
}


/**
 Creates an image, with shadows, from the point of view of a camera located at C (with a field of view fov in radian)  
 of a Scene Sc (consisting of Spheres) illuminated by the light source located at L with an intensity I.
 To calculate shadows, the intersection point between the ray coming from the camera and a sphere is shifted by eps.
*/
std::vector<unsigned char> create_scene_shadow_img(int W, int H, double fov, Vector C, Scene Sc, Vector L, double I, double eps){
	std::vector<unsigned char> img(W*H * 3, 0);
	Vector vect(0,0,10);
	double intensity;

	for(int i = 0; i < H; i++){
		for(int j = 0; j < W; j++){
			Ray rayij(C, i, j, H, W, fov);
			std::vector<Vector> res = Sc.inters_scene_point(rayij); // intersection between rayij and Sc
			Vector P = res[0]; // intersection point
			if(P==(C+vect)){
				// no intersection between rayij and the scene
				img[((H-i-1) * W + j)*3] =  0.;
				img[((H-i-1) * W + j)*3+1] = 0.;
				img[((H-i-1) * W + j)*3+2] = 0.;	
				
			}

			else{
				// calculating the intersection point between rayij and the scene
				Vector n = res[1]; // normal
				Vector c = res[2]; // color of the sphere that is intersected
				Vector Peps = P + n*eps; // point from which we will compute the shadows
				Vector D = L-Peps; // vector going from Peps to the light source L
				double d = D.norm(); // distance between Peps and L
				Vector l = D.normalize(); // unitary vector going from Peps to L

				//creating shadows
				Ray rayP(Peps, l); // creation of a ray of origin the intersection point (shifted by eps) Peps and direction l
				std::vector<Vector> resbis = Sc.inters_scene_point(rayP);
				Vector Pbis = resbis[0]; // intersection point between rayP and the scene
				if(Pbis==(Peps+vect)){
					// rayP does not intersect the scene, so there is no shadow
					intensity = std::max(0., std::min(l.scalar_prod(n)*I/pow( d, 2), 1.) ); // intensity of the pixel i,j
				}
				else{
					// rayP intersects the scene
					Vector PPbis = Pbis-Peps; // vector going from Peps to Pbis
					double dbis = PPbis.norm(); // distance between Peps and Pbis
					if(dbis < d ){
						// there is a shadow
						intensity = 0;					}
					else{
						// there is no shadow
						intensity = std::max(0., std::min(l.scalar_prod(n)*I/pow( d, 2), 1.) ); // intensity of the pixel i,j
					}
				}
				img[((H-i-1) * W + j)*3] =  std::max(0., std::min(intensity*c.x()*255., 255.));
				img[((H-i-1) * W + j)*3+1] = std::max(0., std::min(intensity*c.y()*255., 255.));
				img[((H-i-1) * W + j)*3+2] = std::max(0., std::min(intensity*c.z()*255., 255.));
			}
		}
	}
	return img;
}

/** Creates an image, with shadows, from the point of view of a camera located at C (with a field of view fov in radian)  
 of a Scene Sc (consisting of Spheres of different types: diffuse, specular or transparent) illuminated by the light source
 located at L with an intensity I with gamma correction. To calculate shadows, the intersection point between the ray coming 
 from the camera and a sphere is shifted by eps. There might be specular or transparent surfaces.
*/
std::vector<unsigned char> create_scene_shadow_types_img(int W, int H, double fov, Vector C, Scene Sc, Vector L, double I, double eps, int nbrebound, double Nair){
	std::vector<unsigned char> img(W*H * 3, 0);
	
	for(int i = 0; i < H; i++){
		for(int j = 0; j < W; j++){
			Ray rayij(C, i, j, H, W, fov);
			Material color = Sc.get_color( rayij, nbrebound, L, I, eps, Nair); // getting the "color"
			img[((H-i-1) * W + j)*3] =  pow(color.intensity().x(), 1./2.2)*255.; // gamma correction
			img[((H-i-1) * W + j)*3+1] = pow(color.intensity().y(), 1./2.2)*255.;
			img[((H-i-1) * W + j)*3+2] = pow(color.intensity().z(), 1./2.2)*255.;
		}
	}
	return img;
}




/** Path tracer: creates an image, with shadows, from the point of view of a camera located at C (with a field of view fov in radian)  
 of a Scene Sc (consisting of Spheres) illuminated by the light source located at L with an intensity I with gamma correction.
 Indirect lighting is taken into account. Monte carlo integration is done with Kmax the number of samples (Kmax>=1). 
 If Kmax=1, only direct lighting is taken into account. 
 To calculate shadows, the intersection point between the ray coming from the camera and a sphere will be shifted by eps.
 The surfaces might be specular or transparent.
*/
std::vector<unsigned char> path_tracer(int W, int H, double fov, Vector C, Scene Sc, Vector L, double I, double eps, double Nair, int Kmax, int nbrebound){
	std::vector<unsigned char> img(W*H * 3, 0); // output image

	// looping through the pixels
	#pragma omp parallel for schedule(dynamic,1) // parallelization
	for(int i = 0; i < H; i++){
		for(int j = 0; j < W; j++){
			Ray rayij(C, i, j, H, W, fov);
			Material color;
			
			// Monte Carlo integration
			for(int k=0; k<Kmax; k++){
				color = color + Sc.get_color_brdf_param(rayij, nbrebound, L, I, eps, Nair, Kmax);
			}
			color.set_intensity(color.intensity()*(1./Kmax));

			// output pixel
			img[((H-i-1) * W + j)*3] =  std::max(0., std::min(pow(color.intensity().x(), 1./2.2)*255., 255.)); // gamma correction
			img[((H-i-1) * W + j)*3+1] = std::max(0., std::min(pow(color.intensity().y(), 1./2.2)*255., 255.));
			img[((H-i-1) * W + j)*3+2] = std::max(0., std::min(pow(color.intensity().z(), 1./2.2)*255., 255.));
		}
	}
	return img;
}

// ============================================================================
//                                   Main 
// ============================================================================
int main(int argc, char* argv[]){

	// parallelization
	int core= 6; // number of cores
	omp_set_num_threads(core);
	

	// create the center of the camera
	Vector C(0,0,55);

	// create the "color" vectors
	Vector vectm1(1,1,1);
	Material m1(vectm1); // color white, diffuse
	// create a sphere
	Vector vect1(0,0,0);
	Sphere sphere1(vect1, 10, m1); 
	// dimensions of the image
	int W = 500;
	int H = 500;

//// 1. create an image of a white circle and a black background ////
	//std::vector<unsigned char> img = create_bin_img(W, H, 1.0472, C, sphere1); // fov = 60 degree = 1.0472 rad
	//save_image("WhiteDisk.bmp", &img[0], W, H);	


	// create a light source
	Vector L(-10, 20, 40);
	double I = 3000; // intensity of the light source

//// 2. create an image of a white sphere and a black background ////
	//std::vector<unsigned char> img2 = create_vol_img(W, H, 1.0472, C, sphere1, L, I); // fov = 60 degree = 1.0472 rad
	//save_image("WhiteSphere_I1200.bmp", &img2[0], W, H);	

/// 3. Creation of a colored scene ///

	// create the "color" vectors
	Vector vectm2(0,0,1); 
	Material m2(vectm2); // blue, diffuse
	Vector vectm3(0,1,0); 
	Material m3(vectm3); // green, diffuse
	Vector vectm4(1,0,0); 
	Material m4(vectm4); // red, diffuse
	Vector vectm5(1,1,1); 
	Material m5(vectm5); // white, diffuse
	Vector vectm6(1,0,1); 
	Material m6(vectm6); // magenta, diffuse
	Vector vectm7(0,1,1); 
	Material m7(vectm7); // cyan, diffuse
	Material m7s(vectm7, 1); // cyan, specular
	Vector vectm8(1,1,0); 
	Material m8(vectm8); // yellow, diffuse
	Vector vectm9(.98,.82,0.54); 
	Material m9(vectm9); // yellow sandy like, diffuse
	Vector vectm10(.96, .39,.96); 
	Material m10(vectm10); // warm red, diffuse
	Vector vectm11(.23,.5,.51); 
	Material m11(vectm11); // blue teal, diffuse
	Vector vectm12(1,.61,.36); 
	Material m12(vectm12); // light orange, diffuse
	Vector vectm13(.96,.39,.29); 
	Material m13(vectm13); // dark orange, diffuse

	// white, specular
	Material m1s(vectm1, 1); 

	// transparent glass
	double Nsphere = 1.3;
	Material m1t(vectm1, 2, Nsphere); // white transparent glass

	// create the walls (big spheres)
	Vector vect2(0,1000,0); // sphere2 (ceiling and yellow) 
	Sphere sphere2(vect2, 940, m8);
	Vector vect3(0,0,-1000); // sphere3 (behind the main sphere sphere1 and red) 
	Sphere sphere3(vect3, 940, m4);
	Vector vect4(0,-1000, 0); // sphere4 (floor and yellow sandy) 
	Sphere sphere4(vect4, 990, m9);
	Vector vect5(0,0,1000); // sphere5 (behind the camera and blue teal)
	Sphere sphere5(vect5, 940, m11);
	Vector vect6(1000, 0,0); // sphere6 (right-hand side and dark orange) 
	Sphere sphere6(vect6, 940, m13);
	Vector vect7(-1000, 0,0); // sphere7 (left-hand side and dark orange) 
	Sphere sphere7(vect7, 940, m13);

	// create an additional transparent glass sphere
	Vector vect8(10,2,30);
	Sphere sphere8(vect8,3, m1t);
	// create an additional white specular sphere 
	Vector vect9(-10,0,30);
	Sphere sphere9(vect9, 8, m1s);
	// create an additional white sphere 
	Vector vect10(17,7,17);
	Sphere sphere10(vect10, 5, m1);

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
	scene1.add_sphere(sphere8);
	scene1.add_sphere(sphere9);
	scene1.add_sphere(sphere10);
	

	// create an image of the colored scene	
	//std::vector<unsigned char> img3 = create_scene_img(W, H, 1.0472, C, scene1, L, I); // fov = 60° = 1.0472 rad
	//save_image("ColorScenespheres_I1200.bmp", &img3[0], W, H);

/// 4. Creation of a colored scene with shadows ///
	double eps = 0.001;
	//double before = (clock())/(double)CLOCKS_PER_SEC;
	//std::vector<unsigned char> img4 = create_scene_shadow_img(W, H, 1.22173, C, scene1, L, I, eps); // fov = 60° = 1.0472 rad
	//save_image("ShadowColorScene_I1200_eps0pt001.bmp", &img4[0], W, H);
	//double after = (clock())/(double)CLOCKS_PER_SEC;
	//double diff = after - before; // execution time
	//cout << "time " << diff << endl;
				
/// 5. Creation of a colored scene with shadows, both specular and transparent surfaces ( and gamma correction) ///
	double Nair = 1.;
	int nbrebound = 5;
	//std::vector<unsigned char> img6 = create_scene_shadow_types_img(W, H, 1.22173, C, scene1, L, I, eps, nbrebound, Nair); // fov = 60° = 1.0472 rad
	//save_image("Scene_colors_I3000_eps0pt001_Nsphere1pt3_5rebounds_gamma.bmp", &img6[0], W, H);

/// 6. Creation of a colored scene with shadows and diffuse surfaces with indirect lighting ///
	int Kmax = 1;
	//double before = (clock())/(double)CLOCKS_PER_SEC;
	//std::vector<unsigned char> img7 = path_tracer(W, H, 1.22173, C, scene1, L, I, eps, Nair, Kmax, nbrebound); // fov = 60° = 1.0472 rad
	//save_image("Path_tracer_I1200_eps0pt001_Nsphere1pt3_Kmax1_5rebounds_newbrdf_gamma_rand_01_pi_c.bmp", &img7[0], W, H);
	//double after = (clock())/(double)CLOCKS_PER_SEC;
	//double diff = after - before; // execution time
	//cout << "time " << diff << endl;
	
	
	return 0;
}
