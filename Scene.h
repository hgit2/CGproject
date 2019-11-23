#ifndef SCENE_H_  
#define SCENE_H_
#include "Sphere.h"
#include "Ray.h"
#include "Vector.h"
#include <vector>

class Scene{
	public:
	// =====================================================================
	//                               Constructors
	// =====================================================================
	Scene(std::vector<Sphere> Sc, int N); // creates a scene Sc consisting in a vector containing N Spheres

	// =========================================================================
	//                                Destructor
	// =========================================================================
	virtual ~Scene();

	// =========================================================================
	//                                	Getters
	// =========================================================================
	inline  std::vector<Sphere> sc() const;
	inline  int N() const;

	// ===========================================================================
	//                           Public Function members
	// ===========================================================================
	void add_sphere(Sphere S); // add a sphere to the scene
	bool inters_scene(Ray R); // returns whether an intersection point between the ray R and the scene exists
	std::vector<Vector> inters_scene_point(Ray R); // calculates,  if it exists, the first intersection point P between the ray R and the scene as well as the normal to the sphere at P
	Material get_color(Ray R, int nbrebound, Vector L, double I, double eps, double Nair); // get the color of the surface (that might be specular or transparent) intersected by the ray R, with the refractive index of the air Nair
	Material get_color_brdf(Ray R, int nbrebound, Vector L, double I, double eps, double Nair); // get the color of the surface (that might be diffuse, specular or transparent) intersected by the ray R, with the refractive index of the air Nair, using the render equation


	protected:
	// =====================================================================
	//                              Data members 
	// =====================================================================
	std::vector<Sphere> sc_; // the scene (i.e. a vector of spheres)
	int N_; // number of spheres in the scene

};

// ===========================================================================
//                            Getters definitions
// ===========================================================================
inline  std::vector<Sphere> Scene::sc() const{
	return sc_;
}

inline  int Scene::N() const{
	return N_;
}



#endif //SCENE_H_

