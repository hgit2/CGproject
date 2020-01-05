#ifndef VECTOR_H_  
#define VECTOR_H_
#include <math.h> 

class Vector{

	public:
	// =====================================================================
	//                               Constructors
	// =====================================================================

	// creates a vector (x,y,z)
	Vector(double x, double y, double z); 
	
	// default constructor
	Vector(); 

	// =========================================================================
	//                                Destructor
	// =========================================================================
  	virtual ~Vector();

	// =========================================================================
	//                                  Getters
	// =========================================================================
	inline double x() const;
	inline double y() const;
	inline double z() const;

	// =========================================================================
	//                                 Operators
	// =========================================================================
	
	// multiplication by a float	
	Vector operator*(float const& a);  

	// multiplication by a vector, element-wise
	Vector operator*(Vector const& v); 

	// addition of two vectors (this + v)
	Vector operator+(Vector const& v); 

	// additive inverse of the vector
	Vector operator-(); 

	// substraction of two vectors (this - v)
	Vector operator-(Vector const& v); 

	// checks whether this vector is equal to the vector v
	bool operator==(Vector const& v); 

	// checks whether this vector is not equal to the vector v
	bool operator!=(Vector const& v);

	// checks whether this vector is less than or equal to the vector v
	bool operator<=(Vector const& v); 


	// ===========================================================================
	//                           Public Function members
	// ===========================================================================

	// norm
	double norm(); 

	// squared norm
	double sqnorm(); 

	// scalar product
	double scalar_prod(Vector u); 

	// vector product
	Vector vector_prod(Vector v); 
	
	// normalizes the Vector object
	Vector normalize(); 


	protected:
	// =====================================================================
	//                       Data members (the coordinates)
	// =====================================================================
	double x_;
	double y_;
	double z_;

};

// ===========================================================================
//                            Getters definitions
// ===========================================================================
inline double Vector::x() const {
  return x_;
}

inline double Vector::y() const {
  return y_;
}

inline double Vector::z() const {
  return z_;
}


#endif // VECTOR_H_
