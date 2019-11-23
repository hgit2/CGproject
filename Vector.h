#ifndef VECTOR_H_  
#define VECTOR_H_
#include <math.h> 

class Vector{

	public:
// =====================================================================
//                               Constructors
// =====================================================================
	Vector(double x, double y, double z); // creates a vector (x,y,z)
	Vector(); // default constructor

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
	Vector operator*(float const& a); // multiplication by a float
	Vector operator*(Vector const& v); // multiplication by a vector 
	Vector operator+(Vector const& v); // addition of two vectors
	Vector operator-(); // additive inverse of the vector
	Vector operator-(Vector const& v); // substraction of two vectors (this - v)
	bool operator==(Vector const& v); // checks whether this vector is equal to the vector v
	bool operator!=(Vector const& v); // checks whether this vector is not equal to the vector v
	bool operator<=(Vector const& v); // checks whether this vector is less than or equal to the vector v

// ===========================================================================
//                           Public Function members
// ===========================================================================
	double norm(); // norm

	double sqnorm(); // squared norm

	double scalar_prod(Vector u); //Scalar Product

	Vector vector_prod(Vector v);	// vector product

	Vector normalize();	// normalizes the Vector object

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
