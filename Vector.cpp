#include <iostream>
#include "Vector.h"
#include <math.h> 

// =====================================================================
//                               Constructors
// =====================================================================
Vector::Vector(double x=0., double y=0., double z=0.): x_(x), y_(y), z_(z){
  }

Vector::Vector(){
	x_=0.;
	y_=0.;
	z_=0.;
}

// =========================================================================
//                                 Operators
// =========================================================================

Vector Vector::operator*(float const& a){
	// multiplication by a float
	Vector result(x_*a, y_*a, z_*a);	
	return result;
}

Vector Vector::operator*(Vector const& v){
	// multiplication by a vector
	Vector result(x_*v.x(), y_*v.y(), z_*v.z());	
	return result;
}

Vector Vector::operator+(Vector const& v){
	// addition of two vectors
	Vector result(x_ + v.x(), y_ + v.y(), z_ + v.z());	
	return result;
}

Vector Vector::operator-(){
	// additive inverse of the vector
	Vector result(-x_, -y_, -z_);	
	return result;
}

Vector Vector::operator-(Vector const& v){
	// substraction of two vectors (this - v)
	Vector result(x_ - v.x(), y_ - v.y(), z_ - v.z());	
	return result;
}

bool Vector::operator==(Vector const& v){
	// checks whether the two vectors are equal
	if(x_ == v.x() and y_ == v.y() and z_ == v.z()){	
		return true;
	}
	else{
		return false;
	}
}

bool Vector::operator!=(Vector const& v){
	// checks whether the two vectors are not equal
	if(x_ != v.x() and y_ != v.y() and z_ != v.z()){	
		return true;
	}
	else{
		return false;
	}
}

bool Vector::operator<=(Vector const& v){
	// checks whether the vector is less than or equal to v
	if(x_ <= v.x() and y_ <= v.y() and z_ <= v.z()){	
		return true;
	}
	else{
		return false;
	}
}

// ===========================================================================
//                           Public Function members
// ===========================================================================

// =============================== Norm =======================================
double Vector::norm(){
	// norm
	double result = sqrt( pow(x_, 2) + pow(y_,2) + pow(z_,2) );
	return result;
}

// ========================== Squared Norm =====================================
double Vector::sqnorm(){
	// squared norm
	double result = pow(x_, 2) + pow(y_,2) + pow(z_,2);
	return result;
}

// ========================== Scalar Product ==================================
double Vector::scalar_prod(Vector u){
	// scalar product
	double result = u.x()*x_ + u.y()*y_ + u.z()*z_;
	return result;
}

// ========================== Vector Product ==================================
Vector Vector::vector_prod(Vector v){
	// vector product between the Vector object and the Vector u
	Vector result( y_*v.z() - z_*v.y(), z_*v.x()-x_*v.z(), x_*v.y()-y_*v.x());
	return result;
}

// ============================ Normalize =====================================
Vector Vector::normalize(){
	// normalizes the Vector object 
	double n = this->norm();
	Vector result(x_/n, y_/n, z_/n);
	return result;
}

// =====================================================================
//                                 Destructor
// =====================================================================
Vector::~Vector() = default;
