/*
 * Vector.h
 *
 *  Created on: Mar 16, 2015
 *      Author: erick
 */

#ifndef VECTOR_HPP_
#define VECTOR_HPP_
#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

namespace MTX {

template<class T> class Matrix;
template<class T> class Vector;

template<class T>
std::ostream& operator<<(std::ostream& os, const Vector<T>& v);
template<class T> Vector<T> operator*(const Matrix<T>& m, const Vector<T>& v);
template<class T> Vector<T> operator*(const T& s, const Vector<T>& v);
template<class T> Vector<T> project(const T& s, const Vector<T>& v);

template<class T>
class Vector {
protected:
	unsigned int n_elements;
	std::vector<T> storage;
public:
	friend Vector operator*<>(const Matrix<T>&, const Vector&);
	friend Vector operator*<>(const T& s, const Vector&);
	friend Vector project<>(const T& s, const Vector<T>& v);
	Vector(int size);
	Vector(const Vector& v);
	Vector(const std::vector<T>& v);
	Vector<T>& operator=(const Vector<T>& v);
	T& operator[](const unsigned int& idx);
	T operator[](const unsigned int& idx) const;
	virtual ~Vector();
	Vector<T> operator+(const Vector<T>& v);
	Vector<T> operator-(const Vector<T>& v);
	T operator *(const Vector<T>& v);
	Vector<T> operator*(const T& rhs);
	Vector<T> operator/(const T& rhs);
	Vector<T>& operator +=(const Vector<T>& v);
	Vector<T>& operator -=(const Vector<T>& v);
	Vector<T>& operator *=(const T& rhs);
	Vector<T>& operator /=(const T& rhs);
	unsigned int getSize() const;
	void zeros();
	void ones();
	void fill(const T& val);
	void resize(const int& n);
	Vector<T>& modulo(const T& length);
	T length();
	T angle(const Vector<T>& v);
	Vector<T> unit();
	friend std::ostream& operator<<<>(std::ostream&, const Vector<T>& v);
};

template<class T>
inline MTX::Vector<T>::Vector(int size) :
		n_elements(size), storage(size) {
}

template<class T>
inline MTX::Vector<T>::~Vector() {
	std::vector<T>().swap(storage);
}

template<class T>
inline unsigned int MTX::Vector<T>::getSize() const {
	return n_elements;
}

template<class T>
inline Vector<T> MTX::Vector<T>::operator +(const Vector<T>& v) {
	if (v.n_elements != n_elements) {
		throw std::logic_error("Invalid vector operation: sizes don't match");
	}
	Vector<T> out = Vector(n_elements);
	for (unsigned int i = 0; i < n_elements; ++i) {
		out[i] = storage[i] + v[i];
	}
	return out;
}

template<class T>
inline Vector<T> MTX::Vector<T>::operator -(const Vector<T>& v) {
	if (v.n_elements != n_elements) {
		throw std::logic_error("Invalid vector operation: sizes don't match");
	}
	Vector<T> out = Vector(n_elements);
	for (unsigned int i = 0; i < n_elements; ++i) {
		out[i] = storage[i] - v[i];
	}
	return out;
}

template<class T>
inline T MTX::Vector<T>::operator *(const Vector<T>& v) {
	if (v.n_elements != n_elements) {
		throw std::logic_error("Invalid vector operation: sizes don't match");
	}
	T dot = 0;
	for (unsigned int i = 0; i < n_elements; ++i) {
		dot += storage[i] * v[i];
	}
	return dot;
}

template<class T>
inline Vector<T>& MTX::Vector<T>::operator +=(const Vector<T>& v) {
	if (v.n_elements != n_elements) {
		throw std::logic_error("Invalid vector operation: sizes don't match");
	}
	for (unsigned int i = 0; i < n_elements; ++i) {
		storage[i] += v[i];
	}
	return *this;
}

template<class T>
inline Vector<T>& MTX::Vector<T>::operator -=(const Vector<T>& v) {
	for (unsigned int i = 0; i < n_elements; ++i) {
		storage[i] -= v[i];
	}
	return *this;
}

template<class T>
inline void MTX::Vector<T>::zeros() {
	for (unsigned int i = 0; i < storage.size(); ++i) {
		storage[i] = 0.0;
	}
}

template<class T>
inline void MTX::Vector<T>::ones() {
	for (unsigned int i = 0; i < n_elements; ++i) {
		storage[i] = 1.0;
	}
}

template<class T>
inline void MTX::Vector<T>::fill(const T& val) {
	for (unsigned int i = 0; i < n_elements; ++i) {
		storage[i] = val;
	}
}

template<class T>
inline T MTX::Vector<T>::length() {
//	T out = 0;
//	for (unsigned int i=0;i<n_elements;++i) {
//		out += storage.at(i)*storage.at(i);
//	}
	return sqrt((*this) * (*this));
}

template<class T>
inline Vector<T>& Vector<T>::operator *=(const T& rhs) {
	for (unsigned int i = 0; i < n_elements; ++i) {
		storage[i] *= rhs;
	}
	return *this;
}

template<class T>
inline Vector<T>& Vector<T>::operator /=(const T& rhs) {
	for (unsigned int i = 0; i < n_elements; ++i) {
		storage[i] /= rhs;
	}
	return *this;
}

template<class T>
inline Vector<T>& Vector<T>::modulo(const T& length) {
	//T lengthi = 1.0 / length;
	for (unsigned i = 0; i < storage.size(); ++i) {
		storage[i] -= length * rint(storage[i] / length);
	}
	return *this;
}

template<class T>
inline T Vector<T>::angle(const Vector<T>& v) {
	return acos((v * (*this)) / sqrt((v * v) * ((*this) * (*this))));
}

template<class T>
inline void Vector<T>::resize(const int& n) {
	storage.resize(n);
	storage.shrink_to_fit();
}

template<class T>
inline Vector<T> Vector<T>::operator *(const T& rhs) {
	Vector<T> out(n_elements);
	for (unsigned int i = 0; i < n_elements; ++i) {
		out[i] = storage[i] * rhs;
	}
	return out;
}

template<class T>
inline Vector<T> Vector<T>::operator /(const T& rhs) {
	T inv = 1.0 / rhs;
	return (*this) * inv;
}

template<class T>
inline Vector<T> MTX::Vector<T>::unit() {
	return (*this) / length();
}

template<class T>
inline MTX::Vector<T>::Vector(const Vector& v) {
	*this = v;
}

template<class T>
inline MTX::Vector<T>::Vector(const std::vector<T>& v) {
	n_elements = v.size();
	storage.resize(n_elements);
	storage.shrink_to_fit();
	for (unsigned int i = 0; i < n_elements; ++i) {
		storage[i] = v[i];
	}

}

template<class T>
inline Vector<T>& MTX::Vector<T>::operator =(const Vector<T>& v) {
	n_elements = v.n_elements;
	storage.resize(n_elements);
	storage.shrink_to_fit();
	for (unsigned int i = 0; i < n_elements; ++i) {
		storage[i] = v[i];
	}
	return *this;
}

template<class T>
inline T& MTX::Vector<T>::operator [](const unsigned int& idx) {
	if (idx > n_elements) {
		throw std::out_of_range("Index out of range");
	}
	return storage[idx];
}

template<class T>
inline T MTX::Vector<T>::operator [](const unsigned int& idx) const {
	if (idx > n_elements) {
		throw std::out_of_range("Index out of range");
	}
	return storage[idx];
}

} /* namespace MTX */

template<class T>
inline MTX::Vector<T> MTX::operator *(const T& s, const MTX::Vector<T>& v) {
	MTX::Vector<T> out = v;
	return out * s;
}

template<class T>
std::ostream& MTX::operator <<(std::ostream& os, const MTX::Vector<T>& v) {
	//os.precision(10);
	os.setf(std::ios::fixed, std::ios::floatfield);
	os << "(";
	for (unsigned int i = 0; i < v.n_elements; ++i) {
		os << v[i];
		if (i < v.n_elements - 1) {
			os << ", ";
		}
	}
	os << ")";
	return os;
}

template<class T>
inline MTX::Vector<T> MTX::project(const T& s, const MTX::Vector<T>& v) {
	MTX::Vector<T> out = v;
	return s * out.unit();
}

#endif /* VECTOR_HPP_ */
