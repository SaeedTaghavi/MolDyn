/*
 * Vector3d.hpp
 *
 *  Created on: Apr 11, 2015
 *      Author: erick
 */

#ifndef VECTOR3D_HPP_
#define VECTOR3D_HPP_
#include <stdexcept>
#include <cmath>
#include <iostream>

namespace MolDyn {

template<class T> class Vector3d;

template<class T>
std::ostream& operator<<(std::ostream& os, const Vector3d<T>& v);
template<class T> Vector3d<T> operator*(const T& s, const Vector3d<T>& v);
template<class T> Vector3d<T> project(const T& s, const Vector3d<T>& v);

template<class T>
class Vector3d {
private:
	T storage[3];
public:
	friend Vector3d operator*<>(const T& s, const Vector3d&);
	friend Vector3d project<>(const T& s, const Vector3d<T>& v);
	Vector3d();
	Vector3d(const T& x, const T& y, const T& z);
	Vector3d(const Vector3d<T>& rhs);
	Vector3d<T>& operator=(const Vector3d<T>& v);
	T& operator[](const unsigned char& idx);
	T operator[](const unsigned char& idx) const;
	Vector3d<T> operator+(const Vector3d<T>& v);
	Vector3d<T> operator-(const Vector3d<T>& v);
	T operator *(const Vector3d<T>& v);
	Vector3d<T> operator*(const T& rhs);
	Vector3d<T> operator/(const T& rhs);
	Vector3d<T>& operator +=(const Vector3d<T>& v);
	Vector3d<T>& operator -=(const Vector3d<T>& v);
	Vector3d<T>& operator *=(const T& rhs);
	Vector3d<T>& operator /=(const T& rhs);
	void zeros();
	void ones();
	void fill(const T& val);
	Vector3d<T>& modulo(const T& length);
	T length();
	T angle(const Vector3d<T>& v);
	Vector3d<T> unit();
	friend std::ostream& operator<<<>(std::ostream&, const Vector3d<T>& v);
	virtual ~Vector3d();
};

template<class T>
inline MolDyn::Vector3d<T>::Vector3d(const T& x, const T& y, const T& z) {
	storage[0] = x;
	storage[1] = y;
	storage[2] = z;
}

template<class T>
inline MolDyn::Vector3d<T>::Vector3d(const Vector3d<T>& rhs) {
	for (unsigned char i = 0; i < 3; ++i) {
		storage[i] = rhs[i];
	}
}

template<class T>
inline Vector3d<T>& MolDyn::Vector3d<T>::operator =(const Vector3d<T>& v) {
	for (unsigned char i = 0; i < 3; ++i) {
		storage[i] = v[i];
	}
	return *this;
}

template<class T>
inline T& MolDyn::Vector3d<T>::operator [](const unsigned char& idx) {
	if (idx >= 3) {
		throw std::out_of_range("The vector index is out of range ");
	}
	return storage[idx];
}

template<class T>
inline T MolDyn::Vector3d<T>::operator [](const unsigned char& idx) const {
	if (idx >= 3) {
		throw std::out_of_range("The vector index is out of range ");
	}
	return storage[idx];
}

template<class T>
inline Vector3d<T> MolDyn::Vector3d<T>::operator +(const Vector3d<T>& v) {
	Vector3d<T> out;
	for (unsigned char i = 0; i < 3; ++i) {
		out[i] = storage[i] + v[i];
	}
	return out;
}

template<class T>
inline Vector3d<T> MolDyn::Vector3d<T>::operator -(const Vector3d<T>& v) {
	Vector3d<T> out;
	for (unsigned char i = 0; i < 3; ++i) {
		out[i] = storage[i] - v[i];
	}
	return out;
}

template<class T>
inline T MolDyn::Vector3d<T>::operator *(const Vector3d<T>& v) {
	T dot = 0;
	for (unsigned char i = 0; i < 3; ++i) {
		dot += storage[i] * v[i];
	}
	return dot;
}

template<class T>
inline Vector3d<T> MolDyn::Vector3d<T>::operator *(const T& rhs) {
	Vector3d<T> out;
	for (unsigned char i = 0; i < 3; ++i) {
		out[i] = storage[i] * rhs;
	}
	return out;
}

template<class T>
inline Vector3d<T> MolDyn::Vector3d<T>::operator /(const T& rhs) {
	T rhsi = 1.0 / rhs;
	return (*this) * rhsi;
}

template<class T>
inline Vector3d<T>& MolDyn::Vector3d<T>::operator +=(const Vector3d<T>& v) {
	for (unsigned char i = 0; i < 3; ++i) {
		storage[i] += v[i];
	}
	return *this;
}

template<class T>
inline Vector3d<T>& MolDyn::Vector3d<T>::operator -=(const Vector3d<T>& v) {
	for (unsigned char i = 0; i < 3; ++i) {
		storage[i] -= v[i];
	}
	return *this;
}

template<class T>
inline Vector3d<T>& MolDyn::Vector3d<T>::operator *=(const T& rhs) {
	for (unsigned char i = 0; i < 3; ++i) {
		storage[i] *= rhs;
	}
	return *this;
}

template<class T>
inline Vector3d<T>& MolDyn::Vector3d<T>::operator /=(const T& rhs) {
	T rhsi = 1.0 / rhs;
	for (unsigned char i = 0; i < 3; ++i) {
		storage[i] *= rhsi;
	}
	return *this;
}

template<class T>
inline void MolDyn::Vector3d<T>::zeros() {
	fill(0.0);
}

template<class T>
inline void MolDyn::Vector3d<T>::ones() {
	fill(1.0);

}

template<class T>
inline void MolDyn::Vector3d<T>::fill(const T& val) {
	for (unsigned char i = 0; i < 3; ++i) {
		storage[i] = val;
	}
}

template<class T>
inline Vector3d<T>& MolDyn::Vector3d<T>::modulo(const T& length) {
	T lengthi = 1.0 / length;
	for (unsigned char i = 0; i < 3; ++i) {
		storage[i] -= length * rint(storage[i] * lengthi);
	}
	return *this;
}

template<class T>
inline T MolDyn::Vector3d<T>::length() {
	T sum = 0;
	for (unsigned char i = 0; i < 3; ++i) {
		sum += storage[i] * storage[i];
	}
	return sqrt(sum);
}

template<class T>
inline T MolDyn::Vector3d<T>::angle(const Vector3d<T>& v) {
	return acos((v * (*this)) / sqrt((v * v) * ((*this) * (*this))));
}

template<class T>
inline Vector3d<T>::Vector3d() {

}

template<class T>
inline Vector3d<T>::~Vector3d() {
}

template<class T>
inline Vector3d<T> MolDyn::Vector3d<T>::unit() {
	return (*this) / length();
}

} /* namespace MolDyn */

template<class T>
inline MolDyn::Vector3d<T> MolDyn::operator *(const T& s,
		const MolDyn::Vector3d<T>& v) {
	MolDyn::Vector3d<T> out = v;
	return out * s;
}

template<class T>
inline MolDyn::Vector3d<T> MolDyn::project(const T& s,
		const MolDyn::Vector3d<T>& v) {
	MolDyn::Vector3d<T> out = v;
	return s * out.unit();
}

template<class T>
inline std::ostream& MolDyn::operator <<(std::ostream& os,
		const Vector3d<T>& v) {
	os.setf(std::ios::fixed, std::ios::floatfield);
	os << "(";
	for (unsigned char i = 0; i < 3; ++i) {
		os << v[i];
		if (i < 2) {
			os << ", ";
		}
	}
	os << ")";
	return os;
}

#endif /* VECTOR3D_HPP_ */
