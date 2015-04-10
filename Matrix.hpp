/*
 * Matrix.h
 *
 *  Created on: Mar 16, 2015
 *      Author: erick
 */

#ifndef MATRIX_HPP_
#define MATRIX_HPP_
#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

#include "Vector.hpp"

namespace MTX {

template<class T> class Matrix;
template<class T> class Vector;

template<class T> std::ostream& operator<<(std::ostream& os, const Matrix<T>&);
template<class T> Vector<T> operator*(const Matrix<T>& m, const Vector<T>& v);

template<class T>
class Matrix {
public:
	unsigned int n_rows;
	unsigned int n_cols;
private:
	std::vector<std::vector<T> > storage;
public:
	Matrix(unsigned int rows, unsigned int cols);
	Matrix(const Matrix<T>& m);
	friend Vector<T> operator*<>(const Matrix&, const Vector<T>&);
	T& operator()(unsigned int row, unsigned int col);
	T operator()(unsigned int row, unsigned int col) const;
	Matrix<T> operator+(const Matrix<T>& m);
	Matrix<T> operator-(const Matrix<T>& m);
	Matrix<T> operator*(const Matrix<T>& m);
	Matrix<T> operator%(const Matrix<T>& m);
	Matrix<T>& operator=(const Matrix<T>& m);
	Matrix<T>& operator+=(const Matrix<T>& m);
	Matrix<T>& operator-=(const Matrix<T>& m);
	Matrix<T> operator+(const T&);
	Matrix<T> operator-(const T&);
	Matrix<T> operator*(const T&);
	Matrix<T> operator/(const T&);
	Matrix<T>& operator+=(const T&);
	Matrix<T>& operator-=(const T&);
	Matrix<T>& operator*=(const T&);
	Matrix<T>& operator/=(const T&);
	void zeros();
	void ones();
	void eye();
	void fill(const T& val);
	Matrix<T> reduced(unsigned int i, unsigned int j);
	Matrix<T> t();
	T det();
	friend std::ostream& operator<<<>(std::ostream&, const Matrix<T>&);
	virtual ~Matrix();
};

template<class T>
inline MTX::Matrix<T>::Matrix(unsigned int rows, unsigned int cols) :
		n_rows(rows), n_cols(cols) {
	storage.resize(n_rows);
	storage.shrink_to_fit();
	for (unsigned int i = 0; i < n_rows; ++i) {
		storage[i].resize(n_cols);
		storage[i].shrink_to_fit();
	}
}

template<class T>
inline MTX::Matrix<T>::Matrix(const Matrix<T>& m) {
	*this = m;
}

template<class T>
inline MTX::Matrix<T>::~Matrix() {
	for (unsigned int i = 0; i < n_rows; ++i) {
		std::vector<T>().swap(storage[i]);
	}
	std::vector<std::vector<T> >().swap(storage);
}

template<class T>
inline T& MTX::Matrix<T>::operator ()(unsigned int row, unsigned int col) {
	if (row > n_rows) {
		throw std::out_of_range("Row index out of range");
	}
	if (col > n_cols) {
		throw std::out_of_range("Column index out of range");
	}
	return storage[row][col];
}

template<class T>
inline T MTX::Matrix<T>::operator ()(unsigned int row, unsigned int col) const {
	if (row > n_rows) {
		throw std::out_of_range("Row index out of range");
	}
	if (col > n_cols) {
		throw std::out_of_range("Column index out of range");
	}
	return storage[row][col];
}

template<class T>
inline Matrix<T> MTX::Matrix<T>::operator +(const Matrix<T>& m) {
	Matrix<T> out(m.n_rows, m.n_cols);
	out.zeros();
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = 0; j < n_cols; ++j) {
			out(i, j) = storage[i][j] + m(i, j);
		}
	}
	return out;
}

template<class T>
inline Matrix<T> MTX::Matrix<T>::operator -(const Matrix<T>& m) {
	Matrix<T> out(m.n_rows, m.n_cols);
	out.zeros();
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = 0; j < n_cols; ++j) {
			out(i, j) = storage[i][j] - m(i, j);
		}
	}
	return out;
}

template<class T>
inline Matrix<T> MTX::Matrix<T>::operator *(const Matrix<T>& m) {
	if (n_cols != m.n_rows) {
		throw std::out_of_range("The number of cols of matrix 1 "
				"should match the number of rows of matrix 2");
	}
	Matrix<T> out(n_rows, m.n_cols);
	out.zeros();
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = 0; j < m.n_cols; ++j) {

			for (unsigned int k = 0; k < n_cols; k++) {
				out(i, j) += storage.at(i).at(j) * m(k, j);
			}
		}
	}

	return out;
}

template<class T>
inline Matrix<T> MTX::Matrix<T>::operator %(const Matrix<T>& m) {
	if (n_cols != m.n_rows) {
		throw std::out_of_range("The number of cols of matrix 1 "
				"should match the number of rows of matrix 2");
	}

	Matrix<T> out(n_rows, m.n_cols);
	for (unsigned int i = 0; i < n_cols; ++i) {
		for (unsigned int j = 0; j < m.n_rows; ++j) {
			out(i, j) = storage.at(i).at(j) * m(i, j);
		}
	}

	return out;
}

template<class T>
inline Matrix<T>& MTX::Matrix<T>::operator =(const Matrix<T>& m) {
	n_rows = m.n_rows;
	n_cols = m.n_cols;
	storage.resize(n_rows);
	storage.shrink_to_fit();
	for (unsigned int i = 0; i < n_rows; ++i) {
		storage[i].resize(n_cols);
		storage[i].shrink_to_fit();
	}
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = 0; j < n_cols; ++j) {
			storage[i][j] = m(i, j);
		}
	}
	return *this;
}

template<class T>
inline Matrix<T>& MTX::Matrix<T>::operator +=(const Matrix<T>& m) {
	if (m.n_rows != n_rows || m.n_cols != n_cols) {
		throw std::out_of_range(
				"Sum operator takes matrices of the same size.");
	}
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = 0; j < n_cols; ++j) {
			storage.at(i).at(j) += m(i, j);
		}
	}
	return *this;
}

template<class T>
inline Matrix<T>& MTX::Matrix<T>::operator -=(const Matrix<T>& m) {
	if (m.n_rows != n_rows || m.n_cols != n_cols) {
		throw std::out_of_range(
				"Sum operator takes matrices of the same size.");
	}
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = 0; j < n_cols; ++j) {
			storage.at(i).at(j) -= m(i, j);
		}
	}
	return *this;
}

template<class T>
inline Matrix<T> Matrix<T>::operator +(const T& rhs) {
	Matrix<T> out(n_rows, n_cols);
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = 0; j < n_cols; ++j) {
			out(i, j) = rhs + storage[i][j];
		}
	}
	return out;
}

template<class T>
inline Matrix<T> Matrix<T>::operator -(const T& rhs) {
	Matrix<T> out(n_rows, n_cols);
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = 0; j < n_cols; ++j) {
			out(i, j) = storage[i][j] - rhs;
		}
	}
	return out;
}

template<class T>
inline Matrix<T> Matrix<T>::operator *(const T& rhs) {
	Matrix<T> out(n_rows, n_cols);
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = 0; j < n_cols; ++j) {
			out(i, j) = storage[i][j] * rhs;
		}
	}
	return out;
}

template<class T>
inline Matrix<T> Matrix<T>::operator /(const T& rhs) {
	Matrix<T> out(n_rows, n_cols);
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = 0; j < n_cols; ++j) {
			out(i, j) = storage[i][j] / rhs;
		}
	}
	return out;
}

template<class T>
inline Matrix<T>& Matrix<T>::operator +=(const T& rhs) {
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = 0; j < n_cols; ++j) {
			storage[i][j] += rhs;
		}
	}
	return *this;
}

template<class T>
inline Matrix<T>& Matrix<T>::operator -=(const T& rhs) {
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = 0; j < n_cols; ++j) {
			storage[i][j] -= rhs;
		}
	}
	return *this;
}

template<class T>
inline Matrix<T>& Matrix<T>::operator *=(const T& rhs) {
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = 0; j < n_cols; ++j) {
			storage[i][j] *= rhs;
		}
	}
	return *this;
}

template<class T>
inline Matrix<T>& Matrix<T>::operator /=(const T& rhs) {
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = 0; j < n_cols; ++j) {
			storage[i][j] /= rhs;
		}
	}
	return *this;
}

template<class T>
inline void MTX::Matrix<T>::zeros() {
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = 0; j < n_cols; ++j) {
			storage[i][j] = 0;
		}
	}
}

template<class T>
inline void MTX::Matrix<T>::ones() {
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = 0; j < n_cols; ++j) {
			storage[i][j] = 1;
		}
	}
}

template<class T>
inline void MTX::Matrix<T>::eye() {
	if (n_cols != n_rows) {
		throw std::logic_error("Not a square matrix");
	}
	zeros();
	for (unsigned int i = 0; i < n_rows; ++i) {
		storage[i][i] = 1;
	}
}

template<class T>
inline void MTX::Matrix<T>::fill(const T& val) {
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = 0; j < n_cols; ++j) {
			storage[i][j] = val;
		}
	}
}

template<class T>
inline Matrix<T> MTX::Matrix<T>::reduced(unsigned int row, unsigned int col) {
	Matrix<T> red = Matrix(row - 1, col - 1);
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = 0; j < n_cols; ++j) {
			if (i != row && j != col) {
				red(i, j) = storage[i][j];
			}
		}
	}
	return red;
}

template<class T>
inline Matrix<T> MTX::Matrix<T>::t() {
	Matrix<T> trans = Matrix(n_cols, n_rows);
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = 0; j < n_cols; ++j) {
			trans(j, i) = storage[i][j];
		}
	}
	return trans;
}

template<class T>
inline T MTX::Matrix<T>::det() {
	return determinant(*this);
}

template<class T>
inline T determinant(Matrix<T> in) {
	if (in.n_cols != in.n_rows) {
		throw std::logic_error("Matrix must be square.");
	}
	double sign = 1.0;
	T sum = 0.0;
	if (in.n_rows > 2) {
		for (unsigned int i = 0; i < in.n_cols; ++i) {
			sign = (i % 2 == 0) ? 1.0 : -1.0;
			sum = sum + sign * determinant(in.reduced(0, 1));
		}
	} else {
		sum = sum + in(0, 0) * in(1, 1) - in(0, 1) * in(1, 0);
	}
	return sum;

}

template<class T>
inline Vector<T> operator*(const Matrix<T>& m, const Vector<T>& v) {
	if (m.n_cols != v.size()) {
		throw std::out_of_range("The number of columns in the matrix "
				"does not match the elements of the vector");
	}
	Vector<T> out(m.n_rows);
	out.zeros();
	for (unsigned int i = 0; i < m.n_rows; ++i) {
		for (unsigned int j = 0; j < v.size(); ++j) {
			out[i] += m(i, j) * v[j];
		}
	}
	return out;
}

} /* namespace MTX */

template<class T>
std::ostream& MTX::operator <<(std::ostream& os, const MTX::Matrix<T>& m) {
	//os.precision(10);
	os.setf(std::ios::fixed, std::ios::floatfield);
	for (unsigned int i = 0; i < m.n_rows; ++i) {
		os << "[";
		for (unsigned int j = 0; j < m.n_cols; ++j) {
			os << m(i, j);
			if (j < m.n_cols - 1) {
				os << ", ";
			}
		}
		os << "]";
		os << std::endl;
	}
	os << std::endl;
	return os;
}

#endif /* MATRIX_HPP_ */
