/*
 * SymmetricMatrix.hpp
 * 
 * A sparse symmetric with basic operations on the vector space
 *
 *  Created on: Apr 17, 2015
 *      Author: Erick Martinez Loran
 */

#ifndef SYMMETRICMATRIX_HPP_
#define SYMMETRICMATRIX_HPP_

#include <iostream>
#include <vector>
#include <stdexcept>

namespace MolDyn {

template<class T> class SymmetricMatrix;
template<class T> std::ostream& operator<<(std::ostream& os,
		const SymmetricMatrix<T>&);

template<class T>
class SymmetricMatrix {
public:
	unsigned int n_rows;
	unsigned int n_cols;
	unsigned int n_elements;
private:
	std::vector<T> storage;
	/*
	 * The array of index pointer to the data which
	 * contains only the addresses to the matrix entries.
	 * Since the matrix is symmetric, only entries from the
	 * upper diagonal are stored.
	 *
	 * Matrix entries can now be accessed using the addresses
	 * pointing to the actual data
	 */
	std::vector<std::vector<T*>> nodes;
public:
	SymmetricMatrix(const unsigned int& size);
	SymmetricMatrix(const SymmetricMatrix<T>& rhs);
	T& operator()(unsigned int row, unsigned int col);
	T operator()(unsigned int row, unsigned int col) const;
	SymmetricMatrix<T> operator+(const SymmetricMatrix<T>& m);
	SymmetricMatrix<T> operator-(const SymmetricMatrix<T>& m);
	SymmetricMatrix<T> operator%(const SymmetricMatrix<T>& m);
	SymmetricMatrix<T>& operator=(const SymmetricMatrix<T>& m);
	SymmetricMatrix<T>& operator+=(const SymmetricMatrix<T>& m);
	SymmetricMatrix<T>& operator-=(const SymmetricMatrix<T>& m);
	SymmetricMatrix<T> operator+(const T&);
	SymmetricMatrix<T> operator-(const T&);
	SymmetricMatrix<T> operator*(const T&);
	SymmetricMatrix<T> operator/(const T&);
	SymmetricMatrix<T>& operator+=(const T&);
	SymmetricMatrix<T>& operator-=(const T&);
	SymmetricMatrix<T>& operator*=(const T&);
	SymmetricMatrix<T>& operator/=(const T&);
	void fill(const T& val);
	friend std::ostream& operator<<<>(std::ostream&, const SymmetricMatrix<T>&);
	virtual ~SymmetricMatrix();
private:
	unsigned int sparseCount(const unsigned int& size) const;
};

template<class T>
inline MolDyn::SymmetricMatrix<T>::SymmetricMatrix(const unsigned int& size = 1) :
		storage(sparseCount(size)), nodes(size, std::vector<T*>(size)), n_elements(
				size * size), n_rows(size), n_cols(size) {
	unsigned int counter = 0;
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = i; j < n_cols; ++j) {
			nodes[i][j] = &storage[counter];
			if (i != j) {
				nodes[j][i] = &storage[counter];
			}
			++counter;
		}
	}
}

template<class T>
inline SymmetricMatrix<T>::SymmetricMatrix(const SymmetricMatrix<T>& rhs) :
		n_elements(rhs.n_elements), n_rows(rhs.n_rows), n_cols(rhs.n_cols) {
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = i; j < n_cols; ++j) {
			nodes[i][j] = rhs(i, j);
		}
	}
}

template<class T>
inline T& SymmetricMatrix<T>::operator ()(unsigned int row, unsigned int col) {
	if (row > n_rows) {
		throw std::out_of_range("Row index out of range");
	}
	if (col > n_cols) {
		throw std::out_of_range("Column index out of range");
	}
	return *nodes[row][col];
}

template<class T>
inline T SymmetricMatrix<T>::operator ()(unsigned int row,
		unsigned int col) const {
	if (row > n_rows) {
		throw std::out_of_range("Row index out of range");
	}
	if (col > n_cols) {
		throw std::out_of_range("Column index out of range");
	}
	return *nodes[row][col];
}

template<class T>
inline SymmetricMatrix<T> SymmetricMatrix<T>::operator +(
		const SymmetricMatrix<T>& m) {
	if (!equalSize((*this), m)) {
		throw std::out_of_range("Trying to sum matrices of different sizes");
	}
	SymmetricMatrix<T> out(*this);
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = i; j < n_cols; ++j) {
			out(i, j) += m(i, j);
		}
	}
	return out;
}

template<class T>
inline SymmetricMatrix<T> SymmetricMatrix<T>::operator -(
		const SymmetricMatrix<T>& m) {
	if (!equalSize((*this), m)) {
		throw std::out_of_range(
				"Trying to subtract matrices of different sizes");
	}
	SymmetricMatrix<T> out(*this);
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = i; j < n_cols; ++j) {
			out(i, j) -= m(i, j);
		}
	}
	return out;
}

template<class T>
inline SymmetricMatrix<T> SymmetricMatrix<T>::operator %(
		const SymmetricMatrix<T>& m) {
	if (!equalSize((*this), m)) {
		throw std::out_of_range(
				"Trying to do Schur product of matrices of different sizes");
	}
	SymmetricMatrix<T> out(*this);
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = i; j < n_cols; ++j) {
			out(i, j) *= m(i, j);
		}
	}
	return out;
}

template<class T>
inline SymmetricMatrix<T>& SymmetricMatrix<T>::operator =(
		const SymmetricMatrix<T>& m) {
	n_rows = m.n_rows;
	n_cols = m.n_cols;
	n_elements = m.n_elements;
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = i; j < n_cols; ++j) {
			nodes[i][j] = m(i, j);
		}
	}
	return *this;
}

template<class T>
inline SymmetricMatrix<T>& SymmetricMatrix<T>::operator +=(
		const SymmetricMatrix<T>& m) {
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = i; j < n_cols; ++j) {
			nodes[i][j] += m(i, j);
		}
	}
	return *this;
}

template<class T>
inline SymmetricMatrix<T>& SymmetricMatrix<T>::operator -=(
		const SymmetricMatrix<T>& m) {
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = i; j < n_cols; ++j) {
			nodes[i][j] -= m(i, j);
		}
	}
	return *this;
}

template<class T>
inline SymmetricMatrix<T> SymmetricMatrix<T>::operator +(const T& rhs) {
	SymmetricMatrix<T> out(*this);
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = i; j < n_cols; ++j) {
			out(i, j) += rhs;
		}
	}
	return out;
}

template<class T>
inline SymmetricMatrix<T> SymmetricMatrix<T>::operator -(const T& rhs) {
	SymmetricMatrix<T> out(*this);
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = i; j < n_cols; ++j) {
			out(i, j) -= rhs;
		}
	}
	return out;
}

template<class T>
inline SymmetricMatrix<T> SymmetricMatrix<T>::operator *(const T& rhs) {
	SymmetricMatrix<T> out(*this);
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = i; j < n_cols; ++j) {
			out(i, j) *= rhs;
		}
	}
	return out;
}

template<class T>
inline SymmetricMatrix<T> SymmetricMatrix<T>::operator /(const T& rhs) {
	T rhsi = 1.0 / rhs;
	return (*this) * rhsi;
}

template<class T>
inline SymmetricMatrix<T>& SymmetricMatrix<T>::operator +=(const T& rhs) {
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = i; j < n_cols; ++j) {
			nodes[i][j] += rhs;
		}
	}
	return *this;
}

template<class T>
inline SymmetricMatrix<T>& SymmetricMatrix<T>::operator -=(const T& rhs) {
	return (*this) += (-1.0 * rhs);
}

template<class T>
inline SymmetricMatrix<T>& SymmetricMatrix<T>::operator *=(const T& rhs) {
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = i; j < n_cols; ++j) {
			nodes[i][j] *= rhs;
		}
	}
	return *this;
}

template<class T>
inline SymmetricMatrix<T>& SymmetricMatrix<T>::operator /=(const T& rhs) {
	return (*this) * (1.0 / rhs);
}

template<class T>
inline void SymmetricMatrix<T>::fill(const T& val) {
	for (unsigned int i = 0; i < n_rows; ++i) {
		for (unsigned int j = i; j < n_cols; ++j) {
			*nodes[i][j] = val;
		}
	}
}

template<class T>
inline SymmetricMatrix<T>::~SymmetricMatrix() {
	std::vector<T>().swap(storage);
}

template<class T>
inline unsigned int MolDyn::SymmetricMatrix<T>::sparseCount(
		const unsigned int& size) const {
	return 0.5 * size * (size + 1);
}

} /* namespace MolDyn */

template<class T>
std::ostream& MolDyn::operator <<(std::ostream& os,
		const MolDyn::SymmetricMatrix<T>& m) {
	os.setf(std::ios::scientific, std::ios::floatfield);
	os << std::endl;
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

template<class T>
bool equalSize(const MolDyn::SymmetricMatrix<T>& m1,
		const MolDyn::SymmetricMatrix<T>& m2) {
	return (m1.n_rows == m2.n_rows && m1.n_cols == m2.n_cols);
}

#endif /* SYMMETRICMATRIX_HPP_ */
