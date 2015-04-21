/*
 * Component.hpp
 *
 * Provides a class envelope for a simulation component holding
 * the required parameters
 *
 *  Created on: Apr 20, 2015
 *      Author: Erick R Martinez Loran
 */

#ifndef COMPONENT_HPP_
#define COMPONENT_HPP_

#include <string>
#include "constants.hpp"

namespace MolDyn {

template<class T>
class Component {
public:
	T sigma;
	T epsilon;
	T epsilonKB;
	T mass;
	T densityGCM;
	T density;
	unsigned int nParticles;
	std::string identifier;
public:
	Component(const T& _sigma, const T& _epsilon, const T& _mass,
			const T& _density, const unsigned int& _npart,
			const std::string& _identifier);
	virtual ~Component();
};

template<class T>
inline MolDyn::Component<T>::Component(const T& _sigma, const T& _epsilon,
		const T& _mass, const T& _density, const unsigned int& _npart,
		const std::string& _identifier) :
		sigma(_sigma), epsilonKB(_epsilon), epsilon(_epsilon * KB), mass(_mass), density(
				_density * N0 * CM2A / _mass), densityGCM(_density), nParticles(
				_npart), identifier(_identifier) {
}

template<class T>
inline MolDyn::Component<T>::~Component() {
}

} /* namespace MolDyn */

#endif /* COMPONENT_HPP_ */
