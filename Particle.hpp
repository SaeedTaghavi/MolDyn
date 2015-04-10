/*
 * Particle.hpp
 *
 *  Created on: Apr 9, 2015
 *      Author: erick
 */

#ifndef PARTICLE_HPP_
#define PARTICLE_HPP_
#include "Vector.hpp"

namespace MolDyn {

template<class T>
class Particle {
private:
	unsigned char species;
public:
	MTX::Vector<T> position;
	MTX::Vector<T> velocity;
	MTX::Vector<T> force;
public:
	Particle(unsigned int& z);
	Particle(const MolDyn::Particle<T>& rhs);
	void operator=(const MolDyn::Particle<T>& rhs);
	unsigned int getSpecies() const;
	virtual ~Particle();
};

} /* namespace MolDyn */

template<class T>
inline MolDyn::Particle<T>::Particle(unsigned int& z) :
		species(z), position(3), velocity(3), force(3) {
}

template<class T>
inline MolDyn::Particle<T>::Particle(const MolDyn::Particle<T>& rhs) {
	species = rhs.getSpecies();
	position = rhs.position;
	velocity = rhs.velocity;
	force = rhs.force;
}

template<class T>
inline void MolDyn::Particle<T>::operator =(const MolDyn::Particle<T>& rhs) {
	species = rhs.getSpecies();
	position = rhs.position;
	velocity = rhs.velocity;
	force = rhs.force;
}

template<class T>
inline unsigned int MolDyn::Particle<T>::getSpecies() const {
	return species;
}

#endif /* PARTICLE_HPP_ */
