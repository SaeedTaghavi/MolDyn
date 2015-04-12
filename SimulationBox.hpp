/*
 * SimulationBox.hpp
 *
 *  Created on: Apr 2, 2015
 *      Author: erick
 */

#ifndef SIMULATIONBOX_HPP_
#define SIMULATIONBOX_HPP_
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <bits/random.h>
#include <ctime>
#include <stdexcept>
#include "Vector3d.hpp"
#include "LJPotential.hpp"
#include "VelocityVerlet.hpp"
#include "Particle.hpp"
#include "constants.hpp"

namespace MolDyn {

template<class T>
class SimulationBox {
protected:
	unsigned int n_atoms;
	std::vector<Vector3d<T> > positions;
	std::vector<Vector3d<T> > velocities;
	std::vector<Vector3d<T> > forces;
	Integrator<T>& integrator;
	T ePot;
	T eKin;
	T eTot;
	T pressure;
	T density;
	T length;
	T halfLength;
	T sigmaCut;
	T ELRC;
	T PLRC;
	unsigned int seed;
	//std::default_random_engine generator;
	std::mt19937 generator;
	std::uniform_real_distribution<T> distribution;
public:
	SimulationBox(const unsigned int& n, const T& d, Integrator<T>& in,
			const T& temp, const T& l, const T& scut);
	void initialize();
	void step(bool rescale);
	virtual ~SimulationBox();

	T getEKin() const {
		return eKin;
	}

	T getEPot() const {
		return ePot;
	}

	T getETot() const {
		return eTot;
	}

	const std::vector<Vector3d<T> >& getForces() const {
		return forces;
	}

	const std::vector<Vector3d<T> >& getPositions() const {
		return positions;
	}

	T getPressure() const {
		return pressure;
	}

	const std::vector<Vector3d<T> >& getVelocities() const {
		return velocities;
	}

protected:
	Vector3d<T> randomVector();
	void addPositions();
	void addAtom(const unsigned int& i);
	bool overlaps(const Vector3d<T>& v, const unsigned int& idx);
};

template<class T>
inline MolDyn::SimulationBox<T>::SimulationBox(const unsigned int& n,
		const T& d, Integrator<T>& in, const T& temp, const T& l,
		const T& scut = 0.9) :
		n_atoms(n), density(d), integrator(in), sigmaCut(scut), length(l), positions(
				n, Vector3d<T>()), velocities(n, Vector3d<T>()), forces(n,
				Vector3d<T>()) {
	if (sigmaCut <= 0.0 || sigmaCut >= 1.0) {
		throw std::out_of_range(
				"The overlap cutoff distance is out of range (0:1).");
	}
	halfLength = 0.5 * length;
	pressure = 0;

	std::cout << "Simulation box length = " << length << std::endl;

	ELRC = integrator.potential.getELRC(n_atoms, density);
	PLRC = integrator.potential.getPLRC(density);

	std::cout << "ELRC = " << ELRC << std::endl;
	std::cout << "PLRC = " << PLRC << std::endl;

	sigmaCut *= integrator.potential.getSigma();
	std::cout << "Sigma cut = " << sigmaCut << std::endl;

	seed = time(NULL);
	//generator = std::default_random_engine(seed);
	generator = std::mt19937(seed);
	distribution = std::uniform_real_distribution<T>(0.0, 1.0);
	integrator.setTemp(temp);

}

template<class T>
inline void SimulationBox<T>::initialize() {
	if (!integrator.potential.validateCutoff(halfLength)) {
		std::cout << "******** WARNING: LJCUT > L/2 ********" << std::endl;
		std::cout << "LJCUT = " << integrator.potential.getCutOff()
				<< ",   L/2  = " << halfLength << std::endl;
		std::cout << "******** STOPPING  SIMULATION ********" << std::endl;
		throw std::logic_error("Invalid LJCUT");
	}
	std::cout << "STARTING INITIALIZATION" << std::endl;
	// ASSIGN POSITION OF FIRST ATOM.
	Vector3d<T> p;
	p = randomVector();
	positions[0] = p;
	std::cout << "INSTERTION OF MOLECULE 1: " << positions[0] << std::endl;
	addPositions();
}

template<class T>
inline MolDyn::SimulationBox<T>::~SimulationBox() {
	// Free memory
	std::vector<Vector3d<T> >().swap(positions);
	std::vector<Vector3d<T> >().swap(velocities);
	std::vector<Vector3d<T> >().swap(forces);

}

template<class T>
inline Vector3d<T> MolDyn::SimulationBox<T>::randomVector() {
	Vector3d<T> rdm;
	for (unsigned int i = 0; i < 3; ++i) {
		rdm[i] = distribution(generator) * length;
		rdm[i] -= halfLength;
	}
	return rdm;
}

template<class T>
inline void MolDyn::SimulationBox<T>::step(bool rescale = false) {
	// UPDATE POSITIONS AND VELOCITIES OF THE ATOMS
	integrator.advance(positions, velocities, forces, eKin, ePot, eTot,
			pressure, length, rescale);
	// ADD LONG-RANGE CORRECTIONS TO ENERGY AND PRESSURE
	ePot += ELRC;
	eTot += ELRC;
	pressure += PLRC;
}

template<class T>
inline void MolDyn::SimulationBox<T>::addPositions() {
	Vector3d<T> p, dr, sum, ones;
	ones.ones();
	sum.zeros();
	for (unsigned int i = 1; i < n_atoms; ++i) {
		addAtom(i);
		std::cout.setf(std::ios::fixed, std::ios::floatfield);
		std::cout.precision(5);
		std::cout << "INSTERTION OF MOLECULE " << i + 1 << ": " << positions[i];
		std::cout << std::endl;
	}

	// CALCULATE FORCES
	integrator.forces(positions, forces, ePot, pressure, length);
	// RESCALE VELOCITIES TO (-1,1) & CALCULATE SUM OF X,Y,Z VELOCITY
	for (unsigned int i = 0; i < n_atoms; ++i) {
		p = randomVector();
		velocities[i] = p * 2 - ones;
		sum += velocities[i];
	}
	// ZERO THE TOTAL LINEAR MOMENTUM BY REMOVING/ADDING MOMENTUM
	sum /= n_atoms;
	for (unsigned int i = 0; i < n_atoms; ++i) {
		velocities[i] -= sum;
	}
	// SCALE VELOCITIES TO SET-POINT TEMPERATURE
	integrator.rescale(velocities);
}

template<class T>
inline void MolDyn::SimulationBox<T>::addAtom(const unsigned int& ith) {
	Vector3d<T> p;
	bool retry = true;
	while (retry == true) {
		p = randomVector();
		retry = overlaps(p, ith);
	}
	positions[ith] = p;
}

template<class T>
inline bool MolDyn::SimulationBox<T>::overlaps(const Vector3d<T>& v,
		const unsigned int& ith) {
	Vector3d<T> dr;
	bool result = false;
	T sigmaCut2 = sigmaCut * sigmaCut;
	// CALCULATE DISTANCES FROM PREVIOUSLY INSTERTED ATOMS
	for (unsigned int i = 0; i < ith; ++i) {
		dr = v;
		dr -= positions[i];
		dr.modulo(length);
		if ((dr * dr) < sigmaCut2) {
			result = true;
			break;
		}
	}
	return result;
}

template<typename T>
inline T atomicDensity(const T& dens, const T& mass) {
	return dens * N0 * CM2A / mass;
}

template<typename T>
inline T getSimBoxLength(const unsigned int& atoms, const T& dens) {
	return pow(atoms / dens, 1.0 / 3.0);
}

} /* namespace MolDyn */

#endif /* SIMULATIONBOX_HPP_ */
