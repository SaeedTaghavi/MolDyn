/*
 * SimulationBox.hpp
 *
 * A Molecular Dynamics simulation box
 *
 *  Created on: Apr 2, 2015
 *      Author: Erick R Martinez Loran
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
#include "constants.hpp"
#include "Component.hpp"

namespace MolDyn {

template<class T>
class SimulationBox {
protected:
	/*
	 * The numner of particles in the simulation box
	 */
	unsigned int nParticles;
	/*
	 * The positions of the particles in the simulation
	 * box
	 */
	std::vector<Vector3d<T> > positions;

	/*
	 * The velocities of the particles in the simulation
	 * box
	 */
	std::vector<Vector3d<T> > velocities;

	/*
	 * The force component in each of the particles
	 */
	std::vector<Vector3d<T> > forces;

	/*
	 * Number of different particles in the simulation box
	 */
	unsigned int maxComponents;

	/*
	 * The array of vector species if none is provided, then all
	 * particles are set to species 0
	 */
	std::vector<unsigned int> species;

	std::vector<T> sigma;
	unsigned int componentCounter;

	/*
	 * The masses of the corresponding species
	 */
	std::vector<T> mass;
	/*
	 * The integrator to be used (e.g. Velocity Verlet)
	 * Used by reference to decrease memory usage
	 */
	Integrator<T>& integrator;
	/*
	 * The potential energy of the simulation box
	 */
	T ePot;
	/*
	 * The kinetic energy of the simulation box
	 */
	T eKin;
	/*
	 * The total energy of the simulation box
	 */
	T eTot;
	/*
	 * The computed pressure
	 */
	T pressure;
	/*
	 * The length of the simulation box
	 */
	T length;
	/*
	 * Half the length of the simulation box
	 * (stored for some calculations)
	 */
	T halfLength;
	/*
	 * The sigma cutoff length. Used to determine the
	 * overlap between partiles.
	 */
	T sigmaCut;
	/*
	 * The long range corrections to the energy
	 */
	T ELRC;
	/*
	 * The long range corrections to the potential
	 */
	T PLRC;
	/*
	 * The seed of the random number generator
	 */
	unsigned int seed;
	//std::default_random_engine generator;
	/*
	 * The Mersenne Twister random number generation
	 */
	std::mt19937 generator;
	/*
	 * The uniform distribution generator
	 * Used to generate a uniform distribution of real
	 * numbers in the range (0.0-1.0]
	 */
	std::uniform_real_distribution<T> distribution;
public:
	SimulationBox(const unsigned int& _nparticles, Integrator<T>& _integrator,
			const T& _temperature, const T& _length, const T& scut,
			const unsigned int& _maxComponents);
	void initialize();
	void step(bool rescale);
	void addComponent(const MolDyn::Component<T>& _component);
	void setSpecies(const unsigned int& idx, const unsigned int& type);
	void setSpecies(const unsigned int& start, const unsigned int& end,
			const unsigned int& type);
	void setMass(const unsigned int& idx, const T& type);
	void setPosition(const unsigned int& idx, const Vector3d<T>& v);
	void setVelocity(const unsigned int& idx, const Vector3d<T>& v);
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

	const std::vector<unsigned int>& getSpecies() const {
		return species;
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
inline MolDyn::SimulationBox<T>::SimulationBox(const unsigned int& _nparticles,
		Integrator<T>& _integrator, const T& _temperature, const T& _length,
		const T& scut, const unsigned int& _maxComponents = 1) :
		nParticles(_nparticles), integrator(_integrator), sigmaCut(scut), length(
				_length), positions(_nparticles, Vector3d<T>()), velocities(
				_nparticles, Vector3d<T>()), forces(_nparticles, Vector3d<T>()), species(
				_nparticles, 0), maxComponents(_maxComponents), sigma(
				_maxComponents), componentCounter(0), pressure(0.0), halfLength(
				0.5 * _length), eTot(0.0), eKin(0.0), ePot(0.0) {
	if (sigmaCut <= 0.0 || sigmaCut >= 1.0) {
		throw std::out_of_range(
				"The overlap cutoff distance is out of range (0:1).");
	}
	std::cout << "Simulation box length = " << length << std::endl;

	// Calculate and store the long range corrections
	ELRC = integrator.potential.getELRC();
	PLRC = integrator.potential.getPLRC();

	std::cout << "ELRC = " << ELRC << std::endl;
	std::cout << "PLRC = " << PLRC << std::endl;

	std::cout << "Sigma cut = " << sigmaCut << "sigma " << std::endl;

	// Get the seed as the epoch time in milliseconds
	seed = time(NULL);
	//generator = std::default_random_engine(seed);

	// Feed the seed to the Mersenne Twister generator
	generator = std::mt19937(seed);

	// Feed the generator to the uniform distribution
	distribution = std::uniform_real_distribution<T>(0.0, 1.0);

	// Store the initial temperature
	integrator.setTemperature(_temperature);

}

template<class T>
inline void SimulationBox<T>::initialize() {
	/*
	 * Determine if the size of the system is smaller than the
	 * allowed potential cutoff length and throw an exception if
	 * it does.
	 */
	if (!integrator.potential.validateCutoff(halfLength)) {
		std::cout << "******** WARNING: LJCUT > L/2 ********" << std::endl;
		std::cout << "LJCUT = " << integrator.potential.cutoffBlowsAt
				<< ",   L/2  = " << halfLength << std::endl;
		std::cout << "******** STOPPING  SIMULATION ********" << std::endl;
		throw std::logic_error("Invalid LJCUT");
	}

	/*
	 * Insert the first atom randomly
	 */
	std::cout << "STARTING INITIALIZATION" << std::endl;
	// ASSIGN POSITION OF FIRST ATOM.
	Vector3d<T> p;
	p = randomVector();
	positions[0] = p;
	std::cout << "INSTERTION OF MOLECULE 1: " << positions[0] << std::endl;

	/*
	 * Then add the rest using the overlap check
	 */
	addPositions();
}

template<class T>
inline Vector3d<T> MolDyn::SimulationBox<T>::randomVector() {
	Vector3d<T> rdm;
	/*
	 * Generate the random vector with coordinates x,y,z
	 * in the range [-halfLength:halfLength]
	 */
	for (unsigned int i = 0; i < 3; ++i) {
		rdm[i] = distribution(generator) * length;
		rdm[i] -= halfLength;
	}
	return rdm;
}

template<class T>
inline void MolDyn::SimulationBox<T>::step(bool rescale = false) {
	// UPDATE POSITIONS AND VELOCITIES OF THE ATOMS
	integrator.advance(positions, velocities, forces, species, eKin, ePot, eTot,
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
	/*
	 * Add the positions of the remaining particles
	 */
	for (unsigned int i = 1; i < nParticles; ++i) {
		addAtom(i);
		std::cout.setf(std::ios::fixed, std::ios::floatfield);
		std::cout.precision(5);
		std::cout << "INSTERTION OF MOLECULE " << i + 1 << ": " << positions[i];
		std::cout << std::endl;
	}

	// CALCULATE FORCES
	integrator.forces(positions, forces, species, ePot, pressure, length);
	// RESCALE VELOCITIES TO (-1,1) & CALCULATE SUM OF X,Y,Z VELOCITY
	for (unsigned int i = 0; i < nParticles; ++i) {
		p = randomVector();
		velocities[i] = p * 2.0 - ones;
		sum += velocities[i];
	}
	// ZERO THE TOTAL LINEAR MOMENTUM BY REMOVING/ADDING MOMENTUM
	sum /= nParticles;
	for (unsigned int i = 0; i < nParticles; ++i) {
		velocities[i] -= sum;
	}
	// SCALE VELOCITIES TO SET-POINT TEMPERATURE
	integrator.rescale(velocities, species);
}

template<class T>
inline void MolDyn::SimulationBox<T>::addAtom(const unsigned int& ith) {
	Vector3d<T> p;
	bool retry = true;
	/*
	 * If the random position overlaps get a new random position
	 */
	while (retry == true) {
		p = randomVector();
		retry = overlaps(p, ith);
	}
	/*
	 * Else store the vector in the desired position
	 */
	positions[ith] = p;
}

template<class T>
inline void SimulationBox<T>::setSpecies(const unsigned int& idx,
		const unsigned int& type) {
	if (type >= maxComponents) {
		throw std::out_of_range("Trying to add an unspecified type of atom.");
	}
	if (idx >= nParticles) {
		throw std::out_of_range("Particle index out of range.");
	}
	species[idx] = type;
}

template<class T>
inline void SimulationBox<T>::setPosition(const unsigned int& idx,
		const Vector3d<T>& v) {
	if (idx >= nParticles) {
		throw std::out_of_range("Particle index out of range.");
	}
	positions[idx] = v;
}

template<class T>
inline void SimulationBox<T>::setVelocity(const unsigned int& idx,
		const Vector3d<T>& v) {
	if (idx >= nParticles) {
		throw std::out_of_range("Particle index out of range.");
	}
	velocities[idx] = v;
}

template<class T>
inline void SimulationBox<T>::setSpecies(const unsigned int& start,
		const unsigned int& end, const unsigned int& type) {
	if (type >= maxComponents) {
		throw std::out_of_range("Trying to add an unspecified type of atom.");
	}
	if (start >= nParticles || end >= nParticles) {
		throw std::out_of_range("Particle index out of range.");
	}

	unsigned int imin = (start <= end) ? start : end;
	unsigned int imax = (start <= end) ? end : start;
	for (unsigned int i = imin; i < imax; ++i) {
		species[i] = type;
	}
}

template<class T>
inline void SimulationBox<T>::setMass(const unsigned int& idx, const T& m) {
	if (idx >= maxComponents) {
		throw std::out_of_range("Particle index out of range.");
	}
	mass[idx] = m;
}

template<class T>
inline void SimulationBox<T>::addComponent(
		const MolDyn::Component<T>& _component) {
	if (componentCounter < maxComponents) {
		sigma[componentCounter] = _component.sigma;
		setSpecies(componentCounter,
				componentCounter + _component.nParticles - 1, componentCounter);
		++componentCounter;
	}
}

template<class T>
inline bool MolDyn::SimulationBox<T>::overlaps(const Vector3d<T>& v,
		const unsigned int& ith) {
	Vector3d<T> dr;
	bool result = false;
	T sigmaCut2 = sigmaCut * sigmaCut;
	T scut2;
	// CALCULATE DISTANCES FROM PREVIOUSLY INSTERTED ATOMS
	for (unsigned int i = 0; i < ith; ++i) {
		dr = v;
		dr -= positions[i];
		dr.modulo(length);
		/*
		 * If the (squared) magnitude is less than
		 * the (squared) sigma cutoff then the particle overlaps
		 * with a previously inserted one
		 */
		scut2 =  mixSigma(sigma[species[i]], sigma[species[ith]]);
		scut2 *= scut2*sigmaCut2;
		if ((dr * dr) < scut2) {
			result = true;
			break;
		}
	}
	return result;
}

template<class T>
inline MolDyn::SimulationBox<T>::~SimulationBox() {
	// Free memory
	std::vector<Vector3d<T> >().swap(positions);
	std::vector<Vector3d<T> >().swap(velocities);
	std::vector<Vector3d<T> >().swap(forces);

}

template<typename T>
inline T atomicDensity(const T& dens, const T& mass) {
	return dens * N0 * CM2A / mass;
}

template<typename T>
inline T getSimBoxLength(const unsigned int& atoms, const T& dens) {
	return pow(atoms / dens, 1.0 / 3.0);
}

template<typename T>
inline T getSimBoxLength(const std::vector<MolDyn::Component<T> >& components) {
	T _density;
	_density = 0;
	unsigned int totalParticles = 0;
	/*
	 * Assume a random distribution with all species mixed
	 * and average over all the densities
	 */
	for (unsigned int i = 0; i < components.size(); ++i) {
		_density += components[i].density;
		totalParticles += components[i].nParticles;
	}
	return getSimBoxLength(totalParticles, _density);
}

} /* namespace MolDyn */

#endif /* SIMULATIONBOX_HPP_ */
