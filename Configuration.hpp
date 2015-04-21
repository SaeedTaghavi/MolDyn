/*
 * Configuration.hpp
 *
 * Parses the configuration XML file
 *
 *  Created on: Apr 19, 2015
 *      Author: Erick R Martinez Loran
 */

#ifndef CONFIGURATION_HPP_
#define CONFIGURATION_HPP_
#include <string>
#include <vector>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "Component.hpp"

/* Default configuration values */
#define LJCUT 		2.5
#define DT 			0.002
#define TEMP		150
#define EQSTEPS		10000
#define PRODSTEPS	50000
#define IENERGY		1
#define ITAPE		100
#define IANIMATE	100
#define OVERLAP		0.9

namespace MolDyn {

template<class T>
class Configuration {
public:
	bool writePositions;
	bool writeVelocities;
	bool writeMovie;
	bool writeEnergies;
	bool writePressure;
	bool writeDensityDistribution;
	unsigned int maxComponents;
	unsigned int totalParticles;
	std::string strTapePosFilename;
	std::string strTapeVelFilename;
	std::string strMovieFilename;
	std::string strEnergiesFilename;
	std::string strPressureFilename;
	std::string strDensityDistFilename;
	unsigned int outputFileCount;
	std::vector<MolDyn::Component<T> > components;
	T temperature;
	T dt;
	T LJCutoff;
	T sigmaCut;
	unsigned int equilibrationSteps;
	unsigned int productionSteps;
	unsigned int energyInterval;
	unsigned int tapeInterval;
	unsigned int animationInterval;
public:
	/*
	 * Loads simulation_settings structure from the specified XML file
	 * @param[in]	The XML filename
	 */
	Configuration(const std::string &filename);
	virtual ~Configuration();
};

template<class T>
inline Configuration<T>::Configuration(const std::string &filename) :
		outputFileCount(0), writePositions(false), writeVelocities(false), writeEnergies(
				false), writePressure(false), writeMovie(false), writeDensityDistribution(
				false), totalParticles(0) {
	// Create an empty property tree object
	using boost::property_tree::ptree;
	ptree pt;

	// Define a component counter;
	unsigned int _compCount = 0;

	T sigma, epsilon, mass, density, nparticles;
	std::string identifier;

	// Load the XML file into the property tree. If reading fails
	// (cannot open file, parse error), an exception is thrown.
	read_xml(filename, pt);

	maxComponents = 0;//pt.get<unsigned int>("config.species");

	/* Try to read every output file name */
	try {
		strTapePosFilename = pt.get<std::string>("config.files.positions");
		writePositions = true;
		++outputFileCount;

	} catch (std::exception& e) {
		std::cout << "Did not find 'files.positions'" << std::endl;
	}

	try {
		strTapeVelFilename = pt.get<std::string>("config.files.velocities");
		writeVelocities = true;
		++outputFileCount;

	} catch (std::exception& e) {
		std::cout << "Did not find 'files.velocities'" << std::endl;
	}

	try {
		strMovieFilename = pt.get<std::string>("config.files.movie");
		writeMovie = true;
		++outputFileCount;

	} catch (std::exception& e) {
		std::cout << "Did not find 'files.movie'" << std::endl;
	}

	try {
		strEnergiesFilename = pt.get<std::string>("config.files.energies");
		writeEnergies = true;
		++outputFileCount;

	} catch (std::exception& e) {
		std::cout << "Did not find 'files.energies'" << std::endl;
	}

	try {
		strPressureFilename = pt.get<std::string>("config.files.pressure");
		writePressure = true;
		++outputFileCount;

	} catch (std::exception& e) {
		std::cout << "Did not find 'files.pressure'" << std::endl;
	}

	try {
		strDensityDistFilename = pt.get<std::string>(
				"config.files.densityDistribution");
		writeDensityDistribution = true;
		++outputFileCount;

	} catch (std::exception& e) {
		std::cout << "Did not find 'files.densityDistribution'" << std::endl;
	}

	/* Load each component values into a vector */
	for (auto& comp : pt.get_child("config.components")) {
		if (comp.first == "component") {
			sigma = comp.second.get<T>("sigma");
			epsilon = comp.second.get<T>("epsilon");
			nparticles = comp.second.get<unsigned int>("particles");
			mass = comp.second.get<T>("mass");
			density = comp.second.get<T>("density");
			identifier = comp.second.get<std::string>("identifier");
			components.push_back(
					MolDyn::Component<T>(sigma, epsilon, mass, density,
							nparticles, identifier));
			totalParticles += nparticles;
			++maxComponents;
		}
	}

	/* Read the general parameters */
	temperature = pt.get("config.parameters.temperature", TEMP);
	dt = pt.get("config.parameters.dt", DT);
	LJCutoff = pt.get("config.parameters.LJCut", LJCUT);
	equilibrationSteps = pt.get("config.parameters.equilibrationSteps", EQSTEPS);
	productionSteps = pt.get("config.parameters.productionSteps", PRODSTEPS);
	energyInterval = pt.get("config.parameters.energyInterval", IENERGY);
	tapeInterval = pt.get("config.parameters.tapeInterval", ITAPE);
	animationInterval = pt.get("config.parameters.animationInterval", IANIMATE);
	sigmaCut = pt.get("config.parameters.overlap",OVERLAP);
}

template<class T>
inline MolDyn::Configuration<T>::~Configuration() {
	// Free memory
	std::vector<MolDyn::Component<T> >().swap(components);
}

} /* namespace MolDyn */

#endif /* CONFIGURATION_HPP_ */
