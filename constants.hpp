#ifndef CONSTANTS_HPP_
#define CONSTANTS_HPP_
#include <cmath>
namespace MolDyn {

/* PI */
static const double PI = 4*atan(1.0);
/* Avogadro's constant */
static const double N0 = 6.0221413E+23;
/* 1 cm^3 = 1.0 A^3 × 10−24 factor */
static const double CM2A = 1.0E-24;
/* 1 AMU = 1.660538921(73)×10−27 kg */
static const double AMU = 1.660538921E-27;
/* KB in amu, ps, A */
static const double KB = 0.83144621;

}

#endif /* CONSTANTS_HPP_ */
