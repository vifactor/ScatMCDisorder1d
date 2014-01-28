/*
 * MCMarkovConfigurationBuilder.h
 *
 *  Created on: 5 бер. 2013
 *      Author: kopp
 */

#ifndef MCMARKOVCONFIGURATIONBUILDER_H_
#define MCMARKOVCONFIGURATIONBUILDER_H_

#include "MCCalculator.h"
#include "Structures.h"
#include "MarkovChain.h"
#include "ElementScatteringFactor.h"

class AtomicScatterer1d : public MCScatterer1d
{
public:
	AtomicScatterer1d(const std::string& type, double pos, double p, double sigma = 0.0);
	virtual ~AtomicScatterer1d();
	virtual std::complex<double> getScatteringFactor(double q) const;
protected:
	static ElementScatteringFactorFactory m_ESFFactory;
	const ElementScatteringFactor * m_scatteringFactor;
	double m_sigma;
	double m_p;

	/*variables for acceleration*/
	mutable double m_prevQ;
	mutable std::complex<double> m_prevF;
};

class MCMarkovAtomicConfigurationBuilder: public MCConfigurationBuilder
{
public:
	MCMarkovAtomicConfigurationBuilder(const std::vector<UnitCell1d>& ucells, const gsl_matrix * Q, double totalLength);
	virtual ~MCMarkovAtomicConfigurationBuilder();
	virtual void build(std::vector<MCScatterer1d *>& cfg) const;
protected:
	double m_totalLength;
	std::vector<UnitCell1d> m_unitCells;
	MarkovChain * m_markovChain;
};

#endif /* MCMARKOVCONFIGURATIONBUILDER_H_ */
