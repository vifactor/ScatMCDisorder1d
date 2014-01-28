/*
 * MCMarkovCellConfigurationBuilder.h
 *
 *  Created on: 13 бер. 2013
 *      Author: kopp
 */

#ifndef MCMARKOVCELLCONFIGURATIONBUILDER_H_
#define MCMARKOVCELLCONFIGURATIONBUILDER_H_

#include "MCCalculator.h"
#include "Structures.h"
#include "MarkovChain.h"
#include "ElementScatteringFactor.h"

class CellScatterer1d : public MCScatterer1d
{
public:
	CellScatterer1d(const UnitCell1d * cell, double pos);
	virtual ~CellScatterer1d();
	virtual std::complex<double> getScatteringFactor(double q) const;
protected:
	static ElementScatteringFactorFactory m_ESFFactory;
	const UnitCell1d * m_cell;
	std::vector<const ElementScatteringFactor * > m_scatteringFactors;

	/*variables for acceleration*/
	mutable double m_prevQ;
	mutable std::complex<double> m_prevF;
};

class MCMarkovCellConfigurationBuilder: public MCConfigurationBuilder
{
public:
	MCMarkovCellConfigurationBuilder(const std::vector<UnitCell1d>& ucells, const gsl_matrix * Q, double totalLength);
	virtual ~MCMarkovCellConfigurationBuilder();
	virtual void build(std::vector<MCScatterer1d *>& cfg) const;
protected:
	double m_totalLength;
	std::vector<UnitCell1d> m_unitCells;
	MarkovChain * m_markovChain;
};

#endif /* MCMARKOVCELLCONFIGURATIONBUILDER_H_ */
