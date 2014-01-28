/*
 * MCMarkovConfigurationBuilder.cpp
 *
 *  Created on: 5 mar. 2013
 *      Author: kopp
 */

#include "MCMarkovAtomicConfigurationBuilder.h"

ElementScatteringFactorFactory AtomicScatterer1d::m_ESFFactory;

AtomicScatterer1d::AtomicScatterer1d(const std::string& type, double pos, double p, double sigma) : MCScatterer1d(pos)
{
	m_scatteringFactor = m_ESFFactory[type];
	m_sigma = sigma;
	m_p = p;

	m_prevQ = -1;
	m_prevF = 0.0;
}

AtomicScatterer1d::~AtomicScatterer1d()
{
}

std::complex<double> AtomicScatterer1d::getScatteringFactor(double q) const
{
	static std::complex<double> result;
	static double k;

	if(m_prevQ != q)
	{
		/*scattering factor*/
		k = q / (4 * M_PI);

		result = m_scatteringFactor->f0(k);

		/*dispersion corrections*/
		/*TODO*/

		/*take an account for Debye-Waller factor*/
		result *= exp(-m_sigma * q * q);

		/*take into account occupation probability*/
		result *= m_p;

		m_prevQ = q;
		m_prevF = result;
	}
	else
	{
		result = m_prevF;
	}

	return result;
}

MCMarkovAtomicConfigurationBuilder::MCMarkovAtomicConfigurationBuilder(
		const std::vector<UnitCell1d>& ucells, const gsl_matrix * Q, double totalLength)
{
	m_unitCells = ucells;
	m_markovChain = new MarkovChain(Q);
	m_totalLength = totalLength;
}

MCMarkovAtomicConfigurationBuilder::~MCMarkovAtomicConfigurationBuilder()
{
	if(m_markovChain)
		delete m_markovChain;
}

void MCMarkovAtomicConfigurationBuilder::build(std::vector<MCScatterer1d *>& cfg) const
{
	static double length, currentLength;
	MarkovChain::MarkovState state;
	UnitCell1d::const_iterator it;

	currentLength = 0.0;
	m_markovChain->update();
	while(currentLength < m_totalLength)
	{
		state = m_markovChain->getNextState();
		length = m_unitCells[state].getLength();
		for(it = m_unitCells[state].begin(); it != m_unitCells[state].end(); ++it)
		{
			cfg.push_back(new AtomicScatterer1d(it->m_type, currentLength + it->m_relpos * length, it->m_p, it->m_sigma));
		}
		currentLength += length;
	}
}

