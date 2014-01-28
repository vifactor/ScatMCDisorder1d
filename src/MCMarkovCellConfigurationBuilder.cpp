/*
 * MCMarkovCellConfigurationBuilder.cpp
 *
 *  Created on: 13 бер. 2013
 *      Author: kopp
 */

#include "MCMarkovCellConfigurationBuilder.h"

ElementScatteringFactorFactory CellScatterer1d::m_ESFFactory;

CellScatterer1d::CellScatterer1d(const UnitCell1d * cell, double pos) : MCScatterer1d(pos)
{
	m_cell = cell;
	for(size_t i = 0; i < m_cell->size(); ++i)
	{
		m_scatteringFactors.push_back(m_ESFFactory[(*m_cell)[i].m_type]);
	}

	m_prevQ = -1;
	m_prevF = 0.0;
}

CellScatterer1d::~CellScatterer1d()
{

}

std::complex<double> CellScatterer1d::getScatteringFactor(double q) const
{
	static std::complex<double> result, resulti;
	static const std::complex<double> I = std::complex<double>(0.0, 1.0);
	static double k;

	if (m_prevQ != q)
	{
		k = q / (4 * M_PI);

		result = 0.0;
		for(size_t i = 0; i < m_cell->size(); ++i)
		{
			/*scattering factor*/
			resulti = m_scatteringFactors[i]->f0(k);

			/*dispersion corrections*/
			/*TODO*/

			/*take an account for Debye-Waller factor*/
			resulti *= exp(-(*m_cell)[i].m_sigma * q * q);

			/*take into account occupation probability*/
			resulti *= (*m_cell)[i].m_p;

			result += resulti * exp(I * q * (*m_cell)[i].m_relpos * m_cell->getLength());
		}

		m_prevQ = q;
		m_prevF = result;
	}
	else
	{
		result = m_prevF;
	}

	return result;
}

MCMarkovCellConfigurationBuilder::MCMarkovCellConfigurationBuilder(const std::vector<UnitCell1d>& ucells, const gsl_matrix * Q, double totalLength)
{
	m_unitCells = ucells;
	m_markovChain = new MarkovChain(Q);
	m_totalLength = totalLength;
}

MCMarkovCellConfigurationBuilder::~MCMarkovCellConfigurationBuilder()
{
	if(m_markovChain)
		delete m_markovChain;
}

void MCMarkovCellConfigurationBuilder::build(std::vector<MCScatterer1d *>& cfg) const
{
	static double currentLength;
	MarkovChain::MarkovState state;

	currentLength = 0.0;
	m_markovChain->update();
	while(currentLength < m_totalLength)
	{
		state = m_markovChain->getNextState();
		std::cout << "state:\t" << state << std::endl;
		cfg.push_back(new CellScatterer1d(&m_unitCells[state], currentLength));
		currentLength += m_unitCells[state].getLength();
	}
}
