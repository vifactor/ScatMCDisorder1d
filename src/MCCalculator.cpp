/*
 * MCCalculator.cpp
 *
 *  Created on: 5 бер. 2013
 *      Author: kopp
 */

#include "MCCalculator.h"
#include <iostream>

MCCalculator::MCCalculator(std::size_t nb, MCConfiguration * configuration)
{
	m_nbMCSteps = nb;
	m_configuration = configuration;
}

double MCCalculator::getIntensity(double q) const
{
	static double result;

	result = 0.0;
	for(size_t i = 0; i < m_nbMCSteps; ++i)
	{
		m_configuration->update();
		result += makeStep(q);
	}
	result /= m_nbMCSteps;

	return result;
}

void MCCalculator::getIntensity(std::vector<double> qs,
		std::vector<double>& intensities) const
{
	intensities.resize(qs.size(), 0.0);

	for (size_t i = 0; i < m_nbMCSteps; ++i)
	{
		m_configuration->update();
		std::cout << i << "-th configuration size:\t" << m_configuration->size() << std::endl;
		for (size_t i = 0; i < qs.size(); ++i)
		{
			intensities[i] += makeStep(qs[i]);
		}
	}
	for (size_t i = 0; i < qs.size(); ++i)
	{
		intensities[i] /= m_nbMCSteps;
	}
}

double MCCalculator::makeStep(double q) const
{
	static std::complex<double> result;
	/*imaginary unity*/
	static const std::complex<double> I = std::complex<double>(0.0, 1.0);

	result = 0.0;
	/*offdiagonal element*/
	for(size_t i = 0; i < m_configuration->size(); ++i)
	{
		for(size_t j = i + 1; j < m_configuration->size(); ++j)
		{
			result += 2.0 * (*m_configuration)[i]->getScatteringFactor(q) * conj((*m_configuration)[j]->getScatteringFactor(q))
					* exp(I * q * ((*m_configuration)[i]->getCoordinate() - (*m_configuration)[j]->getCoordinate()));
		}
	}

	/*diagonal element*/
	for(size_t i = 0; i < m_configuration->size(); ++i)
	{
		result += (*m_configuration)[i]->getScatteringFactor(q) * conj((*m_configuration)[i]->getScatteringFactor(q));
	}

	return result.real();
}
