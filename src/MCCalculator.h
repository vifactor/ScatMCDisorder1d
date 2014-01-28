/*
 * MCCalculator.h
 *
 *  Created on: 5 бер. 2013
 *      Author: kopp
 */

#ifndef MCCALCULATOR_H_
#define MCCALCULATOR_H_

#include <string>
#include <vector>
#include <complex>

class MCScatterer1d
{
public:
	MCScatterer1d(double pos) {m_pos = pos;}
	virtual ~MCScatterer1d() {}
	virtual double getCoordinate() const {return m_pos;};
	virtual std::complex<double> getScatteringFactor(double q) const = 0;
protected:
	double m_pos;
};

class MCConfigurationBuilder
{
public:
	MCConfigurationBuilder() {}
	virtual ~MCConfigurationBuilder() {}
	virtual void build(std::vector<MCScatterer1d *>& cfg) const = 0;
	virtual void clear(std::vector<MCScatterer1d *>& cfg) const
	{
		std::vector<MCScatterer1d *>::iterator it;

		for(it = cfg.begin(); it != cfg.end(); ++it)
		{
			delete (*it);
		}
		cfg.clear();
	}
protected:
};

class MCConfiguration : public std::vector<MCScatterer1d *>
{
public:
	MCConfiguration(const MCConfigurationBuilder * configurationBuilder)
	{
		m_configurationBuilder = configurationBuilder;
	}
	virtual ~MCConfiguration()
	{
		m_configurationBuilder->clear(*this);
	}
	void update()
	{
		m_configurationBuilder->clear(*this);
		m_configurationBuilder->build(*this);
	}
protected:
	const MCConfigurationBuilder * m_configurationBuilder;
};

class MCCalculator
{
public:
	MCCalculator(std::size_t nbMCSteps, MCConfiguration * configuration);
	virtual ~MCCalculator() {}
	double getIntensity(double q) const;
	void getIntensity(std::vector<double> qs, std::vector<double>& intensities) const;
protected:
	double makeStep(double q) const;

	MCConfiguration * m_configuration;
	std::size_t m_nbMCSteps;
};

#endif /* MCCALCULATOR_H_ */
