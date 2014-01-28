/*
 * Structures.h
 *
 *  Created on: 5 бер. 2013
 *      Author: kopp
 */

#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include <string>
#include <vector>

struct AtomicSite1d
{
	AtomicSite1d(std::string type = "Dummy", double relpos = 0.0, double p = 1.0, double sigma = 0.0)
	{m_type = type; m_relpos = relpos; m_p = p; m_sigma = sigma;}
	/*element*/
	std::string m_type;
	/*coordinate*/
	double m_relpos;
	/*occupation probability*/
	double m_p;
	/*Debye-Waller factor*/
	double m_sigma;
};

class UnitCell1d: public std::vector<AtomicSite1d>
{
public:
	double getLength() const {return m_length;}
	void setLength(double length) {m_length = length;}
protected:
	double m_length;
};

struct Range
{
	double m_min;
	double m_max;
	size_t m_sampling;

	double getStep() const
	{
		return (m_sampling > 1) ? (m_max - m_min) / (m_sampling - 1) : 0.0;
	}
	void toVector(std::vector<double>& vec) const
	{
		double step = getStep();
		for (size_t i = 0; i < m_sampling; ++i)
		{
			vec.push_back(m_min + step * i);
		}
	}
	bool good() const
	{
		if ((m_min <= m_max) && (m_sampling > 0))
		{
			return true;
		}
		else
		{
			return false;
		}
	}
};

#endif /* STRUCTURES_H_ */
