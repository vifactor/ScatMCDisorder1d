/*
 * ElementScatteringFactor.h
 *
 *  Created on: 6 бер. 2013
 *      Author: kopp
 */

#ifndef ELEMENTSCATTERINGFACTOR_H_
#define ELEMENTSCATTERINGFACTOR_H_

#include <map>
#include <cmath>
#include <complex>
#include <string>
#include <fstream>
#include <StringTools.h>
#include <gsl/gsl_spline.h>

struct F0CalculationData
{
	enum {NB_COEFS = 5};
	/*
	 * f0[k] = c + [SUM a_i*EXP(-b_i*(k^2)) ]
	 *              i=1,5
	 *k = sin(theta) / lambda
	 */
	double m_A[NB_COEFS];
	double m_B[NB_COEFS];
	double m_C;
};

struct F1CalculationData
{
	gsl_interp_accel * m_acc_reF;
	gsl_interp_accel * m_acc_imF;
	gsl_spline * m_spline_reF;
	gsl_spline * m_spline_imF;
};

/*Database class here is singleton (Meyers)*/
class F0DataBase
{
public:
	static F0DataBase& Instance()
	{
		static F0DataBase theDataBase;
		return theDataBase;
	}
	void setPath(const std::string& path = "");
	void setFilename(const std::string& filename = "");
	const F0CalculationData& getCoeffs(const std::string& element);
protected:
	std::map<std::string, F0CalculationData> m_dbMap;
	std::string m_dbPath;
	std::string m_dbFilename;

	void clear() {m_dbMap.clear();}

	bool readF0WaasmaierKirfel(const std::string& element);
private:
	F0DataBase();								// Private constructor
	~F0DataBase() {clear();}					// Hidden destructor
	F0DataBase(const F0DataBase& );				// No copy construction
	F0DataBase& operator=(const F0DataBase& );	// No assignment
};

/*Database class here is singleton (Meyers)*/
class F1DataBase
{
public:
	enum F1DataBaseType {f1BRENNAN_COWAN, f1HENKE, f1SASAKI, f1SHADOW, f1WINDT};
	static F1DataBase& Instance()
	{
		static F1DataBase theDataBase;
		return theDataBase;
	}
	void setPath(const std::string& path = "");
	void setFilename(const std::string& filename, F1DataBaseType dbtype);
	std::complex<double> getF1(const std::string& element, double en);
protected:
	std::map<std::string, F1CalculationData> m_dbMap;

	std::string m_dbPath;
	std::string m_dbFilename;
	F1DataBaseType m_dbType;

	void clear();

	bool readF1Database(const std::string& element);

	bool readF1BrennanCowan(const std::string& element);
	bool readF1Sasaki(const std::string& element);
	bool readF1Henke(const std::string& element);
	bool readF1Shadow(const std::string& element);
	bool readF1Windt(const std::string& element);

	void updateMap(const std::string& element, std::vector<double>&en, std::vector<double>&ref, std::vector<double>&imf);
private:
	F1DataBase();								// Private constructor
	~F1DataBase();								// Hidden destructor
	F1DataBase(const F1DataBase& );				// No copy construction
	F1DataBase& operator=(const F1DataBase& );	// No assignment
};

class ElementScatteringFactor
{
public:
	ElementScatteringFactor(const std::string& element);
	std::complex<double> f0(double k) const;
	std::complex<double> f1(double en) const;
	std::complex<double> f(double k, double en) const;
	virtual ~ElementScatteringFactor() {}
protected:
	std::string m_elementName;
};

class ElementScatteringFactorFactory
{
public:
	ElementScatteringFactorFactory() {}
	~ElementScatteringFactorFactory();
	const ElementScatteringFactor * operator[](const std::string& type);
protected:
	std::map<std::string, ElementScatteringFactor *> m_scatteringFactors;
};

#endif /* ELEMENTSCATTERINGFACTOR_H_ */
