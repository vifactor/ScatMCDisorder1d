/*
 * ProgramSettings.h
 *
 *  Created on: 7 mar. 2013
 *      Author: kopp
 */

#ifndef PROGRAMSETTINGS_H_
#define PROGRAMSETTINGS_H_

#include "Log.h"
#include <libconfig.h++>
#include "Structures.h"
#include <StringTools.h>

Range readRange(const libconfig::Setting& stg);
UnitCell1d readUnitCell1d(const libconfig::Setting& stg);
AtomicSite1d readAtomicSite1d(const libconfig::Setting& stg);
void printRange(const Range& r);
void printUnitCell1d(const UnitCell1d& uc);
void printAtomicSite1d(const AtomicSite1d & as);

class ProgramSettings
{
public:
	class Exception: public std::exception
	{
	public:
		Exception(std::string m)
		{
			msg = "Configuration::" + m;
		}
		~Exception() throw ()
		{
		}
		const char* what() const throw ()
		{
			return msg.c_str();
		}
	private:
		std::string msg;
	};
	struct SampleSettings
	{
		struct
		{
			size_t nbRows;
			size_t nbColumns;
			double * data;
		} probabilityMatrix;

		std::vector<UnitCell1d> unitCells;
		double totalLength;
	};
	struct CalculatorSettings
	{
		size_t nbMCSteps;
		Range qrange;
		double background;
		double scale;
	};
	struct EngineSettings
	{
		std::string outfile;
	};
	ProgramSettings();
	virtual ~ProgramSettings();
	const SampleSettings& getSampleSettings() const
	{
		return sampleSettings;
	}
	const CalculatorSettings& getCalculatorSettings() const
	{
		return calculatorSettings;
	}
	const EngineSettings& getEngineSettings() const
	{
		return engineSettings;
	}
	void read(const std::string& cfg);
	void print() const;
protected:
	void readCalculatorSettings(const libconfig::Setting& root);
	void readSampleSettings(const libconfig::Setting& root);
	void readEngineSettings(const libconfig::Setting& root);

	void readProbabilityMatrix(const libconfig::Setting& stg);
	void readCells(const libconfig::Setting& stg);

	void printSampleSettings() const;
	void printCalculatorSettings() const;
	void printEngineSettings() const;

	void printProbabilityMatrix() const;
	void printCells() const;

	SampleSettings sampleSettings;
	CalculatorSettings calculatorSettings;
	EngineSettings engineSettings;
};

#endif /* PROGRAMSETTINGS_H_ */
