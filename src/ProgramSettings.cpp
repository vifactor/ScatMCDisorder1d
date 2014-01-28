/*
 * ProgramSettings.cpp
 *
 *  Created on: 7 бер. 2013
 *      Author: kopp
 */

#include "ProgramSettings.h"

Range readRange(const libconfig::Setting& range)
{
	Range temp;

	temp.m_min = range[0][0];
	temp.m_max = range[0][1];
	temp.m_sampling = range[1];

	return temp;
}

UnitCell1d readUnitCell1d(const libconfig::Setting& stg)
{
	UnitCell1d ucell;

	ucell.setLength(double(stg["length"]));

	const libconfig::Setting &atoms = stg["atoms"];

	if(atoms.isList())
	{
		size_t nbAtoms = atoms.getLength();
		//loop over atoms
		for(size_t i = 0; i < nbAtoms; ++i)
		{
			ucell.push_back(readAtomicSite1d(atoms[i]));
		}
	}
	else
	{
		//todo
	}

	return ucell;
}

AtomicSite1d readAtomicSite1d(const libconfig::Setting& stg)
{
	AtomicSite1d asite;
	if(stg.isList())
	{
		asite.m_type = stg[0].c_str();
		asite.m_relpos = stg[1];
		asite.m_p = stg[2];
		asite.m_sigma = stg[3];
	}
	else
	{
		//todo
	}

	return asite;
}

void printRange(const Range& r)
{
	LOG(logINFO) << "[" << r.m_min << "," << r.m_max <<"]" << ":" << r.m_sampling << std::endl;
}

void printUnitCell1d(const UnitCell1d& uc)
{
	LOG(logINFO) << "length = " << uc.getLength() << std::endl;

	for(size_t i = 0; i < uc.size(); ++i)
	{
		printAtomicSite1d(uc[i]);
	}
}

void printAtomicSite1d(const AtomicSite1d & as)
{
	LOG(logINFO) << "(" << as.m_type <<  ", " << as.m_relpos << ", " << as.m_p << ", " << as.m_sigma << ")" << std::endl;
}

ProgramSettings::ProgramSettings()
{
	sampleSettings.probabilityMatrix.nbRows = 0;
	sampleSettings.probabilityMatrix.nbColumns = 0;
	sampleSettings.probabilityMatrix.data = NULL;
}

ProgramSettings::~ProgramSettings()
{
	if(sampleSettings.probabilityMatrix.data)
		delete 	sampleSettings.probabilityMatrix.data;
}

void ProgramSettings::read(const std::string& cfgfile)
{
	libconfig::Config cfg;
	// Read the file. If there is an error, report it
	try
	{
		cfg.readFile(cfgfile.c_str());
		cfg.setAutoConvert(true);
		const libconfig::Setting& root = cfg.getRoot();

		readSampleSettings(root);
		readCalculatorSettings(root);
		readEngineSettings(root);
	} catch (const libconfig::FileIOException &fioex)
	{
		throw Exception(toString(fioex.what()) + " in\t" + cfgfile);
	} catch (const libconfig::ParseException &pex)
	{
		throw Exception(
				toString(pex.what()) + " in\t" + cfgfile + ":"
						+ toString(pex.getLine()) + " - "
						+ toString(pex.getError()));
	} catch (const libconfig::SettingNotFoundException &nfex)
	{
		throw Exception(
				toString(nfex.what()) + "\t" + toString(nfex.getPath())
						+ " in\t" + cfgfile);
	} catch (libconfig::SettingTypeException& tex)
	{
		throw Exception(
				toString(tex.what()) + "\t" + toString(tex.getPath()) + " in\t"
						+ cfgfile);
	}
}

void ProgramSettings::print() const
{
	LOG(logINFO) << "---Sample settings---" << std::endl;
	printSampleSettings();
	LOG(logINFO) << "---Calculator settings---" << std::endl;
	printCalculatorSettings();
	LOG(logINFO) << "---Engine settings---" << std::endl;
	printEngineSettings();
	LOG(logINFO) << "---------------------" << std::endl;
}

void ProgramSettings::readCalculatorSettings(const libconfig::Setting& root)
{
	const libconfig::Setting &calculator = root["Calculator"];

	calculatorSettings.nbMCSteps = calculator["nbMCSteps"];
	calculatorSettings.background = calculator["background"];
	calculatorSettings.scale = calculator["scale"];
	calculatorSettings.qrange = readRange(calculator["qrange"]);
}

void ProgramSettings::printCalculatorSettings() const
{
	LOG(logINFO) << "Nb MC Steps:\t" << calculatorSettings.nbMCSteps << std::endl;
	LOG(logINFO) << "qrange:" << std::endl;
	printRange(calculatorSettings.qrange);
}

void ProgramSettings::readSampleSettings(const libconfig::Setting& root)
{
	const libconfig::Setting &sample = root["Sample"];

	const libconfig::Setting &matrix = sample["probabilityMatrix"];
	readProbabilityMatrix(matrix);
	const libconfig::Setting &cells = sample["cells"];
	readCells(cells);
	const libconfig::Setting &length = sample["length"];
	sampleSettings.totalLength = length;
}

void ProgramSettings::readProbabilityMatrix(const libconfig::Setting& stg)
{
	size_t nbRows, nbColumns;
	double aij;

	/*check whether it is a matrix*/
	if(!stg.isList())
	{
		//todo
	}
	nbRows = stg.getLength();
	for(size_t irow = 0; irow < nbRows; ++irow)
	{
		if(!stg[irow].isList())
		{
			//todo
		}
		else
		{
			nbColumns = stg[irow].getLength();
			/*check whether it is a square matrix*/
			if(nbColumns != nbRows)
			{
				//todo
			}
		}
		for(size_t icol = 0; icol < nbColumns; ++icol)
		{
			aij = stg[irow][icol];
			if((aij < 0) || (aij > 1.0))
			{

			}
		}
	}

	sampleSettings.probabilityMatrix.nbRows = nbRows;
	sampleSettings.probabilityMatrix.nbColumns = nbColumns;

	sampleSettings.probabilityMatrix.data = new double[nbRows * nbColumns];
	for(size_t irow = 0; irow < nbRows; ++irow)
	{
		for(size_t icol = 0; icol < nbColumns; ++icol)
		{
			sampleSettings.probabilityMatrix.data[nbColumns * irow + icol] = stg[irow][icol];
		}
	}
}

void ProgramSettings::printProbabilityMatrix() const
{
	std::ostringstream os;
	size_t nbRows, nbColumns;

	nbRows = sampleSettings.probabilityMatrix.nbRows;
	nbColumns = sampleSettings.probabilityMatrix.nbColumns;

	for(size_t irow = 0; irow < nbRows; ++irow)
	{
		os << "||";
		for(size_t icol = 0; icol < nbColumns - 1; ++icol)
		{
			os << sampleSettings.probabilityMatrix.data[nbColumns * irow + icol] << "\t";
		}
		os << sampleSettings.probabilityMatrix.data[nbColumns * irow + (nbColumns - 1)] << "||\n";
	}
	LOG(logINFO) << "Probability matrix:\n" << os.str();
}

void ProgramSettings::readCells(const libconfig::Setting& stg)
{
	if(stg.isList())
	{
		/*loop over cells*/
		size_t nbCells = stg.getLength();
		if(nbCells != sampleSettings.probabilityMatrix.nbColumns)
		{
			//todo
		}
		for (size_t i = 0; i < nbCells; ++i)
		{
			sampleSettings.unitCells.push_back(readUnitCell1d(stg[i]));
		}
	}
	else
	{
		//todo
	}
}

void ProgramSettings::printCells() const
{
	for(size_t i = 0; i < sampleSettings.unitCells.size(); ++i)
	{
		LOG(logINFO) << "Cell "<< i << std::endl;
		printUnitCell1d(sampleSettings.unitCells[i]);
	}
}

void ProgramSettings::printSampleSettings() const
{
	printProbabilityMatrix();
	printCells();
	LOG(logINFO) << "Total length:\t" << sampleSettings.totalLength << std::endl;
}

void ProgramSettings::readEngineSettings(const libconfig::Setting& root)
{
	const libconfig::Setting &engine = root["Engine"];

	engineSettings.outfile = engine["outfile"].c_str();
}

void ProgramSettings::printEngineSettings() const
{
	LOG(logINFO) << "Output file:\t" << engineSettings.outfile << std::endl;
}
