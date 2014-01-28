/*
 * ElementScatteringFactor.cpp
 *
 *  Created on: 6 бер. 2013
 *      Author: kopp
 */

#include "ElementScatteringFactor.h"

F0DataBase::F0DataBase()
{
	F0CalculationData dummyCoeffs;

	for(size_t i = 0; i < F0CalculationData::NB_COEFS; ++i)
	{
		dummyCoeffs.m_A[i] = 0.0;
		dummyCoeffs.m_B[i] = 0.0;
	}
	dummyCoeffs.m_C = 0.0;

	m_dbFilename = "f0_WaasKirf.db";
	m_dbPath = ".\\DataBases\\";
	m_dbMap["Dummy"] = dummyCoeffs;
}

void F0DataBase::setPath(const std::string& path)
{
	m_dbPath = path;
	clear();
}

void F0DataBase::setFilename(const std::string& filename)
{
	m_dbFilename = filename;
	clear();
}

bool F0DataBase::readF0WaasmaierKirfel(const std::string& element)
{
	std::string filename;
	std::string line;
	std::ifstream fin;
	std::vector<std::string> substrings;
	F0CalculationData f0Coefs;

	filename = m_dbPath + m_dbFilename;

	fin.open(filename.c_str());
	if(!fin)
	{
		std::cout << "Warning: no database found!\t" << filename << std::endl;
		return false;
	}

	//read data
	while (true)
	{
		//check whether we reached the end of the file
		if(fin.eof())
		{
			std::cout << "Warning: no record for "<< element <<" found!" << filename << std::endl;
			fin.close();
			return false;
		}

		//read line
		getline(fin, line);
		//split line into parts
		substrings = split(line, "\t \n\r");

		if (substrings[0].compare("#S") == 0)
		{
			/*
			 * check the current element (ion):
			 *
			 * substrings[0] - '#S'
			 * substrings[1] - 'atomic number'
			 * substrings[2] - 'name of the element'
			 */
			std::string elname = substrings[2];

			if (element.compare(elname) == 0)
			{
				//skip 2 lines
				getline(fin, line);
				getline(fin, line);

				//read atomic parameters
				getline(fin, line);
				//std::cout<<"Text(data):\t"<< line <<std::endl;

				std::istringstream is(line);
				is >> f0Coefs.m_A[0] >> f0Coefs.m_A[1] >> f0Coefs.m_A[2]
						>> f0Coefs.m_A[3] >> f0Coefs.m_A[4] >> f0Coefs.m_C
						>> f0Coefs.m_B[0] >> f0Coefs.m_B[1] >> f0Coefs.m_B[2]
						>> f0Coefs.m_B[3] >> f0Coefs.m_B[4];
				m_dbMap[element] = f0Coefs;

				fin.close();
				return true;
			}
		}
	}
	fin.close();
	return false;
}

const F0CalculationData& F0DataBase::getCoeffs(const std::string& element)
{
	std::map<std::string, F0CalculationData>::const_iterator it;

	it = m_dbMap.find(element);

	if(it != m_dbMap.end())//element found in local database
	{
		return it->second;
	}
	else if(readF0WaasmaierKirfel(element))//element found in the database
	{
		it = m_dbMap.find(element);
		return it->second;
	}
	else//element not found
	{
		return m_dbMap["Dummy"];
	}
}

F1DataBase::F1DataBase()
{
	m_dbType = f1BRENNAN_COWAN;
	m_dbFilename = "f1f2_BrennanCowan.db";
	m_dbPath = ".\\DataBases\\";
}

F1DataBase::~F1DataBase()
{
	clear();
}

void F1DataBase::clear()
{
	std::map<std::string, F1CalculationData>::iterator it;

	for(it = m_dbMap.begin(); it != m_dbMap.end(); ++it)
	{
		gsl_interp_accel_free (it->second.m_acc_reF);
		gsl_interp_accel_free (it->second.m_acc_imF);
		gsl_spline_free (it->second.m_spline_reF);
		gsl_spline_free (it->second.m_spline_imF);
	}
	m_dbMap.clear();
}

void F1DataBase::setPath(const std::string& path)
{
	clear();
	m_dbPath = path;
}

void F1DataBase::setFilename(const std::string& filename, F1DataBaseType dbType)
{
	clear();
	m_dbFilename = filename;
	m_dbType = dbType;
}

void F1DataBase::updateMap(const std::string& element, std::vector<double>&en, std::vector<double>&ref, std::vector<double>&imf)
{
	F1CalculationData data;
	double *arr_en, *arr_ref, *arr_imf;

	arr_en = new double[en.size()];
	arr_ref = new double[ref.size()];
	arr_imf = new double[imf.size()];
	data.m_acc_reF = gsl_interp_accel_alloc ();
	data.m_acc_imF = gsl_interp_accel_alloc ();
	data.m_spline_reF = gsl_spline_alloc (gsl_interp_cspline, en.size());
	data.m_spline_imF = gsl_spline_alloc (gsl_interp_cspline, en.size());

	std::copy(en.begin(), en.end(), arr_en);
	std::copy(ref.begin(), ref.end(), arr_ref);
	std::copy(imf.begin(), imf.end(), arr_imf);

	gsl_spline_init (data.m_spline_reF, arr_en, arr_ref, en.size());
	gsl_spline_init (data.m_spline_imF, arr_en, arr_imf, en.size());

	m_dbMap[element] = data;

	delete[] arr_en;
	delete[] arr_imf;
	delete[] arr_ref;
}

bool F1DataBase::readF1Database(const std::string& element)
{
	bool result;

	switch (m_dbType)
	{
	case f1BRENNAN_COWAN:
		result = readF1BrennanCowan(element);
		break;
	case f1HENKE:
		result = readF1Henke(element);
		break;
	case f1SASAKI:
		result = readF1Sasaki(element);
		break;
	case f1SHADOW:
		result = readF1Shadow(element);
		break;
	case f1WINDT:
		result = readF1Windt(element);
		break;
	default:
		result = false;
		break;
	}

	return result;
}

bool F1DataBase::readF1BrennanCowan(const std::string& element)
{
	std::string filename;
	std::string line;
	std::ifstream fin;
	std::vector<std::string> substrings;

	filename = m_dbPath + m_dbFilename;
	fin.open(filename.c_str());

	if(!fin)
	{
		std::cout << "Warning: no database found!" << filename << std::endl;
		return false;
	}

	//read data
	while (true)
	{
		if(fin.eof())
		{
			std::cout << "Warning: no record for "<< element <<" found!" << std::endl;
			fin.close();
			return false;
		}
		//read line
		getline (fin, line);
		//split line into parts
		substrings = split(line, "\t \n");

		if(substrings[0].compare("#S") == 0)
		{
			/*
			 * check the current element (ion):
			 *
			 * substrings[0] - '#S'
			 * substrings[1] - 'atomic number'
			 * substrings[2] - 'name of the element'
			 */
			std::string elname = substrings[2];

			if(element.compare(elname) == 0)
			{
				std::istringstream is;
				double en, ref, imf;
				std::vector<double> vec_en, vec_ref, vec_imf;

				//skip 3 lines
				getline(fin, line);
				getline(fin, line);
				getline(fin, line);

				//read all energies and correction values for this element
				getline(fin, line);
				while((line[0] != '#') && (!fin.eof()))
				{
					is.str(line);
					is >> en >> ref >> imf;
					vec_en.push_back(en);
					vec_ref.push_back(ref);
					vec_imf.push_back(imf);
					is.clear();
					getline(fin, line);
				}
				updateMap(element, vec_en, vec_ref, vec_imf);
				fin.close();
				return true;
			}
		}
	}
	fin.close();
	return false;
}

bool F1DataBase::readF1Sasaki(const std::string& element)
{
	std::string filename;
	std::string line;
	std::ifstream fin;
	std::vector<std::string> substrings;

	filename = m_dbPath + m_dbFilename;
	fin.open(filename.c_str());

	if(!fin)
	{
		std::cout << "Warning: no database found!\t" << filename << std::endl;
		return false;
	}

	//read data
	while (true)
	{
		if(fin.eof())
		{
			std::cout << "Cannot find " << element<< " in " << filename << std::endl;
			break;
		}
		//read line
		getline (fin, line);
		//split line into parts
		substrings = split(line, "\t \n");

		if(substrings[0].compare("#S") == 0)
		{
			/*
			 * check the current element (ion):
			 *
			 * substrings[0] - '#S'
			 * substrings[1] - 'atomic number'
			 * substrings[2] - 'name of the element'
			 */
			std::string elname = substrings[2];

			if(element.compare(elname) == 0)
			{
				std::istringstream is;
				double atNumb;
				double en, ref, imf;
				std::vector<double> vec_en, vec_ref, vec_imf;

				//find atomic number
				is.str(substrings[1]);
				is >> atNumb;
				is.clear();

				//skip 4 lines
				for(size_t i = 0; i < 4; ++i) getline(fin, line);

				//read all energies and correction values for this element
				getline(fin, line);
				while(line[0] != '#' && (!fin.eof()))
				{
					is.str(line);
					is >> en >> ref >> imf;
					//corr.ref1 -= atNumb;
					vec_en.push_back(en);
					vec_ref.push_back(ref);
					vec_imf.push_back(imf);
					is.clear();
					getline(fin, line);
				}
				updateMap(element, vec_en, vec_ref, vec_imf);
				fin.close();
				return true;
			}
		}
	}
	fin.close();
	return false;
}

bool F1DataBase::readF1Henke(const std::string& element)
{
	std::string filename;
	std::string line;
	std::ifstream fin;
	std::vector<std::string> substrings;

	filename = m_dbPath + m_dbFilename;
	fin.open(filename.c_str());
	if(!fin)
	{
		std::cout << "Warning: no database found!\t" << filename << std::endl;
	}

	//read data
	while (true)
	{
		if(fin.eof())
		{
			std::cout << "Cannot find " << element << " in " << filename << std::endl;
			break;
		}
		//read line
		getline (fin, line);
		//split line into parts
		substrings = split(line, "\t \n");

		if(substrings[0].compare("#S") == 0)
		{
			/*
			 * check the current element (ion):
			 *
			 * substrings[0] - '#S'
			 * substrings[1] - 'atomic number'
			 * substrings[2] - 'name of the element'
			 */
			std::string elname = substrings[2];

			if(element.compare(elname) == 0)
			{
				std::istringstream is;
				double atNumb;
				double en, ref, imf;
				std::vector<double> vec_en, vec_ref, vec_imf;

				//find atomic number
				is.str(substrings[1]);
				is >> atNumb;
				is.clear();

				//skip 5 lines
				for(size_t i = 0; i < 5; ++i) getline(fin, line);

				//read all energies and correction values for this element
				getline(fin, line);
				while(line[0] != '#' && (!fin.eof()))
				{
					is.str(line);
					is >> en >> ref >> imf;
					ref -= atNumb;
					vec_en.push_back(en);
					vec_ref.push_back(ref);
					vec_imf.push_back(imf);
					is.clear();
					getline(fin, line);
				}
				updateMap(element, vec_en, vec_ref, vec_imf);
				fin.close();
				return true;
			}
		}
	}
	fin.close();
	return false;
}

bool F1DataBase::readF1Shadow(const std::string& element)
{
	std::string filename;
	std::string line;
	std::ifstream fin;
	std::vector<std::string> substrings;

	filename = m_dbPath + m_dbFilename;
	fin.open(filename.c_str());

	if(!fin)
	{
		std::cout << "Warning: no database found!\t" << filename << std::endl;
	}

	//read data
	while (true)
	{
		if(fin.eof())
		{
			std::cout << "Cannot find " << element << " in " << filename << std::endl;
			break;
		}
		//read line
		getline (fin, line);
		//split line into parts
		substrings = split(line, "\t \n");

		if(substrings[0].compare("#S") == 0)
		{
			/*
			 * check the current element (ion):
			 *
			 * substrings[0] - '#S'
			 * substrings[1] - 'atomic number'
			 * substrings[2] - 'name of the element'
			 */
			std::string elname = substrings[2];

			if(element.compare(elname) == 0)
			{
				std::istringstream is;
				double atNumb;
				double en, ref, imf;
				std::vector<double> vec_en, vec_ref, vec_imf;

				//find atomic number
				is.str(substrings[1]);
				is >> atNumb;
				is.clear();

				//skip 8 lines
				for(size_t i = 0; i < 8; ++i) getline(fin, line);

				//read all energies and correction values for this element
				getline(fin, line);
				while(line[0] != '#' && (!fin.eof()))
				{
					is.str(line);
					is >> en >> ref >> imf;
					ref -= atNumb;
					vec_en.push_back(en);
					vec_ref.push_back(ref);
					vec_imf.push_back(imf);
					is.clear();
					getline(fin, line);
				}
				updateMap(element, vec_en, vec_ref, vec_imf);
				fin.close();
				return true;
			}
		}
	}

	fin.close();
	return false;
}

bool F1DataBase::readF1Windt(const std::string& element)
{
	std::string filename;
	std::string line;
	std::ifstream fin;
	std::vector<std::string> substrings;

	filename = m_dbPath + m_dbFilename;
	fin.open(filename.c_str());
	if(!fin)
	{
		std::cout << "Warning: no database found!\t" << filename << std::endl;
	}

	//read data
	while (true)
	{
		if(fin.eof())
		{
			std::cout << "Cannot find " << element << " in " << filename << std::endl;
			break;
		}
		//read line
		getline (fin, line);
		//split line into parts
		substrings = split(line, "\t \n");

		if(substrings[0].compare("#S") == 0)
		{
			/*
			 * check the current element (ion):
			 *
			 * substrings[0] - '#S'
			 * substrings[1] - 'atomic number'
			 * substrings[2] - 'name of the element'
			 */
			std::string elname = substrings[2];

			if(element.compare(elname) == 0)
			{
				std::istringstream is;
				double atNumb;
				double en, ref, imf;
				std::vector<double> vec_en, vec_ref, vec_imf;

				//find atomic number
				is.str(substrings[1]);
				is >> atNumb;
				is.clear();

				//skip 11 lines
				for(size_t i = 0; i < 11; ++i) getline(fin, line);

				//read all energies and correction values for this element
				getline(fin, line);
				while(line[0] != '#' && (!fin.eof()))
				{
					is.str(line);
					is >> en >> ref >> imf;
					ref -= atNumb;
					vec_en.push_back(en);
					vec_ref.push_back(ref);
					vec_imf.push_back(imf);
					is.clear();
					getline(fin, line);
				}
				updateMap(element, vec_en, vec_ref, vec_imf);
				fin.close();
				return true;
			}
		}
	}

	fin.close();
	return false;
}

std::complex<double> F1DataBase::getF1(const std::string& element, double en)
{
	std::map<std::string, F1CalculationData>::const_iterator it;
	double re_F, im_F;

	it = m_dbMap.find(element);

	if(it != m_dbMap.end())
	{
		re_F = gsl_spline_eval (it->second.m_spline_reF, en, it->second.m_acc_reF);
		im_F = gsl_spline_eval (it->second.m_spline_imF, en, it->second.m_acc_imF);
	}
	else if(readF1Database(element))
	{
		it = m_dbMap.find(element);
		re_F = gsl_spline_eval (it->second.m_spline_reF, en, it->second.m_acc_reF);
		im_F = gsl_spline_eval (it->second.m_spline_imF, en, it->second.m_acc_imF);
	}
	else//element not found
	{
		re_F = 0.0;
		im_F = 0.0;
	}

	return std::complex<double>(re_F, im_F);
}

ElementScatteringFactor::ElementScatteringFactor(const std::string& element)
{
	m_elementName = element;
}

std::complex<double> ElementScatteringFactor::f0(double k) const
{
	static double result;
	const F0CalculationData& coeffs = F0DataBase::Instance().getCoeffs(m_elementName);

	result = coeffs.m_C;
	for (size_t i = 0; i < F0CalculationData::NB_COEFS; ++i)
	{
		result += coeffs.m_A[i] * exp(-coeffs.m_B[i] * k * k);
	}
	return result;
}

std::complex<double> ElementScatteringFactor::f1(double en) const
{
	return F1DataBase::Instance().getF1(m_elementName, en);
}

std::complex<double> ElementScatteringFactor::f(double k, double en) const
{
	return f0(k) + f1(en);
}

const ElementScatteringFactor * ElementScatteringFactorFactory::operator[](const std::string& type)
{
	std::map<std::string, ElementScatteringFactor *>::const_iterator it;

	it = m_scatteringFactors.find(type);
	if(it == m_scatteringFactors.end())
	{
		m_scatteringFactors[type] = new ElementScatteringFactor(type);
		return m_scatteringFactors[type];
	}
	else
	{
		return it->second;
	}
}

ElementScatteringFactorFactory::~ElementScatteringFactorFactory()
{
	std::map<std::string, ElementScatteringFactor *>::iterator it;
	for(it = m_scatteringFactors.begin(); it != m_scatteringFactors.end(); ++it)
	{
		delete it->second;
	}
}
