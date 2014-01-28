//============================================================================
// Name        : Disorder1dMonteCarlo.cpp
// Author      : Viktor Kopp
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include "MCMarkovAtomicConfigurationBuilder.h"
#include "MCMarkovCellConfigurationBuilder.h"
#include "ElementScatteringFactor.h"
#include "ProgramSettings.h"
using namespace std;


double getPLGfactor(double theta);
double getTheta(double Q);

int main(int argc, char ** argv)
{
	OutputTO::Stream()=stdout;
	Log::ReportingLevel()=logDEBUG;

	ProgramSettings settings;

	try
	{
		settings.read("default.cfg");
	}
	catch(ProgramSettings::Exception& ex)
	{
		LOG(logERROR) << "Fail:\t" << ex.what();
		return -1;
	}

	settings.print();
	MCConfiguration * configuration;
	MCConfigurationBuilder * builder;
	MCCalculator * calculator;
	std::vector<double> qs, intensities;
	double background, scale;
	std::ofstream fout;

	gsl_matrix_const_view Qview = gsl_matrix_const_view_array(settings.getSampleSettings().probabilityMatrix.data,
			settings.getSampleSettings().probabilityMatrix.nbRows,
			settings.getSampleSettings().probabilityMatrix.nbColumns);

	builder = new MCMarkovCellConfigurationBuilder(settings.getSampleSettings().unitCells, &Qview.matrix, settings.getSampleSettings().totalLength);
	configuration = new MCConfiguration(builder);
	calculator = new MCCalculator(settings.getCalculatorSettings().nbMCSteps, configuration);

	/*range2vector*/
	settings.getCalculatorSettings().qrange.toVector(qs);
	background = settings.getCalculatorSettings().background;
	scale = settings.getCalculatorSettings().scale;

	LOG(logINFO) << "Calculating..." << std::endl;
	/*calculate*/
	calculator->getIntensity(qs, intensities);

	LOG(logINFO) << "Saving..." << std::endl;
	double plg;
	double theta;
	fout.open(settings.getEngineSettings().outfile.c_str());
	fout <<"2theta[deg]\tQ\tintensity\tPLG * intensity" << std::endl;
	for(double i = 0; i < qs.size(); ++i)
	{
		theta = getTheta(qs[i]);
		plg = getPLGfactor(theta);
		fout << 2 * theta * 180 / M_PI << "\t" << qs[i] << "\t" << scale * intensities[i] + background
				<< "\t" << scale * plg * intensities[i] + background << std::endl;
	}
	fout.close();

	delete calculator;
	delete configuration;
	delete builder;

	LOG(logINFO) << "Done." << std::endl;

	return 0;
}

double getPLGfactor(double theta)
{
	double sintheta = sin(theta);
	double costheta = cos(theta);
	double cos2theta = cos(2 *theta);

	return (1 + cos2theta * cos2theta)/(sintheta * sintheta * costheta);
}

double getTheta(double Q)
{
	const double lambda = 1.54;
	return asin(lambda * Q / (4 * M_PI));
}
