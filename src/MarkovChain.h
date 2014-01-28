/*
 * MarkovChain.h
 *
 *  Created on: 5 бер. 2013
 *      Author: kopp
 */

#ifndef MARKOVCHAIN_H_
#define MARKOVCHAIN_H_

#include <ctime>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>

class MarkovChain
{
public:
	typedef size_t MarkovState;
	MarkovChain(const gsl_matrix * Q);
	virtual ~MarkovChain();
	MarkovState getNextState() const;
	void update() const;
	const gsl_matrix * getCDF_Q() const {return m_cdfQ;}
	const gsl_vector * getCDF_S() const {return m_cdfS;}
protected:
	mutable MarkovState m_state;
	gsl_matrix * m_cdfQ;
	gsl_vector * m_cdfS;
	gsl_rng * m_rng;

	void setupRNG();
	void setupCDFS(const gsl_matrix * Q);
	void setupCDFQ(const gsl_matrix * Q);
};

#endif /* MARKOVCHAIN_H_ */
