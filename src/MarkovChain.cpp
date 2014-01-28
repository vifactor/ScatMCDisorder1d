/*
 * MarkovChain.cpp
 *
 *  Created on: 5 бер. 2013
 *      Author: kopp
 */

#include "MarkovChain.h"
#include <iostream>

MarkovChain::MarkovChain(const gsl_matrix * Q)
{
	m_cdfQ = gsl_matrix_alloc(Q->size1, Q->size2);
	m_cdfS = gsl_vector_alloc(Q->size1);
	m_rng = gsl_rng_alloc (gsl_rng_taus);

	setupCDFS(Q);
	setupCDFQ(Q);
	setupRNG();
	update();
}

MarkovChain::~MarkovChain()
{
	gsl_rng_free(m_rng);
	gsl_matrix_free(m_cdfQ);
	gsl_vector_free(m_cdfS);
}

void MarkovChain::setupRNG()
{
	 gsl_rng_set (m_rng, time(0));
}

void MarkovChain::update() const
{
	double r;
	MarkovState state;

	state = 0;
	r = gsl_rng_uniform (m_rng);
	while(r > gsl_vector_get(m_cdfS, state))
	{
		++state;
	}
	m_state = state;
}

MarkovChain::MarkovState MarkovChain::getNextState() const
{
	double r;
	MarkovState nextState, thisState;

	thisState = m_state;
	nextState = 0;
	r = gsl_rng_uniform (m_rng);
	while (r > gsl_matrix_get(m_cdfQ, nextState, thisState))
	{
		++nextState;
	}
	m_state = nextState;

	return thisState;
}

void MarkovChain::setupCDFS(const gsl_matrix * Q)
{
	double cdf, norm;
	gsl_vector_complex *eval;
	gsl_matrix_complex *evec;
	gsl_eigen_nonsymmv_workspace * w;
	gsl_vector_complex_view S;

	gsl_matrix_memcpy (m_cdfQ, Q);
	eval = gsl_vector_complex_alloc (Q->size1);
	evec = gsl_matrix_complex_alloc (Q->size1, Q->size2);
	w = gsl_eigen_nonsymmv_alloc(Q->size1);

	gsl_eigen_nonsymmv (m_cdfQ, eval, evec, w);
    gsl_eigen_nonsymmv_sort (eval, evec,
                             GSL_EIGEN_SORT_ABS_DESC);

    /*vector of stationary probabilities corresponding to the eigenvalue 1 */
	S = gsl_matrix_complex_column(evec, 0);
	/*sum of vector elements*/
	norm = 0.0;
	for(size_t i = 0; i < Q->size1; ++i)
	{
		norm += GSL_REAL(gsl_vector_complex_get(&S.vector, i));
	}

	/*cdfs*/
	cdf = 0.0;
	for(size_t i = 0; i < Q->size1; ++i)
	{
		cdf += GSL_REAL(gsl_vector_complex_get(&S.vector, i)) / norm;
		gsl_vector_set(m_cdfS, i, cdf);
	}

	gsl_eigen_nonsymmv_free (w);
    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);
}

void MarkovChain::setupCDFQ(const gsl_matrix * Q)
{
	double cdf;

	/*filling of cdf matrix, which columns represent the discrete cdf*/
	for(size_t icolumn = 0; icolumn < Q->size2; ++icolumn)
	{
		cdf = 0.0;
		for(size_t irow = 0; irow < Q->size1; ++irow)
		{
			cdf += gsl_matrix_get (Q, irow, icolumn);
			gsl_matrix_set (m_cdfQ, irow, icolumn, cdf);
		}
	}
}
