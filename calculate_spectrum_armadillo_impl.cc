/* -*- c++ -*- */
/*
 * Copyright 2010 Communications Engineering Lab, KIT
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */


#include <calculate_spectrum_armadillo_impl.h>
#include <specest_correst_impl.h>
#include <armadillo>
#include <complex>
#include <gr_complex.h>
#include <cmath>
#include <iostream>
using namespace arma;

calculate_spectrum_armadillo_impl::calculate_spectrum_armadillo_impl(unsigned int n, unsigned int m) : d_n(n), d_m(m)
{
	if (n > m)
			throw std::invalid_argument("specest_beamform_armadillo_impl: n cannot exceed m in length.");
}

calculate_spectrum_armadillo_impl::~calculate_spectrum_armadillo_impl()
{
}

void
calculate_spectrum_armadillo_impl::caculate_scm_spectrum(const gr_complexd (*data)[1024],
						unsigned int data_len,
						unsigned int spec_len,
						unsigned int snapshots,
						unsigned int sen_num,
						unsigned int src_num,
						double* pspectrum)
	arma::cx_mat R(sen_num, sen_num);
	specest_impl::multi_correst(data, data_len, sen_num, &R);
	arma::cx_mat Total_Sscm;
	arma::cx_vec Sscm;
	unsigned int i;
	arma::cx_vec v = arma::zeros<arma::cx_vec>(R.n_rows);
	for(i = 0; i < snapshots; i++){
		if (!arma::sum(1 - (R.col(i) == v))){
			Total_Sscm = Total_Sscm + 0;
		} else {
			Sscm = R.col(i) / norm(R.col(i), 2);
			Total_Sscm = Total_Sscm + Sscm * arma::trans(Sscm);
		}
	}

	arma::cx_mat Cscm = Total_Sscm / snapshots;

	calculate_spectrum(Cscm, spec_len, src_num, sen_num, pspectrum);
}

calculate_spectrum_armadillo_impl::caculate_tcm_spectrum(const gr_complexd (*data)[1024],
						unsigned int data_len,
						unsigned int spec_len,
						unsigned int snapshots,
						unsigned int sen_num,
						unsigned int src_num,
						double* pspectrum)
{
	arma::cx_mat Total_Sscm;
	arma::cx_vec v = arma::zeros<arma::cx_vec>(R.n_rows);
	arma::cx_mat R(sen_num, sen_num);
	specest_impl::multi_correst(data, data_len, sen_num, &R);
	arma::cx_vec Rtcm,Stcm;

	for(int i = 0; i < snapshots; i++){
		for(int j = 0; j < snapshots; j++){
			Rtcm = R.col(i) - R.col(j);
			if(!arma::sum(1 - (Rtcm == v))){
				Total_Sscm = Total_Sscm + 0;
			} else {
				Stcm = Rtcm / norm(Rtcm, 2);
				Total_Sscm = Total_Sscm + Stcm * arma::trans(Stcm);
			}
		}
	}

	arma::cx_mat Ctcm = Total_Sscm / (snapshots * (snapshots - 1));

	calculate_spectrum(Ctcm, spec_len, src_num, sen_num, pspectrum);
}

inline void calculate_spectrum(arma::cx_mat R,
						unsigned int spec_len,
						unsigned int src_num,
						unsigned int sen_num,
						double *pspectrum){
	//resu_vals <-> R_Ind
	//eigvec <-> R_V
	//eigvals <-> R_D(lambda)
	complex<double> pure_i(0, 1);
	arma::cx_mat eigvec, R_Gn;
	arma::mat tmpmat = arma::zeros<arma::mat>(R.n_rows, 1);
	arma::vec eigvals;
	arma::uvec resu_vals;

	arma::eig_sym(eigvals, eigvec, R);
	resu_vals = arma::sort_index(eigvals, "descend");

	arma::vec tmpvec = arma::zeros<arma::vec>(resu_vals.n_rows - src_num);
	for (unsigned int i = src_num; i < resu_vals.n_rows; i++){
		tmpvec(i - src_num, 0) = resu_vals(i, 0);
	}

	R_Gn = arma::zeros<arma::cx_mat>(eigvec.n_rows, tmpvec.n_rows);
	for (unsigned int i = 0; i < tmpvec.n_rows; i++){
		R_Gn.col(i) = eigvec.col(tmpvec(i, 0));
	}

	arma::cx_mat a_theata;
	arma::mat lin = arma::linspace<arma::mat>(0, sen_num - 1, sen_num);
	arma::vec est_theta = arma::linspace<arma::vec>(-90, 90, 181);
	for (int i = 0; i < spec_len; i++){
		a_theata = (arma::exp(-pure_i * arma::datum::pi * arma::trans(lin) * sin(est_theta(i) / 180 * arma::datum::pi)));
		arma::mat tmpmat = arma::abs(arma::conj(a_theata) * R_Gn * arma::trans(R_Gn) * a_theata.st());
		pspectrum[i] = 1 / tmpmat(0, 0);
	}
}
