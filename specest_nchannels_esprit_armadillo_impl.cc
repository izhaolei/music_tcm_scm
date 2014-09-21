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


#include <specest_nchannels_esprit_armadillo_impl.h>
#include <specest_correst_impl.h>
#include <armadillo>
#include <complex>
#include <gr_complex.h>
#include <cmath>
#include <iostream>
using namespace arma;

specest_nchannels_esprit_armadillo_impl::specest_nchannels_esprit_armadillo_impl(unsigned int n, unsigned int m) : d_n(n), d_m(m)
{
	if (n > m)
			throw std::invalid_argument("specest_beamform_armadillo_impl: n cannot exceed m in length.");
}

specest_nchannels_esprit_armadillo_impl::~specest_nchannels_esprit_armadillo_impl()
{
}

void
specest_nchannels_esprit_armadillo_impl::calculate_nchannels_esprit_pspectrum(const gr_complexd(*data)[1024],
	unsigned int data_len,unsigned int channels,double*spectrum,unsigned int spectrum_len)
{
	arma::cx_mat R(channels,channels);
	arma::cx_mat G, G1, G2, G12(channels-1,2*d_n), U11(channels-d_n,channels-d_n), U21(channels-d_n,channels-d_n), phi_tls;
	arma::cx_vec search_steering_vector(channels);
	specest_impl::multi_correst(data,data_len,channels,&R);
	arma::vec eigvals, eigvals1;
	arma::cx_mat eigvec1, eigvec2, eigvec3;
	arma::cx_vec eigvals2;

        arma::eig_sym(eigvals,eigvec1,R);   

        G = eigvec1.cols(channels-d_n,channels-1);
	G1 = G.rows(0,channels-2);
	G2 = G.rows(1, channels-1);	

	G12.cols(0,d_n-1) = G1;
	G12.cols(d_n,2*d_n-1) = G2;

	arma::eig_sym(eigvals1,eigvec2,arma::trans(G12)*G12);

	U11 = eigvec2(span(0,d_n-1),span(0,d_n-1));

	U21 = eigvec2(span(d_n,2*d_n-1),span(0,d_n-1));

	phi_tls = (-1)*U11*arma::inv(U21);

	arma::eig_gen(eigvals2,eigvec3,phi_tls);
	
	for (int i=0;i<d_n;i++)
	{
   		spectrum[i]=std::asin(-std::arg(eigvals2[i])/(M_PI))*180/M_PI;
		std::cout<<"angle="<<spectrum[i]<<std::endl;
	}
}
