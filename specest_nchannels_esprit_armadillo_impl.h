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

#include <specest_api.h>
#include <specest_nchannels_esprit_impl.h>
#include <armadillo>

class SPECEST_API specest_nchannels_esprit_armadillo_impl : virtual public specest_nchannels_esprit_impl
{
	public:
		specest_nchannels_esprit_armadillo_impl(unsigned int n, unsigned int m);
		~specest_nchannels_esprit_armadillo_impl();

		void calculate_nchannels_esprit_pspectrum(const gr_complexd (*data)[1024],
				unsigned int data_len,unsigned int channels,
				double*spectrum,unsigned int spectrum_len);

	private:
		unsigned int d_n;
		unsigned int d_m;

};
