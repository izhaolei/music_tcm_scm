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
#ifndef CALCULATE_SPECTRUM_IMPL_H
#define CALCULATE_SPECTRUM_IMPL_H

#include <gr_complex.h>
#include <specest_api.h>

class SPECEST_API calculate_spectrum_impl
{
	public:
		
		virtual void caculate_tcm_spectrum(const gr_complexd (*data)[1024],
						unsigned int data_len,
						unsigned int spec_len,
						unsigned int snapshots,
						unsigned int sen_num,
						unsigned int src_num,
						double* pspectrum)=0;

		virtual void caculate_scm_spectrum(const gr_complexd (*data)[1024],
						unsigned int data_len,
						unsigned int spec_len,
						unsigned int snapshots,
						unsigned int sen_num,
						unsigned int src_num,
						double* pspectrum)=0;
#endif /* INCLUDED_SPECESTI_BEAMFORM_H */
