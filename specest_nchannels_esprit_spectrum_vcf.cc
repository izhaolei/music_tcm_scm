#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include <gr_io_signature.h>
#include <specest_nchannels_esprit_spectrum_vcf.h>
#include <specest_nchannels_esprit_armadillo_impl.h>
#include <iostream>
#include <math.h>


specest_nchannels_esprit_spectrum_vcf_sptr
specest_make_nchannels_esprit_spectrum_vcf (unsigned int n, 
	unsigned int m, unsigned int nsamples,unsigned int spectrum_len)
{
	return specest_nchannels_esprit_spectrum_vcf_sptr (
		new specest_nchannels_esprit_spectrum_vcf (
		n, m, nsamples,spectrum_len));
}

specest_nchannels_esprit_spectrum_vcf::specest_nchannels_esprit_spectrum_vcf (
	unsigned int n, unsigned int m, unsigned int nsamples,unsigned int spectrum_len)
	: gr_sync_block ("nchannels_esprit_spectrum_vcf",
		gr_make_io_signature (1, -1, sizeof (gr_complex) * nsamples),
		gr_make_io_signature (1, 1, sizeof (float)*n)),
		dm_m(m),
		dm_n(n),
		dm_nsamples(nsamples),
		dm_spectrum_len(spectrum_len),
		dm_impl(new specest_nchannels_esprit_armadillo_impl(n,m))
{				  
}

specest_nchannels_esprit_spectrum_vcf::~specest_nchannels_esprit_spectrum_vcf ()
{
	delete dm_impl;
}

int
specest_nchannels_esprit_spectrum_vcf::work (int noutput_items,
                                  gr_vector_const_void_star &input_items,
                                  gr_vector_void_star &output_items)
{  
	unsigned int channels=input_items.size();
      
      	const gr_complex *in[channels];//指针数组,存储数据通道指针
      	float* out = static_cast<float*> (output_items[0]);
      	double* tmpout = new double[dm_n];   

      	if(channels>1)
         {
            gr_complexd double_input[channels][1024];//最大采样单元1024*4
            for(int k = 0; k < dm_nsamples; ++k)
                  {
                  for(int j=0;j<channels;++j)
		         {
                       double_input[j][k]=((gr_complex *) input_items[j])[k];
                         }    
                   } 
           
            dm_impl->calculate_nchannels_esprit_pspectrum(
			double_input,dm_nsamples,channels, 
			tmpout, dm_spectrum_len);
        }
       	   		
	for(int i = 0; i < dm_spectrum_len; i++)
		{
                out[i] = (float) tmpout[i];	
                }
       float max=out[0];
       for(int m=1; m<dm_spectrum_len; m++)
		{
                if (out[m]>=max) max=out[m];	
                }
       for(int l = 0; l< dm_spectrum_len; l++)
		{
                out[l]=out[l]/max;
                }
	delete[] tmpout;//删除内存数据

	return noutput_items;
}
