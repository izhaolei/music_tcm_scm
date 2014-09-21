#ifndef INCLUDED_CACULATE_SPECTRUM_VCF_H
#define INCLUDED_CACULATE_SPECTRUM_VCF_H

#include <specest_api.h>
#include <gr_sync_block.h>
#include <calculate_spectrum_impl.h>

class calculate_spectrum_vcf;
typedef boost::shared_ptr<calculate_spectrum_vcf> 
	calculate_spectrum_vcf_sptr;

SPECEST_API calculate_spectrum_vcf_sptr 
	make_calculate_spectrum_vcf (
	unsigned int n, unsigned int m, unsigned int nsamples,unsigned int spectrum_len);

/*!
 * \brief
 *
 */
class SPECEST_API calculate_spectrum_vcf : 
	public gr_sync_block{
	friend SPECEST_API calculate_spectrum_vcf_sptr 
	make_calculate_spectrum_vcf (unsigned int n, 
	unsigned int m, unsigned int nsamples,unsigned int spectrum_len);

	calculate_spectrum_vcf (unsigned int n, unsigned int m, 
	unsigned int nsamples,unsigned int spectrum_len);

 public:
	~calculate_spectrum_vcf ();

       
	int work (int noutput_items,
	          gr_vector_const_void_star &input_items,
	          gr_vector_void_star &output_items);

 private:
       unsigned int dm_nsamples;
	calculate_spectrum_impl* dm_impl;
	unsigned int dm_n;
	unsigned int dm_m;
	unsigned int dm_spectrum_len;	
};

#endif /* INCLUDED_specest_beamform_nchannels_spectrum_vcf_H */
