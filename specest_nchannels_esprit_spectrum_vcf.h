#ifndef INCLUDED_SPECEST_NCHANNELS_ESPRIT_SPECTRUM_VCF_H
#define INCLUDED_SPECEST_NCHANNELS_ESPRIT_SPECTRUM_VCF_H

#include <specest_api.h>
#include <gr_sync_block.h>
#include <specest_nchannels_esprit_impl.h>

class specest_nchannels_esprit_spectrum_vcf;
typedef boost::shared_ptr<specest_nchannels_esprit_spectrum_vcf> 
	specest_nchannels_esprit_spectrum_vcf_sptr;

SPECEST_API specest_nchannels_esprit_spectrum_vcf_sptr 
	specest_make_nchannels_esprit_spectrum_vcf (
	unsigned int n, unsigned int m, unsigned int nsamples,unsigned int spectrum_len);

/*!
 * \brief
 *
 */
class SPECEST_API specest_nchannels_esprit_spectrum_vcf : 
	public gr_sync_block{
	friend SPECEST_API specest_nchannels_esprit_spectrum_vcf_sptr 
	specest_make_nchannels_esprit_spectrum_vcf (unsigned int n, 
	unsigned int m, unsigned int nsamples,unsigned int spectrum_len);

	specest_nchannels_esprit_spectrum_vcf (unsigned int n, unsigned int m, 
	unsigned int nsamples,unsigned int spectrum_len);

 public:
	~specest_nchannels_esprit_spectrum_vcf ();

       
	int work (int noutput_items,
	          gr_vector_const_void_star &input_items,
	          gr_vector_void_star &output_items);

 private:
       unsigned int dm_nsamples;
	specest_nchannels_esprit_impl* dm_impl;
	unsigned int dm_n;
	unsigned int dm_m;
	unsigned int dm_spectrum_len;	
};

#endif /* INCLUDED_specest_beamform_nchannels_spectrum_vcf_H */
