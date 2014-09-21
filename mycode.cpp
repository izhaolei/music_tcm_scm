arma::cx_mat SCM_music(const gr_complexd(*data)[1024],unsigned int data_len,
	unsigned int snapshots)
{

	arma::cx_mat R(channels,channels);
	specest_impl::multi_correst(data,data_len,channels,&R);
	arma::cx_mat Total_Sscm;
	arma::cx_vec Sscm;
	unsigned int i;
	arma::cx_vec v = arma::zeros<arma::cx_vec>(R.n_rows);
	for(i=1;i<=snapshots;i++)
	{
		if (!arma::sum(1 - (R.col(i) == v)))
			Total_Sscm=Total_Sscm+0;
		else
		{
			Sscm=R.col(i)/norm(R.col(i),2);
			Total_Sscm=Total_Sscm+Sscm*arma::trans(Sscm);
		}
	}
	return Total_Sscm/snapshots;
}

arma::cx_mat TCM_music(const gr_complexd(*data)[1024],unsigned int data_len,
	unsigned int snapshots)
{
	arma::cx_mat Total_Sscm;
	unsigned int i,j;
	arma::cx_vec v = arma::zeros<arma::cx_vec>(R.n_rows);
	arma::cx_mat R(channels,channels);
	specest_impl::multi_correst(data,data_len,channels,&R);
	arma::cx_vec Rtcm,Stcm;
	for(i=1;i<=snapshots;i++)
		for(j=1;j<=snapshots;j++)
		{
			Rtcm=R.col(i)-R.col(j);
			if(!arma::sum(1 - (Rtcm == v)))
				Total_Sscm=Total_Sscm+0;
			else
			{
				Stcm=Rtcm/norm(Rtcm,2);
				Total_Sscm=Total_Sscm+Stcm*arma::trans(Stcm);
			}
		}
		return Total_Sscm/(snapshots*(snapshots-1));
}

void calculate(arma::cx_mat R, unsigned int src_num, unsigned int spec_len, unsigned int sen_num)
{
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
	arma::mat R_spec(1, spec_len);
	for (int i = 0; i < spec_len; i++){
		a_theata = (arma::exp(-pure_i * arma::datum::pi * arma::trans(lin) * sin(est_theta(i) / 180 * arma::datum::pi)));
		arma::mat tmpmat = arma::abs(arma::conj(a_theata) * R_Gn * arma::trans(R_Gn) * a_theata.st());
		R_spec(i) = 1 / tmpmat(0, 0);
	}
}