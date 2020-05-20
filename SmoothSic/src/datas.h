#ifndef FUNCTIONS_H
#define FUNCTIONS_H

# include <iostream>
# include <fstream>
# include <sstream>
# include <iomanip>

# include <gsl/gsl_sf.h>
# include <gsl/gsl_randist.h>
# include <gsl/gsl_cdf.h>
//# include <gsl/gsl_statistics_double.h>
//# include <gsl/gsl_fit.h>
//# include <gsl/gsl_multifit.h>
//#include <gsl/gsl_multifit_nlin.h>
//#include <gsl/gsl_math.h>
# include <math.h>
# include <cmath>
# include <vector>
# include <string>
# include <string.h>
# include <map>
#include <unistd.h>
# include <typeinfo>

using namespace std;



class datas
{

public:
	double  _prior_mg_mean=10.0;
	int _MAX_Z, _NMax;
	double _discretize;
	int ** _z,**_x, **_dk;
	double **_s;
	double *_pn;
	double _a;
	double apratio;
	bool _nb_flag, _zetak_flag, _kappan_flag;
	int ** _counts;

	double ** _predictiveposterior_cdf;
	double ** _predictiveposterior_pdf;
	double ** _predictiveposterior_tail;

	double **_kappa_n_comp; // 3x_N table comp=i,w,z
	int Ncomp=4;
	bool _window_adjustment_flag; // vector size 3 comp=i,w,z
	bool *_kappa_n_comp_flag; // vector size 3 comp=i,w,z
	double *_prior_a_kappa_n_comp;
	double *_prior_b_kappa_n_comp;
	int *_naccept_kappa_n_comp_moves;
	int *_nreject_kappa_n_comp_moves;
	long _nchanges_kg_move;
	long _nnochanges_kg_move;

	double *_kappa_g;
	bool _kappa_g_flag;

	int _naccept_kappag_move;
	int _nreject_kappag_move;

	int _naccept_an_move;
	int _nreject_an_move;

	double **_mu;
	double **_nu;
	int **_map_mu_nu;

	int  *_window;
	int **_kg;
	double** _gamma;
	int** _count_kg;
	double* _relweight_fcoolingtemp;
	double _coolingtemp;

	double* _ocur_pkn;// the prior on k for each cell n
	double _epsilon;
	int _GG, _G, _N, _K;
	float** _table_emis;
	float** _table_emis_an;
	double _prior_m_mu, _prior_b_mu,  _prior_d0, _prior_bg, _lambda_p,_lambda_a;
	double _var_kappa_g;

	double* _mg_mu;
	double _prior_b_shape, _prior_b_scale, _prior_mg_shape;

	double  _prior_zetap, _prior_zetaa, _prior_mp, _prior_ma, _prior_zetas;
	double* _prior_zetask;
	int  ** _temp;
	int* _sum_xk, *_sum_zn, *_sum_xn, *_sum_zk;
	double * _probas_k, *_sum_countk, *_sum_mun, *_probas_z, *_temp1, *_temp2;
	double *_an, *_bg;

	int* _d;/// for each gene indicator if it is the same among clustrers or not
	float * _fact_sum_zk; /// a vector used to store factorials for each cluster
	string* _cells, *_genes,*_genes_frag, *_clusters;

	double _mhstep;

	datas(int G, int N, int K);
	datas();
	virtual ~datas();

	void read_data(string file);
	void read_metadata(string file);
	void set_metadata(int G,int N, int K);
	void allocate_data_kg();
	void initialize_data_kg(gsl_rng* rgen);
	void set_flags(bool nbfl,bool zetafl,bool kappan_flag,bool kappa_i_flag,bool kappa_w_flag,bool kappa_w2_flag,bool kappa_z_flag, bool kappa_g_flag, bool window_adjustment_flag);
	double emission2 (int counts, int z, double a, double kappa_i, double kappa_w,double kappa_w2, int w, double kappa_z, double kappa_g ,bool verbose=false);
	double emission2_cdf (int counts, int z, double a, double kappa_i, double kappa_w,double kappa_w2, int w, double kappa_z, double kappa_g );
	void set_coolingtemp(double temp);

	void update_epsilon(gsl_rng* rgen);

	void update_z(gsl_rng* rgen);

	void update_kg(gsl_rng* rgen);
	void update_kg_z(gsl_rng* rgen);
	void update_kg_marginalized_mu(gsl_rng* rgen);
	void shuffle_kg(gsl_rng* rgen);
	void shuffle_kg_pair(gsl_rng* rgen);
	void update_gamma_kg(gsl_rng* rgen);


	void update_mu_nu(gsl_rng* rgen);
	void update_b(gsl_rng* rgen);
	void update_mg(gsl_rng* rgen);
	void update_mg_shape(gsl_rng* rgen);
	void deterministic_rescale_mg();
	void update_rescale_mg_mu_pn(gsl_rng* rgen, bool pn_allequal);

	void update_d(gsl_rng* rgen);
	void update_map_mu_nu_crp(gsl_rng* rgen);
	void update_param_crp(gsl_rng* rgen);
	double proba_config_crp(int *countsets, double theta);

	void update_prior_zetaa(gsl_rng* rgen);
	void update_prior_ma(gsl_rng* rgen);

	void update_pn(gsl_rng* rgen);
	void update_pn_allequal(gsl_rng* rgen);
	void update_prior_zetap(gsl_rng* rgen);
	void update_prior_mp(gsl_rng* rgen);

	void update_kappa_n_comp_an(gsl_rng* rgen, bool all_kappa_n_equal, bool all_an_equal, int nselgenes, int nloops,bool update_also_an);

	void update_kappa_g(gsl_rng* rgen);
	void update_var_kappa_g(gsl_rng* rgen);

	void update_counts_an(gsl_rng* rgen);
	double compute_likelihood();
	double* compute_likelihood_n();
	double* compute_likelihood_marginalized_z_s_n();
	double  compute_likelihood_marginalized_z_s();
	double compute_likelihood_marginalized_z();
	double compute_likelihood_marginalized_z_k(gsl_rng* rgen, bool predictiveposterior);

	void update_s(gsl_rng* rgen);
	void update_prior_zetas(gsl_rng* rgen);

	void update_prior_mg_mean(gsl_rng* rgen); // in principle should not be useful but introduced for validation purposes.
	void run_iteration(gsl_rng* rgen,int iter,int iter_allnequal,int iter_update_kappa,string* updates,int n_updates,bool update_also_an);
	void update_an_pn_ratio(gsl_rng* rgen,int iter_allnequal,int iter_update_kappa);
	void update_an_pn_ratio_n(gsl_rng* rgen, int iter_allnequal,int iter_update_kappa);
};

void simulate_data (gsl_rng* rgen, datas* dat);
#endif // FUNCTIONS_H
