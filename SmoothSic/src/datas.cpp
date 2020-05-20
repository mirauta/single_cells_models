
#include "datas.h"
#include "additional_functions.h"


datas::datas(){ }

datas::~datas() {}

double datas::emission2 (int counts, int z, double a,double kappa_i, double kappa_w, double kappa_w2, int w, double kappa_z, double kappa_g,bool verbose){


	if(z==0&&counts==0) return 1;
	else if (z==0&&counts>0) {
		//		if (verbose) cout<<"z=0 :"<<counts<<" "<<z<<endl;
		return 0;
	}
	else {
		double overdisp=(kappa_i + kappa_w/w+ kappa_w2*w+kappa_z/(z+1))*kappa_g;
		if(overdisp<0.00001)
			return gsl_ran_poisson_pdf(counts, w*a*z);
		if (_window_adjustment_flag)
			return gsl_ran_negative_binomial_pdf(counts, 1.0/(w*a*z*overdisp+1.0), 1.0/overdisp);
		else
			return gsl_ran_negative_binomial_pdf(counts, 1.0/(a*z*overdisp+1.0), 1.0/overdisp);
	}
}

double datas::emission2_cdf (int counts, int z, double a,double kappa_i, double kappa_w, double kappa_w2, int w, double kappa_z, double kappa_g){
	if(z==0&&counts==0) return 1;
	else if (z==0&&counts>0) return 0;
	else {
		double overdisp=(kappa_i + kappa_w/w+ kappa_w2*w+kappa_z/(z+1))*kappa_g;
		if (_window_adjustment_flag)
			return gsl_cdf_negative_binomial_P(counts, 1.0/(w*a*z*overdisp+1.0), 1.0/overdisp);
		else
			return gsl_cdf_negative_binomial_P(counts, 1.0/(a*z*overdisp+1.0), 1.0/overdisp);
	}
}

void datas::read_metadata(string file){

	ifstream file1;
	file1.open(file.c_str());
	if(!file1) { 
		cerr<<"Did not found the file: "<< file <<endl;
		exit(EXIT_FAILURE);
	}
	file1>>_GG;file1>>_G; file1>>_N;
	file1.close();
	cout<<"Read metadata. No. genes:"<<_G<<" No. cells: "<<_N<<endl;
}

void datas::set_metadata( int G,int N, int K){
	this->_G=G;
	this->_N=N;
	this->_K=K;
}

void datas::read_data(string file){

	string dummy;
	int temp;
	ifstream file1;
	file1.open(file.c_str());
	if(!file1) {
		cerr<<"Did not found the file: "<< file <<endl;
		exit(EXIT_FAILURE);
	}

	file1>>dummy; file1>>dummy;file1>>dummy;	file1>>dummy; file1>>dummy;file1>>dummy;
	for(int n=0; n<_N; n++){file1>>_cells[n];}
	for(int n=0; n<_N; n++){cout<<_cells[n]<<" ";}
	cout<<endl;
	for(int g=0; g<_G; g++){
		file1>>_genes[g];
		file1>>_genes_frag[g];
		file1>>temp;
		_window[g]=ceil(log(temp));
		for(int n=0; n<_N; n++){file1>>_counts[g][n];}
		//		cout<<";"<<g<<":"<<_genes[g]<<": "<<_genes_frag[g]<<"/ "<<temp<<endl;
		//		for(int n=0; n<_N; n++){cout<<_counts[g][n]<<" ";}
		//		cout<<endl;
	}
	file1.close();
	cout<<"Finished reading the counts. No. genes:"<<_G<<" No. cells: "<<_N<<endl;




	//	exit (EXIT_FAILURE);
}



void datas::update_an_pn_ratio_n(gsl_rng* rgen, int iter_allnequal,int iter_update_kappa){



	/* keep old values:*/
	double *oldan= new double [_N];
	double *oldpn= new double [_N];
	int** oldz = new int*[_G];
	//	double oldprior_zetap=_prior_zetap, oldprior_zetaa=_prior_zetaa, oldprior_mp=_prior_mp, oldprior_ma=_prior_ma ;

	for (int n=0; n<_N;n++){
		oldan[n]=_an[n];
		oldpn[n]=_pn[n];
	}
	for (int g=0; g<_G;g++) {
		oldz[g]= new int[_N];
		for (int n=0; n<_N;n++)
			oldz[g][n]=_z[g][n];
	}

	/*keep the old likelihood*/
	double* oldlike= compute_likelihood_marginalized_z_s_n() ;
	cout<<"OLD likelihood_marg zs: "<<oldlike[0]<<endl;
	/* sample  a proposed rescaling of a  and update an and pn:*/
	double proposed= (rand() > RAND_MAX/2) ? 1.1 : 0.9;
	cout<<"proposed: "<<proposed<<endl;

	for (int j=0;j<5;j++)		cout<<_an[j]<<" ";	cout<<endl;
	for (int n=0; n<_N;n++){
		_an[n]=_an[n]*proposed;
		_pn[n]=_pn[n]/proposed;	}
	for (int j=0;j<5;j++)		cout<<_an[j]<<" ";	cout<<endl;
	/*run 100 interations to update the other parameters:*/
	string  updates[]={"update_z"};
	//	bool update_also_an=false;
	//	int iter_allnequal=0;
	//	int iter_update_kappa=1;

	for (int i=99999;i<99999+100;i++)
		run_iteration(rgen,i,0,1,updates,sizeof(updates)/sizeof(updates[0]),false);

	//	for (int i=0;i<100;i++){
	//		string  updates[]={"update_b","update_mg","update_mg_shape","update_rescale_mg_mu_pn","update_pn","update_prior_mp","update_prior_zetap",\
	//				"update_z","update_prior_ma","update_prior_zetaa","update_kappa_an","update_s"};
	//		//run_iteration(rgen,i,iter_allnequal,iter_update_kappa,updates,sizeof(updates)/sizeof(updates[0]));
	//	}
	double* newlike=compute_likelihood_marginalized_z_s_n() ;
	cout<<"NEW likelihood"<<newlike[0]<<endl;
	cout<<"COMPARE likelihood" <<endl;
	for (int j=0;j<5;j++)		cout<<oldlike [j]<<" ";	cout<<endl;
	for (int j=0;j<5;j++)		cout<<newlike [j]<<" ";	cout<<endl;
	cout<<"ACCEPTED: ";
	for (int n=0; n<_N;n++){
		if (oldlike[n]>newlike[n]){
			for (int g=0; g<_G;g++)
				_z[g][n]=oldz[g][n];
			_an[n]=oldan[n];
			_pn[n]=oldpn[n];
		}
		else	cout<<n<<"; ";
	}
	cout<<endl;
	for (int j=0;j<5;j++)		cout<<_an[j]<<" ";	cout<<endl;
	cout<<"<<<<"<<endl;
	//	exit(EXIT_FAILURE);
}

void datas::update_an_pn_ratio(gsl_rng* rgen, int iter_allnequal,int iter_update_kappa){



	/* keep old values:*/
	double *oldan= new double [_N];
	double *oldpn= new double [_N];
	int** oldz = new int*[_G];
	int** olds = new int*[_G];
	double** oldmu = new double*[_G];
	double** oldkappa_n_comp = new double*[3];
	double* oldkappa_g = new double[_G];
	double oldprior_b_shape=_prior_b_shape, oldprior_b_scale=_prior_b_scale, oldprior_mg_shape=_prior_mg_shape, oldprior_mg_mean=_prior_mg_mean, oldprior_zetap=_prior_zetap, oldprior_zetaa=_prior_zetaa, oldprior_mp=_prior_mp, oldprior_ma=_prior_ma, oldprior_zetas=_prior_zetas;

	for (int icomp=0; icomp<3; icomp++) {
		oldkappa_n_comp[icomp] = new double [_N];
		for (int n=0; n<_N;n++)
			oldkappa_n_comp[icomp][n]=_kappa_n_comp[icomp][n];

	}

	for (int g=0; g<_G;g++){
		oldmu[g]=new double[_K];
		oldkappa_g[g]=_kappa_g[g];
		for (int k=0; k<_K;k++)
			oldmu[g][k]=_mu[g][k];
	}
	for (int n=0; n<_N;n++){
		oldan[n]=_an[n];
		oldpn[n]=_pn[n];
	}
	for (int g=0; g<_G;g++) {
		oldz[g]= new int[_N];
		for (int n=0; n<_N;n++)
			oldz[g][n]=_z[g][n];
	}
	for (int g=0; g<_G;g++) {
		olds[g]= new int[_N];
		for (int n=0; n<_N;n++)
			olds[g][n]=_s[g][n];
	}
	/*keep the old likelihood*/
	double oldlike=compute_likelihood();
	cout<<"OLD likelihood: "<<oldlike<<endl;
	double oldlike_zs= compute_likelihood_marginalized_z_s();
	cout<<"OLD likelihood_marg zs: "<<oldlike_zs<<endl;



	/* sample  a proposed rescaling of a  and update an and pn:*/

	double proposed= (rand() > RAND_MAX/2) ? 1.1 : 0.9;
	cout<<"proposed: "<<proposed<<endl;

	for (int j=0;j<5;j++)		cout<<_an[j]<<" ";	cout<<endl;
	for (int n=0; n<_N;n++){		_an[n]=_an[n]*proposed;		_pn[n]=_pn[n]/proposed;	}
	for (int j=0;j<5;j++)		cout<<_an[j]<<" ";	cout<<endl;
	/*run 100 interations to update the other parameters:*/
	for (int i=99999;i<99999+200;i++){
		string  updates[]={"update_b","update_mu_nu","update_mg","update_mg_shape","update_prior_mp","update_prior_zetap",\
				"update_z","update_prior_ma","update_prior_zetaa","update_kappa_an","update_s","update_pn"};
		bool update_also_an=false;
		int iter_allnequal=0;
		int  iter_update_kappa=1;

		run_iteration(rgen,i,iter_allnequal,iter_update_kappa,updates,sizeof(updates)/sizeof(updates[0]),update_also_an);
	}
	//	for (int i=0;i<100;i++){
	//		string  updates[]={"update_b","update_mg","update_mg_shape","update_rescale_mg_mu_pn","update_pn","update_prior_mp","update_prior_zetap",\
	//				"update_z","update_prior_ma","update_prior_zetaa","update_kappa_an","update_s"};
	//		//run_iteration(rgen,i,iter_allnequal,iter_update_kappa,updates,sizeof(updates)/sizeof(updates[0]));
	//	}
	double newlike=compute_likelihood();
	double newlike_zs= compute_likelihood_marginalized_z_s();
	cout<<"NEW likelihood"<<newlike<<endl;
	cout<<"NEW likelihood_marg zs: "<<newlike_zs<<endl;
	//	exit(EXIT_FAILURE);
	if (oldlike>newlike){
		_prior_b_shape=oldprior_b_shape, _prior_b_scale=oldprior_b_scale, _prior_mg_shape=oldprior_mg_shape, _prior_mg_mean=oldprior_mg_mean, _prior_zetap=oldprior_zetap, _prior_zetaa=oldprior_zetaa, \
				_prior_mp=oldprior_mp, _prior_ma=oldprior_ma, _prior_zetas=oldprior_zetas;

		for (int g=0; g<_G;g++){
			_kappa_g[g]=oldkappa_g[g];
			for (int k=0; k<_K;k++)
				_mu[g][k]=oldmu[g][k];
			for (int n=0; n<_N;n++){
				_z[g][n]=oldz[g][n];
				_s[g][n]=olds[g][n];
			}
		}
		for (int n=0; n<_N;n++){
			_an[n]=oldan[n];
			_pn[n]=oldpn[n];
		}

		for (int icomp=0; icomp<3; icomp++) {
			for (int n=0; n<_N;n++)
				_kappa_n_comp[icomp][n]=oldkappa_n_comp[icomp][n];

		}
		string  updates[]={"update_b"};
		run_iteration(rgen,1,iter_allnequal,iter_update_kappa,updates,sizeof(updates)/sizeof(updates[0]),false);
	}
	else {

		cout<<"ACCEPTED\n"<<endl;
	}
	for (int j=0;j<5;j++)		cout<<_an[j]<<" ";	cout<<endl;
	cout<<"<<<<"<<endl;
	//	exit(EXIT_FAILURE);
}

void datas::update_kappa_n_comp_an(gsl_rng* rgen, bool all_kappa_n_equal, bool all_an_equal, int nselgenes, int nloops,bool update_also_an){

	double *kappa_n_comp_prop= new double [Ncomp];
	double *logpobs_old = new double [_N];
	double *logpobs_prop = new double [_N];
	int* src=new int [_G];
	for (int i=0; i<_G; i++) {
		src[i]=i;
	}
	int* selgenes = new int [nselgenes]; // we propose to select only a subset of all the genes to compute MH acceptation ratios
	gsl_ran_choose(rgen,selgenes,nselgenes,src,_G,sizeof(int));
	//cout << "selgenes done" << endl;

	for (int n=0; n<_N;n++){
		logpobs_old[n]=0;
		for (int ig=0; ig<nselgenes;ig++){
			int g=selgenes[ig];
			logpobs_old[n] += log(emission2(_counts[g][n],_z[g][n],_an[n],_kappa_n_comp[0][n],_kappa_n_comp[1][n],_kappa_n_comp[3][n],_window[g],_kappa_n_comp[2][n],_kappa_g[g] ));
		}
	}
	for (int iloops=0; iloops<nloops; iloops++) {
		//update kappa_n_comp
		for (int icomp=0; icomp<Ncomp; icomp++) {
			if (_kappa_n_comp_flag[icomp]) {
				if (all_kappa_n_equal || !_kappan_flag) {
					for (int jcomp=0; jcomp<Ncomp; jcomp++) {
						if (jcomp==icomp) {
							kappa_n_comp_prop[jcomp]=gsl_ran_gamma(rgen, 1.0/pow(_mhstep,2.0), _kappa_n_comp[jcomp][0]*pow(_mhstep,2.0));
						} else {
							kappa_n_comp_prop[jcomp]=_kappa_n_comp[jcomp][0];
						}
					}
					double proba_prop=log(gsl_ran_gamma_pdf(_kappa_n_comp[icomp][0],1.0/pow(_mhstep,2.0),kappa_n_comp_prop[icomp]*pow(_mhstep,2.0)))+log(gsl_ran_gamma_pdf(kappa_n_comp_prop[icomp],_prior_a_kappa_n_comp[icomp],_prior_b_kappa_n_comp[icomp]));
					double proba_old=log(gsl_ran_gamma_pdf(kappa_n_comp_prop[icomp],1.0/pow(_mhstep,2.0),_kappa_n_comp[icomp][0]*pow(_mhstep,2.0)))+log(gsl_ran_gamma_pdf(_kappa_n_comp[icomp][0],_prior_a_kappa_n_comp[icomp],_prior_b_kappa_n_comp[icomp]));
					double logratioobs=0;
					for (int n=0; n<_N;n++){
						logpobs_prop[n]=0;
						for (int ig=0; ig<nselgenes;ig++){
							int g=selgenes[ig];
							logpobs_prop[n] += log(emission2(_counts[g][n],_z[g][n],_an[n],kappa_n_comp_prop[0],kappa_n_comp_prop[1],kappa_n_comp_prop[3],_window[g],kappa_n_comp_prop[2],_kappa_g[g]));
						}
						logratioobs+=logpobs_prop[n]-logpobs_old[n];
					}
					if(log(gsl_ran_flat(rgen,0,1))<proba_prop-proba_old+logratioobs) {
						for (int n=0; n<_N;n++) {
							_kappa_n_comp[icomp][n]=kappa_n_comp_prop[icomp];
							logpobs_old[n]=logpobs_prop[n];
						}
						_naccept_kappa_n_comp_moves[icomp]++;
					} else {
						_nreject_kappa_n_comp_moves[icomp]++;
					}
				} else { //!all_kappa_n_equal
					for (int n=0; n<_N;n++){
						for (int jcomp=0; jcomp<Ncomp; jcomp++) {
							if (jcomp==icomp) {
								kappa_n_comp_prop[jcomp]=gsl_ran_gamma(rgen, 1.0/pow(_mhstep,2.0), _kappa_n_comp[jcomp][n]*pow(_mhstep,2.0));;
							} else {
								kappa_n_comp_prop[jcomp]=_kappa_n_comp[jcomp][n];
							}
						}
						double proba_prop=log(gsl_ran_gamma_pdf(_kappa_n_comp[icomp][n],1.0/pow(_mhstep,2.0),kappa_n_comp_prop[icomp]*pow(_mhstep,2.0)))+log(gsl_ran_gamma_pdf(kappa_n_comp_prop[icomp],_prior_a_kappa_n_comp[icomp],_prior_b_kappa_n_comp[icomp]));
						double proba_old=log(gsl_ran_gamma_pdf(kappa_n_comp_prop[icomp],1.0/pow(_mhstep,2.0),_kappa_n_comp[icomp][n]*pow(_mhstep,2.0)))+log(gsl_ran_gamma_pdf(_kappa_n_comp[icomp][n],_prior_a_kappa_n_comp[icomp],_prior_b_kappa_n_comp[icomp]));
						logpobs_prop[n]=0;
						for (int ig=0; ig<nselgenes;ig++){
							int g=selgenes[ig];
							logpobs_prop[n] += log(emission2(_counts[g][n],_z[g][n],_an[n],kappa_n_comp_prop[0],kappa_n_comp_prop[1],kappa_n_comp_prop[3],_window[g],kappa_n_comp_prop[2],_kappa_g[g]));
						}
						double logratioobs = logpobs_prop[n]-logpobs_old[n];
						if(log(gsl_ran_flat(rgen,0,1))<proba_prop-proba_old+logratioobs) {
							_kappa_n_comp[icomp][n]=kappa_n_comp_prop[icomp];
							logpobs_old[n]=logpobs_prop[n];
							_naccept_kappa_n_comp_moves[icomp]++;
						} else {
							_nreject_kappa_n_comp_moves[icomp]++;
						}
					}
				}
			}
		}
		//update an
		if (update_also_an){
			if (all_an_equal) {

				double proposed=gsl_ran_gamma(rgen, 1.0/pow(_mhstep,2.0), _an[0]*pow(_mhstep,2.0));
				double proba_prop=log(gsl_ran_gamma_pdf(_an[0],1.0/pow(_mhstep,2.0),proposed*pow(_mhstep,2.0)))+log(gsl_ran_gamma_pdf(proposed,_prior_zetaa,_prior_ma/_prior_zetaa));
				double proba_old=log(gsl_ran_gamma_pdf(proposed,1.0/pow(_mhstep,2.0),_an[0]*pow(_mhstep,2.0)))+log(gsl_ran_gamma_pdf(_an[0],_prior_zetaa,_prior_ma/_prior_zetaa));
				double logratioobs=0;
				for (int n=0; n<_N;n++){
					logpobs_prop[n]=0;
					for (int ig=0; ig<nselgenes;ig++){
						int g=selgenes[ig];
						logpobs_prop[n]+=log(emission2(_counts[g][n],_z[g][n],proposed,_kappa_n_comp[0][n],_kappa_n_comp[1][n],_kappa_n_comp[3][n],_window[g],_kappa_n_comp[2][n],_kappa_g[g]));
					}
					logratioobs+=logpobs_prop[n]-logpobs_old[n];
				}
				if(log(gsl_ran_flat(rgen,0,1))<proba_prop-proba_old+logratioobs) {
					for (int n=0; n<_N;n++){
						_an[n]=proposed;
						logpobs_old[n]=logpobs_prop[n];
					}
					_naccept_an_move++;
				} else {
					_nreject_an_move++;
				}
			} else {
				for (int n=0; n<_N;n++){
					double proposed=gsl_ran_gamma(rgen, 1.0/pow(_mhstep,2.0), _an[n]*pow(_mhstep,2.0));
					double proba_prop=log(gsl_ran_gamma_pdf(_an[n],1.0/pow(_mhstep,2.0),proposed*pow(_mhstep,2.0)))+log(gsl_ran_gamma_pdf(proposed,_prior_zetaa,_prior_ma/_prior_zetaa));
					double proba_old=log(gsl_ran_gamma_pdf(proposed,1.0/pow(_mhstep,2.0),_an[n]*pow(_mhstep,2.0)))+log(gsl_ran_gamma_pdf(_an[n],_prior_zetaa,_prior_ma/_prior_zetaa));
					logpobs_prop[n]=0;
					for (int ig=0; ig<nselgenes;ig++){
						int g=selgenes[ig];
						logpobs_prop[n]+=log(emission2(_counts[g][n],_z[g][n],proposed,_kappa_n_comp[0][n],_kappa_n_comp[1][n],_kappa_n_comp[3][n],_window[g],_kappa_n_comp[2][n],_kappa_g[g]));
					}
					double logratioobs=logpobs_prop[n]-logpobs_old[n];

					if(log(gsl_ran_flat(rgen,0,1))<proba_prop-proba_old+logratioobs) {
						_an[n]=proposed;
						logpobs_old[n]=logpobs_prop[n];
						_naccept_an_move++;
					} else {
						_nreject_an_move++;
					}
				}
			}
		}
	}
	delete [] src;
	delete [] selgenes;
	delete [] kappa_n_comp_prop;
	delete [] logpobs_old;
	delete [] logpobs_prop;
}  

void datas::update_kappa_g(gsl_rng* rgen){
	if (!_kappa_g_flag) {
		return;
	}
	//cout<<"AA"<<endl;
	_naccept_kappag_move=1000;
	for (int g=0; g<_G;g++){
		double kappa=_kappa_g[g];
		double proposed=gsl_ran_gamma(rgen,1.0/pow(_mhstep,2.0), kappa*pow(_mhstep,2.0));
		double proba_prop=log(gsl_ran_gamma_pdf(kappa,1.0/pow(_mhstep,2.0),proposed*pow(_mhstep,2.0)))+log(gsl_ran_gamma_pdf(proposed,1.0/_var_kappa_g,_var_kappa_g));
		double proba_old=log(gsl_ran_gamma_pdf(proposed,1.0/pow(_mhstep,2.0),kappa*pow(_mhstep,2.0)))+log(gsl_ran_gamma_pdf(kappa,1.0/_var_kappa_g,_var_kappa_g));
		double logratioobs=0;
		for (int n=0; n<_N;n++){
			logratioobs+= log(emission2(_counts[g][n],_z[g][n],_an[n],_kappa_n_comp[0][n],_kappa_n_comp[1][n],_kappa_n_comp[3][n],_window[g],_kappa_n_comp[2][n],proposed))-log(emission2(_counts[g][n],_z[g][n],_an[n],_kappa_n_comp[0][n],_kappa_n_comp[1][n],_kappa_n_comp[3][n],_window[g],_kappa_n_comp[2][n],_kappa_g[g]));
		}

		if(log(gsl_ran_flat(rgen,0,1))<proba_prop-proba_old+logratioobs) {
			_kappa_g[g]=proposed;
			_naccept_kappag_move++;
		} else {
			_nreject_kappag_move++;
		}
	}
}

void datas::update_var_kappa_g(gsl_rng* rgen){
	if (!_kappa_g_flag) {
		return;
	}

	double var_proposed=gsl_ran_gamma(rgen,1.0/pow(_mhstep,2.0), _var_kappa_g*pow(_mhstep,2.0));

	double proba_prop=log(gsl_ran_gamma_pdf(_var_kappa_g,1.0/pow(_mhstep,2.0),var_proposed*pow(_mhstep,2.0)))-var_proposed;
	double proba_old=log(gsl_ran_gamma_pdf(var_proposed,1.0/pow(_mhstep,2.0),_var_kappa_g*pow(_mhstep,2.0)))-_var_kappa_g;
	double logratioobs=0;
	for (int g=0; g<_G;g++){
		logratioobs+=log(gsl_ran_gamma_pdf(_kappa_g[g],1.0/var_proposed,var_proposed))-log(gsl_ran_gamma_pdf(_kappa_g[g],1.0/_var_kappa_g,_var_kappa_g));
	}
	if(log(gsl_ran_flat(rgen,0,1))<proba_prop-proba_old+logratioobs) {
		_var_kappa_g=var_proposed;
	}
}

void datas::update_z(gsl_rng* rgen) {  
	for (int n=0; n<_N;n++){
		for (int g=0; g<_G;g++){
			int proposed=gsl_ran_poisson(rgen,_z[g][n]+1);
			//if (g==9 && n==144) {
			//cout << "update_z " << _z[g][n] << "->" << proposed << endl;
			//}
			if (proposed!=_z[g][n]) {
				double proba_prop= gsl_ran_poisson_pdf(_z[g][n], proposed+1)*emission2(_counts[g][n],proposed,_an[n],_kappa_n_comp[0][n],_kappa_n_comp[1][n],_kappa_n_comp[3][n],_window[g],_kappa_n_comp[2][n],_kappa_g[g])*\
						gsl_ran_poisson_pdf(proposed,_pn[n]*_s[g][n]*_mu[g][_kg[g][n]]);
				double proba_old=gsl_ran_poisson_pdf(proposed, _z[g][n]+1)*emission2(_counts[g][n],_z[g][n],_an[n],_kappa_n_comp[0][n],_kappa_n_comp[1][n],_kappa_n_comp[3][n],_window[g],_kappa_n_comp[2][n],_kappa_g[g])*\
						gsl_ran_poisson_pdf(_z[g][n],_pn[n]*_s[g][n]*_mu[g][_kg[g][n]]);

				if (proba_prop>proba_old) {
					_z[g][n]=proposed;
				} else {
					double paccept=proba_prop/proba_old;
					double r=gsl_ran_flat(rgen,0,1);
					if (r<paccept)
						_z[g][n]=proposed;
				}
			}
		}
	}
}


void datas::update_pn(gsl_rng* rgen){
	for (int n=0; n<_N;n++){
		double sum_zn=0;
		double sum_mun=0;
		if(_nb_flag)
			for (int g=0; g<_G;g++){
				sum_zn+=_z[g][n];
				sum_mun+=_mu[g][_kg[g][n]]*_s[g][n];
			}
		else
			for (int g=0; g<_G;g++){
				sum_zn+=_z[g][n];
				sum_mun+=_mu[g][_kg[g][n]];
			}
		_pn[n]=gsl_ran_gamma(rgen, _prior_zetap+sum_zn,1.0/(sum_mun+ _prior_zetap/_prior_mp));
	}
}


void datas::update_pn_allequal(gsl_rng* rgen){

	double sum_z=0, sum_mu=0;
	for (int n=0; n<_N;n++){
		double t_sum_z=0;
		double t_sum_mu=0;
		if(_nb_flag) {
			for (int g=0; g<_G;g++){
				t_sum_z+=_z[g][n];
				t_sum_mu+=_mu[g][_kg[g][n]]*_s[g][n];
			}
		} else {
			for (int g=0; g<_G;g++){
				t_sum_z+=_z[g][n];
				t_sum_mu+=_mu[g][_kg[g][n]];
			}
		}
		sum_z+=t_sum_z;
		sum_mu+=t_sum_mu;
	}
	double temppn=gsl_ran_gamma(rgen, _prior_zetap+sum_z,1/(sum_mu+ _prior_zetap/_prior_mp));
	for (int n=0; n<_N;n++)
		_pn[n]=temppn;
}

/// mug
void datas::update_mu_nu(gsl_rng* rgen){
	///prior_b_mu= shape  and prior_mg_mu[g]=mean
	///***************************************************************************************///
	for (int g=0; g<_G;g++){
		for(int l=0;l<_K;l++){
			_sum_countk[l]=0; _sum_zk[l]=0; _sum_xk[l]=0;
		}
		for (int k=0; k<_K; k++) {
			_sum_xk[_map_mu_nu[g][k]]++;
		}
		for (int n=0; n<_N;n++){
			int q=_map_mu_nu[g][_kg[g][n]];
			_sum_zk[q]+=_z[g][n];
			if (_nb_flag) {
				_sum_countk[q]+=_pn[n]*_s[g][n];
			} else{
				_sum_countk[q]+=_pn[n];
			}
		}
		for(int l=0;l<_K;l++) {
			if (_sum_xk[l]>0) {
				_nu[g][l]=gsl_ran_gamma(rgen, _prior_b_mu+_sum_zk[l],
						1.0/(_sum_countk[l]+_prior_b_mu/_mg_mu[g]));
			} else {
				_nu[g][l]=-1;
			}
		}
		for(int l=0;l<_K;l++) {
			_mu[g][l]=_nu[g][_map_mu_nu[g][l]];
		}
	}
}

void datas::update_counts_an(gsl_rng* rgen){
	for (int n=0; n<_N;n++){
		for (int g=0; g<_G;g++){
			if (_z[g][n]==0) {
				_counts[g][n]=0;
			} else {
				double coverdisp=(_kappa_n_comp[0][n] + _kappa_n_comp[1][n]/_window[g] + _kappa_n_comp[2][n]/_z[g][n])*_kappa_g[g];
				double cmean=_window[g]*_an[n]*_z[g][n];
				_counts[g][n]=gsl_ran_negative_binomial(rgen, 1/(cmean*coverdisp+1), 1.0/coverdisp);
			}
		}
	}
}
double datas::compute_likelihood_marginalized_z_s() {
	double likelog=0;
	for (int n=0; n<_N;n++) {
		for (int g=0; g<_G;g++) {
			double tlike_g=0;
			double meanz=_pn[n] *_mu[g][_kg[g][n]];
			int minz=floor(meanz-5*sqrt(meanz));
			if (minz<0)
				minz=0;
			int maxz=ceil(meanz+5*sqrt(meanz));
			//			cout<<meanz<<","<<minz<<","<<maxz<<endl;;
			for (int z=minz; z<=maxz; z++) {
				tlike_g+=gsl_ran_negative_binomial_pdf(z, 1.0/(_pn[n]*_mu[g][_kg[g][n]]/_prior_zetas+1.0), _prior_zetas)*emission2(_counts[g][n],z, _an[n],_kappa_n_comp[0][n],_kappa_n_comp[1][n],_kappa_n_comp[3][n],_window[g],_kappa_n_comp[2][n],_kappa_g[g]);
			}
			likelog+=log(tlike_g);
		}
	}
	return likelog;
}
double* datas::compute_likelihood_marginalized_z_s_n() {
	double* likelog=new double [_N];
	for (int n=0; n<_N;n++) {
		likelog[n]=0;
		for (int g=0; g<_G;g++) {

			double meanz=_pn[n] *_mu[g][_kg[g][n]];
			int minz=floor(meanz-5*sqrt(meanz));
			if (minz<0) minz=0;
			if (_counts[g][n]>0) minz=1;
			int maxz=ceil(meanz+5*sqrt(meanz));

			for (int z=minz; z<=maxz; z++){
				likelog[n]+=log(gsl_ran_negative_binomial_pdf(z, 1.0/(_pn[n]*_mu[g][_kg[g][n]]/_prior_zetas+1.0), _prior_zetas)*\
						emission2(_counts[g][n],z, _an[n],_kappa_n_comp[0][n],_kappa_n_comp[1][n],_kappa_n_comp[3][n],_window[g],_kappa_n_comp[2][n],_kappa_g[g]));
				//				if (isinf(likelog[n])) {
				//					cout<<temp<<endl;
				//					cout<<likelog[n]<<endl;
				//					cout << "n=" << n << " g=" << g << " _counts[g][n] " << _counts[g][n] << " _z[g][n] " << _z[g][n] << endl;
				//					cout << "params _an[n] " << _an[n] << " _kappa_n_comp[0][n] " << _kappa_n_comp[0][n] << " _kappa_n_comp[1][n] " << _kappa_n_comp[1][n] << " _window[g] " << _window[g] << " _kappa_n_comp[2][n] " << _kappa_n_comp[2][n] << " _kappa_g[g] " << _kappa_g[g] << endl;
				//					cout<<_window[g]<<endl;
				//					double overdisp=(_kappa_n_comp[0][n]+_kappa_n_comp[1][n]/_window[g]+_kappa_n_comp[2][n]/_z[g][n])*_kappa_g[g];
				//					double mean=_window[g]*_an[n]*_z[g][n];
				//					double tprob=gsl_ran_negative_binomial_pdf(_counts[g][n], 1.0/(mean*overdisp+1.0), 1.0/overdisp);
				//					double tl2=emission2(_counts[g][n],z, _an[n],_kappa_n_comp[0][n],_kappa_n_comp[1][n],_kappa_n_comp[3][n],_window[g],_kappa_n_comp[2][n],_kappa_g[g],true);
				//					double tl=gsl_ran_negative_binomial_pdf(z, 1.0/(_pn[n]*_mu[g][_kg[g][n]]/_prior_zeta+1.0), _prior_zeta);
				//					cout << "overdisp " << overdisp << " mean " << mean << " trpob " <<  tprob << " log(tprob) " << log(tprob) << " tl2 " << tl2 << endl;
				//					cout << "pn " << _pn[n] << " zeta " <<_prior_zeta << " trpob " <<  tl <<   endl;
				//					exit(EXIT_FAILURE);
				//
				//					//					double overdisp=(kappa_i + kappa_w/w+kappa_z/z)*kappa_g;
				//					//					return gsl_ran_negative_binomial_pdf(counts, 1.0/(w*a*z*overdisp+1.0), 1.0/overdisp);
				//
				//				}
			}
		}
	}
	return likelog;
}
double* datas::compute_likelihood_n() {
	double* likelog=new double [_N];
	for (int n=0; n<_N;n++) {
		likelog[n]=0; // sub-summing to avoid potential numerical problems
		for (int g=0; g<_G;g++) {
			likelog[n]+=log(emission2(_counts[g][n],_z[g][n], _an[n],_kappa_n_comp[0][n],_kappa_n_comp[1][n],_kappa_n_comp[3][n],_window[g],_kappa_n_comp[2][n],_kappa_g[g]));
			likelog[n]+=log(gsl_ran_negative_binomial_pdf(_z[g][n], 1.0/(_pn[n]*_mu[g][_kg[g][n]]/_prior_zetas+1.0), _prior_zetas));

		}
		//		cout<<n<<": "<<likelog[n]<<endl;
	}
	return likelog;
}
double datas::compute_likelihood() {
	double likelog=0;
	for (int n=0; n<_N;n++) {
		double tlikelog=0; // sub-summing to avoid potential numerical problems
		for (int g=0; g<_G;g++) {
			double told=tlikelog;
			double tl=log(emission2(_counts[g][n],_z[g][n], _an[n],_kappa_n_comp[0][n],_kappa_n_comp[1][n],_kappa_n_comp[3][n],_window[g],_kappa_n_comp[2][n],_kappa_g[g]));
			tlikelog+=tl;
			if (isinf(tl) || (isinf(tlikelog) && !isinf(told))) {
				cout << "n=" << n << " g=" << g << " _counts[g][n] " << _counts[g][n] << " _z[g][n] " << _z[g][n] << endl;
				cout << "params _an[n] " << _an[n] << " _kappa_n_comp[0][n] " << _kappa_n_comp[0][n] << " _kappa_n_comp[1][n] " << _kappa_n_comp[1][n] << " _window[g] " << _window[g] << " _kappa_n_comp[2][n] " << _kappa_n_comp[2][n] << " _kappa_g[g] " << _kappa_g[g] << endl;
				cout<<_window[g]<<endl;
				double overdisp=(_kappa_n_comp[0][n]+_kappa_n_comp[1][n]/_window[g]+_kappa_n_comp[2][n]/_z[g][n])*_kappa_g[g];
				double mean=_window[g]*_an[n]*_z[g][n];
				double tprob=gsl_ran_negative_binomial_pdf(_counts[g][n], 1.0/(mean*overdisp+1.0), 1.0/overdisp);
				double tl2=log(emission2(_counts[g][n],_z[g][n], _an[n],_kappa_n_comp[0][n],_kappa_n_comp[1][n],_kappa_n_comp[3][n],_window[g],_kappa_n_comp[2][n],_kappa_g[g]));
				cout << "overdisp " << overdisp << " mean " << mean << " trpob " <<  tprob << " log(tprob) " << log(tprob)<< " tl " << tl << " tl2 " << tl2 << endl;
				cout << "tlikelog " << tlikelog << " told " << told << endl;
			}
		}
		likelog+=tlikelog;
		//cout << "n=" << n << " tlikelog " << tlikelog << endl;
	}
	return likelog;
}


double datas::compute_likelihood_marginalized_z() {
	double likelog=0;
	for (int n=0; n<_N;n++) {
		double tlikelog=0; // sub-summing to avoid potential numerical problems
		for (int g=0; g<_G;g++) {
			double tlike_g=0;
			double meanz=_pn[n]*_s[g][n]*_mu[g][_kg[g][n]];
			int minz=floor(meanz-5*sqrt(meanz));
			if (minz<0) {
				minz=0;
			}
			int maxz=ceil(meanz+5*sqrt(meanz));
			for (int z=minz; z<=maxz; z++) {
				tlike_g+=gsl_ran_poisson_pdf(z,meanz)*emission2(_counts[g][n],z, _an[n],_kappa_n_comp[0][n],_kappa_n_comp[1][n],_kappa_n_comp[3][n],_window[g],_kappa_n_comp[2][n],_kappa_g[g]);
			}
			tlikelog+=log(tlike_g);
		}
		likelog+=tlikelog;
	}
	return likelog;
}

double datas::compute_likelihood_marginalized_z_k(gsl_rng* rgen, bool predictiveposterior) {
	//	this->update_gamma_kg(rgen);  //warning the values of _gamma are meaningfull only when _coolingtemp==1
	double likelog=0;
	for (int n=0; n<_N;n++) {
		double tlikelog=0;
		for (int g=0; g<_G;g++) {
			_predictiveposterior_cdf[g][n]=0;
			int minz_glob=-1;
			int maxz_glob=-1;
			double predpost_cdf=0;
			for (int k=0; k<_K; k++) {
				double meanz=_pn[n]*_s[g][n]*_mu[g][k];
				int minz=floor(meanz-5*sqrt(meanz));
				if (minz<0) {
					minz=0;
				}
				int maxz=ceil(meanz+5*sqrt(meanz));
				if (minz_glob<0 || minz<minz_glob) {
					minz_glob=minz;
				}
				if (maxz_glob<0 || maxz>maxz_glob) {
					maxz_glob=maxz;
				}
			}
			double tlike_g=0;
			for (int z=minz_glob; z<=maxz_glob; z++) {
				double pz=0;
				for (int k=0; k<_K; k++) {
					pz+=_gamma[n][k]*gsl_ran_poisson_pdf(z,_pn[n]*_s[g][n]*_mu[g][k]);
				}
				tlike_g+=pz*emission2(_counts[g][n],z,_an[n],_kappa_n_comp[0][n],_kappa_n_comp[1][n],_kappa_n_comp[3][n],_window[g],_kappa_n_comp[2][n],_kappa_g[g]);
				if (predictiveposterior) {
					predpost_cdf += pz*emission2_cdf(_counts[g][n],z,_an[n],_kappa_n_comp[0][n],_kappa_n_comp[1][n],_kappa_n_comp[3][n],_window[g],_kappa_n_comp[2][n],_kappa_g[g]);
				}
			}
			if (predictiveposterior) {
				_predictiveposterior_cdf[g][n]=predpost_cdf-0.5*tlike_g; // -0.5*tlike_g is a trick to make it symmetrical in a context of discrete obs
			}
			tlikelog+=log(tlike_g);
		}
		likelog+=tlikelog;
	}
	return likelog;
}

void datas::set_flags(bool nbfl,bool zetafl,bool kappan_flag,bool kappa_i_flag,bool kappa_w_flag,bool kappa_w2_flag,bool kappa_z_flag, bool kappa_g_flag, bool window_adjustment_flag) {
	_nb_flag=nbfl;
	_zetak_flag=zetafl;
	_kappan_flag=kappan_flag;
	_kappa_n_comp_flag[0]=kappa_i_flag;
	_kappa_n_comp_flag[1]=kappa_w_flag;
	_kappa_n_comp_flag[2]=kappa_z_flag;
	_kappa_n_comp_flag[3]=kappa_w2_flag;
	_kappa_g_flag=kappa_g_flag;
	_window_adjustment_flag= window_adjustment_flag;
}

void datas::allocate_data_kg(){

	cout<<"Initialize the data structure."<<endl;

	//_prior_m_mu=0.2; //not used
	_prior_b_mu=10;
	_prior_d0=0.4;
	_prior_b_shape=1;     _prior_b_scale=1;
	_prior_mg_shape=5;
	_prior_zetap=1;_prior_zetaa=1;_prior_mp=1;_prior_ma=1;
	_prior_zetas=1;


	_kappa_n_comp = new double* [Ncomp];
	_kappa_n_comp_flag = new bool [Ncomp];
	for (int icomp=0; icomp<Ncomp; icomp++) {
		_kappa_n_comp[icomp] = new double [_N];
	}
	_naccept_kappa_n_comp_moves = new int [3];
	_nreject_kappa_n_comp_moves = new int [3];
	_prior_a_kappa_n_comp = new double [3];
	_prior_b_kappa_n_comp = new double [3];
	_prior_a_kappa_n_comp[0]=1;_prior_b_kappa_n_comp[0]=1;
	_prior_a_kappa_n_comp[1]=1;_prior_b_kappa_n_comp[1]=100;
	_prior_a_kappa_n_comp[2]=1;_prior_b_kappa_n_comp[2]=1;
	_prior_a_kappa_n_comp[3]=1;_prior_b_kappa_n_comp[3]=0.1;

	_counts= new int*[_G];
	for (int i=0; i<_G;i++)
		_counts[i]= new int[_N];

	_s= new double*[_G];
	for (int i=0; i<_G;i++)
		_s[i]= new double[_N];

	_predictiveposterior_cdf= new double*[_G];
	_predictiveposterior_pdf= new double*[_G];
	_predictiveposterior_tail= new double*[_G];
	for (int i=0; i<_G;i++) {
		_predictiveposterior_cdf[i]= new double[_N];
		_predictiveposterior_pdf[i]= new double[_N];
		_predictiveposterior_tail[i]= new double[_N];
	}

	_z= new int*[_G];
	for (int i=0; i<_G;i++)
		_z[i]= new int[_N];

	_dk= new int*[_G];
	for (int i=0; i<_G;i++)
		_dk[i]= new int[_K];

	_temp= new int*[_G];
	for (int i=0; i<_G;i++)
		_temp[i]= new int[_N];
	_temp1= new double[_K];
	_temp2= new double[_K];
	_mg_mu=new double [_G];

	_kg= new int *[_G];
	for (int i=0; i<_G;i++)
		_kg[i]=new int [_N];

	_gamma= new double*[_N];
	_count_kg= new int*[_N];
	for (int i=0; i<_N;i++) {
		_gamma[i]=new double [_K];
		_count_kg[i]=new int [_K];
	}
	_relweight_fcoolingtemp = new double [_G+1];

	_kappa_g= new double[_G];

	_ocur_pkn=new  double[_K];
	_prior_zetask= new double[_K];
	_mu= new double*[_G];
	for (int i=0; i<_G;i++)
		_mu[i]=new double [_K];
	_nu= new double*[_G];
	for (int i=0; i<_G;i++)
		_nu[i]=new double [_K];
	_map_mu_nu = new int* [_G];
	for (int i=0; i<_G;i++)
		_map_mu_nu[i]=new int [_K];

	_table_emis=new float* [100];
	for (int i=0; i<100;i++)
		_table_emis[i]=new float[100];

	_sum_xk=new int[_K];
	_sum_zk=new int[_K];
	_sum_countk=new double[_K];

	_sum_xn=new int[_N];
	_sum_zn=new int[_N];
	_sum_mun=new double[_N];

	_window= new int[_G];
	_pn= new double[_N];
	_an= new double[_N];
	_bg= new double[_G];
	_probas_k=new double[_K];
	_cells=new string[_N];
	_genes_frag=new string[_G];
	_genes=new string[_G];
	_clusters=new string [_K];
	char* tempc=new char[2];
	for(int l=0;l<_K;l++){
		snprintf(tempc, 2, "%d", l);
		_clusters[l]=string(("cluster_"+string(tempc)));
		cout<<_clusters[l]<<endl;
	}
	_fact_sum_zk=new float[_K];
	_d=new int[_G];
	cout<<"Data structures initialized."<<endl;
}

void datas::set_coolingtemp(double temp)
{
	if (temp>=0.1) { // numerical problem with pow giving +Inf when temp is too low
		_coolingtemp=temp;
		//cout << "_relweight_fcoolingtemp" << endl;
		for (int g=0; g<_G+1;g++) {
			_relweight_fcoolingtemp[g]=pow(g,1.0/_coolingtemp);
			//cout << " " << _relweight_fcoolingtemp[g];
		}
		//cout << endl;
	}
}

void datas::initialize_data_kg(gsl_rng* rgen){
	_mhstep=0.1;
	_epsilon=0.99;
	_a=1;
	apratio=1.00;

	_var_kappa_g=0.01;

	this->set_coolingtemp(1);

	_naccept_kappag_move=0;
	_nreject_kappag_move=0;
	for (int icomp=0; icomp<3; icomp++) {
		_naccept_kappa_n_comp_moves[icomp]=0;
		_nreject_kappa_n_comp_moves[icomp]=0;
	}
	_naccept_an_move=0;
	_nreject_an_move=0;
	_nchanges_kg_move=0;
	_nnochanges_kg_move=0;

	for (int g=0; g<_G;g++)
		_d[g]=1;

	for (int g=0; g<_G;g++)
		_bg[g]=1;

	for (int g=0; g<_G;g++)
		_mg_mu[g]=_prior_mg_mean; //(g+1)/200.0+1;
	for (int l=0; l<_K;l++)
		_prior_zetask[l]=100;
	for (int g=0; g<_G;g++)
		for (int l=0; l<_K;l++)
			_dk[g][l]= gsl_rng_uniform_int (rgen,2);

	for (int g=0; g<_G;g++) {
		for (int l=0; l<_K;l++) {
			_nu[g][l]=_prior_mg_mean;//gsl_ran_gamma(rgen,  prior_m_mu*prior_b_mu, 1.0/prior_b_mu);
		}
		for (int l=0; l<_K;l++) {
			_map_mu_nu[g][l]=l;
		}
		for (int l=0; l<_K;l++) {
			_mu[g][l]=_nu[g][_map_mu_nu[g][l]];
		}
	}

	for (int i=0; i<_N;i++) {
		for (int l=0; l<_K;l++) {
			_gamma[i][l]=1.0/_K; //not used when using _count_kg
		}
	}

	for (int n=0; n<_N;n++) {
		for (int l=0; l<_K;l++) {
			_count_kg[n][l]=0;
		}
		for (int g=0; g<_G;g++) {
			_kg[g][n]=gsl_rng_uniform_int (rgen,_K);
			_count_kg[n][_kg[g][n]]++;
		}
	}

	for (int n=0; n<_N;n++){

		int totcountsn=0;
		int totw=0;
		for (int g=0; g<_G; g++) {
			//			cout<<_counts[g][n]<<endl;
			totcountsn+=_counts[g][n];
			totw+=_window[g];
		}
		if (totcountsn==0) {
			totcountsn=1;
		}
		_pn[n]=pow((double)totcountsn/((double)totw*_prior_mg_mean ),0.5);
		//		cout << totcountsn << "/(" << totw << "*" << _prior_mg_mean << "*" << _pn[n] << ")";
		_an[n]=pow((double)totcountsn/((double)totw*_prior_mg_mean ),0.5);
		//		cout << " -> an[" << n << "]" << "=" << _an[n] << endl;
	}
	//	exit (EXIT_FAILURE);

	for (int g=0; g<_G;g++) {
		_kappa_g[g]=1;
	}

	for (int n=0; n<_N;n++){
		for (int icomp=0; icomp<Ncomp; icomp++) {
			if (_kappa_n_comp_flag[icomp]) {
				_kappa_n_comp[icomp][n]=_prior_a_kappa_n_comp[icomp]*_prior_b_kappa_n_comp[icomp]*5;
			} else {
				_kappa_n_comp[icomp][n]=0;
			}
		}
	}
	_prior_zetas=1;

	for (int g=0; g<_G;g++) {
		for (int n=0; n<_N;n++) {
			_s[g][n]=1;
			_predictiveposterior_cdf[g][n]=-1;
			_predictiveposterior_pdf[g][n]=-1;
			_predictiveposterior_tail[g][n]=0;
		}
	}

	for (int g=0; g<_G;g++) {
		for (int n=0; n<_N;n++) {
			//			_z[g][n]=(int)(_prior_mg_mean*_pn[n])+1;
			_z[g][n]=ceil((double)_counts[g][n]/((double)_window[g]*_an[n]))+1;
		}
	}
	cout<<"an_pn initialized: "<<endl;
	for (int nn=0; nn<_N;nn++) printf ("%f \t", _an[nn]);
	cout<<"\n"<<endl;

	for (int nn=0; nn<_N;nn++) printf ("%f \t", _pn[nn]);
	cout<<"\n"<<endl;

}

///**********************///

void datas::update_b(gsl_rng* rgen){
	double proposed=gsl_ran_gamma(rgen,1.0/pow(_mhstep,2.0),_prior_b_mu*pow(_mhstep,2.0));

	double temp=0, templog=0, templogmg=0, templogmgprop=0 ;
	int nrdata=0; //_K*_G accounting for the same k  for d=0
	for(int l=0;l<_K;l++){
		double t_temp=0, t_templog=0, t_templogmg=0, t_templogmgprop=0;
		int t_nrdata=0;
		for (int g=0; g<_G;g++){
			if (_nu[g][l]>0) { //=-1 when no cluster pointing to l
				t_nrdata++;
				t_templog+=log(_nu[g][l]);
				t_temp+=_nu[g][l]/_mg_mu[g];
				t_templogmg+=_prior_b_mu*log(_prior_b_mu/_mg_mu[g]);
				t_templogmgprop+=proposed*log(proposed/_mg_mu[g]);
			}
		}
		temp += t_temp; templog += t_templog; templogmg += t_templogmg; templogmgprop += t_templogmgprop;
		nrdata += t_nrdata;
	}

	// double proba_old = (_prior_b_mu-1)*templog - _prior_b_mu*temp + templogmg -
	//   gsl_sf_lngamma (_prior_b_mu)*nrdata+ log(gsl_ran_gamma_pdf(proposed,1.0/pow(_mhstep,2.0),_prior_b_mu*pow(_mhstep,2.0)))+
	//   (_prior_b_shape-1) * log(_prior_b_mu)-_prior_b_mu/_prior_b_scale;

	// double proba_prop= (proposed-1)*templog - proposed*temp + templogmgprop -
	//   gsl_sf_lngamma (proposed)*nrdata + log(gsl_ran_gamma_pdf(_prior_b_mu,1.0/pow(_mhstep,2.0),proposed*pow(_mhstep,2.0)))+
	//   (_prior_b_shape-1) * log(proposed)-proposed/_prior_b_scale;

	double logratio = 0;
	logratio += (proposed-_prior_b_mu)*templog;
	logratio += -(proposed-_prior_b_mu)*temp;
	logratio += templogmgprop-templogmg;
	logratio += -(gsl_sf_lngamma(proposed)-gsl_sf_lngamma(_prior_b_mu))*(double)nrdata;
	logratio += log(gsl_ran_gamma_pdf(_prior_b_mu,1.0/pow(_mhstep,2.0),proposed*pow(_mhstep,2.0)))-log(gsl_ran_gamma_pdf(proposed,1.0/pow(_mhstep,2.0),_prior_b_mu*pow(_mhstep,2.0)));
	logratio += (_prior_b_shape-1)*(log(proposed)-log(_prior_b_mu));
	logratio += -(proposed-_prior_b_mu)/_prior_b_scale;

	double r=gsl_ran_flat(rgen,0,1);
	// if ((log(r)<proba_prop-proba_old) != (log(r)<logratio)) {
	//   cout << "orig " << _prior_b_mu << " proposed " << proposed << " update_b init " << (proba_prop-proba_old) << " new " << logratio << endl;
	//   cout << "  proba_prop " << proba_prop << " proba_old " << proba_old << endl;
	//   cout << temp << " " << templog << " " << templogmg << " " << templogmgprop << " " << nrdata << endl;
	//   exit(1);
	// }
	//if(log(gsl_ran_flat(rgen,0,1))<proba_prop-proba_old) {
	//if(log(gsl_ran_flat(rgen,0,1))<logratio) {
	if(log(r)<logratio) {
		//if(log(r)<proba_prop-proba_old) {
		_prior_b_mu=proposed;
	}
}

void datas::deterministic_rescale_mg() {
	cout << "deterministic rescaling of mg" << endl;
	double tsum_mg_mu=0;
	for (int g=0; g<_G;g++){
		tsum_mg_mu+=_mg_mu[g];
	}
	double tmean_mg_mu=tsum_mg_mu/(double)_G;
	double talpha=_prior_mg_mean/tmean_mg_mu;
	cout << " mean of _mg_mu's is " << tmean_mg_mu << endl;
	cout << " applying rescaling factor " << talpha << endl;

	for (int g=0; g<_G;g++){
		_mg_mu[g]*=talpha;
	}
	for (int g=0; g<_G;g++){
		for(int l=0;l<_K;l++){
			if (_nu[g][l]>=0) {
				_nu[g][l]*=talpha;
			}
			_mu[g][l]*=talpha;
		}
	}
	for (int n=0; n<_N; n++) {
		_pn[n]/=talpha;
	}
	_prior_mp/=talpha; // correspond to mean of _pn[n]'s
	cout << "done" << endl;
}

void datas::update_rescale_mg_mu_pn(gsl_rng* rgen, bool pn_allequal) {
	double tmhstep=_mhstep*1;
	//cout << "entering update_rescale_mg_mu_pn" << endl;
	double proposed_alpha=gsl_ran_gamma(rgen,1.0/pow(tmhstep,2.0),pow(tmhstep,2.0));
	//cout << "proposed_alpha " << proposed_alpha << endl;
	double log_proposal_ratio=log(gsl_ran_gamma_pdf(1/proposed_alpha,1.0/pow(tmhstep,2.0),pow(tmhstep,2.0)))-log(gsl_ran_gamma_pdf(proposed_alpha,1.0/pow(tmhstep,2.0),pow(tmhstep,2.0)));

	//cout << "log_proposal_ratio " << log_proposal_ratio << endl;

	double tsum_mg_mu=0;
	double log_lik_ratio_mg=0;
	for (int g=0; g<_G;g++){
		tsum_mg_mu+=_mg_mu[g];
		log_lik_ratio_mg += log(gsl_ran_gamma_pdf(_mg_mu[g]*proposed_alpha,_prior_mg_shape,_prior_mg_mean/_prior_mg_shape));
		log_lik_ratio_mg -= log(gsl_ran_gamma_pdf(_mg_mu[g],_prior_mg_shape,_prior_mg_mean/_prior_mg_shape));
	}
	//cout << "mean of _mg_mu's is " << tsum_mg_mu/(double)_G << endl;
	//cout << "log_lik_ratio_mg " << log_lik_ratio_mg << endl;

	double log_lik_ratio_mu=0;
	int nrnu=0;
	for (int g=0; g<_G;g++){
		for(int l=0;l<_K;l++){
			if (_nu[g][l]>=0) {
				nrnu++;
				log_lik_ratio_mu += log(gsl_ran_gamma_pdf(_nu[g][l]*proposed_alpha,_prior_b_mu,(_mg_mu[g]/_prior_b_mu)*proposed_alpha));
				log_lik_ratio_mu -= log(gsl_ran_gamma_pdf(_nu[g][l],_prior_b_mu,_mg_mu[g]/_prior_b_mu));
			}
		}
	}
	//cout << "log_lik_ratio_mu " << log_lik_ratio_mu << endl;

	double log_lik_ratio_pn=0;
	if (!pn_allequal) {
		for (int n=0; n<_N; n++) {
			log_lik_ratio_pn += log(gsl_ran_gamma_pdf(_pn[n]/proposed_alpha,_prior_zetap,(_prior_mp/_prior_zetap)/proposed_alpha));
			log_lik_ratio_pn -= log(gsl_ran_gamma_pdf(_pn[n],_prior_zetap,(_prior_mp/_prior_zetap)));
		}
	} else {
		log_lik_ratio_pn += log(gsl_ran_gamma_pdf(_pn[0]/proposed_alpha,_prior_zetap,(_prior_mp/_prior_zetap)/proposed_alpha));
		log_lik_ratio_pn -= log(gsl_ran_gamma_pdf(_pn[0],_prior_zetap,(_prior_mp/_prior_zetap)));
	}

	//cout << "log_lik_ratio_pn " << log_lik_ratio_pn << endl;

	double log_lik_ratio_prior_mp=-(_prior_mp/proposed_alpha-_prior_mp);
	//cout << "log_lik_ratio_prior_mp " << log_lik_ratio_prior_mp << endl;

	double log_jacobian=0;
	if (!pn_allequal) {
		log_jacobian += (double)(_G+nrnu-(_N+1+2))*log(proposed_alpha); // +1 is for _prior_mp'<-_prior_mp/proposed_alpha +2 is for proposed_alpha'<-1/proposed_alpha
	} else {
		log_jacobian += (double)(_G+nrnu-(1+1+2))*log(proposed_alpha);
	}

	//cout << "log_jacobian " << log_jacobian << endl;

	double log_accept=log_proposal_ratio+log_lik_ratio_mg+log_lik_ratio_mu+log_lik_ratio_pn+log_lik_ratio_prior_mp+log_jacobian;
	//cout << "p_accept " << exp(log_accept) << endl;
	if(log(gsl_ran_flat(rgen,0,1))<log_accept) {
		//cout << " -> ************ accepted" << endl;
		for (int g=0; g<_G;g++){
			_mg_mu[g]*=proposed_alpha;
		}
		for (int g=0; g<_G;g++){
			for(int l=0;l<_K;l++){
				if (_nu[g][l]>=0) {
					_nu[g][l]*=proposed_alpha;
				}
				_mu[g][l]*=proposed_alpha;
			}
		}
		for (int n=0; n<_N; n++) {
			_pn[n]/=proposed_alpha;
		}
		_prior_mp/=proposed_alpha; // correspond to mean of _pn[n]'s
	}
}

void datas::update_mg(gsl_rng* rgen){
	for (int g=0; g<_G;g++){
		double temp=0;
		int nrdata=0;
		double proposed=gsl_ran_gamma(rgen,1.0/pow(_mhstep,2.0),_mg_mu[g]*pow(_mhstep,2.0));

		for(int l=0;l<_K;l++){
			if (_nu[g][l]>=0) {
				nrdata++;
				temp+=_nu[g][l];
			}
		}

		// double proba_old = -_prior_b_mu/_mg_mu[g]*temp -nrdata*_prior_b_mu*log(_mg_mu[g])+
		//   log(gsl_ran_gamma_pdf(proposed,1.0/pow(_mhstep,2.0),_mg_mu[g]*pow(_mhstep,2.0)))+
		//   (_prior_mg_shape-1) * log(_mg_mu[g])-_mg_mu[g]/(_prior_mg_mean/_prior_mg_shape);

		// double proba_prop= -_prior_b_mu/proposed*temp -nrdata*_prior_b_mu*log(proposed)+
		//   log(gsl_ran_gamma_pdf(_mg_mu[g],1.0/pow(_mhstep,2.0),proposed*pow(_mhstep,2.0)))+
		//   (_prior_mg_shape-1) * log(proposed)-proposed/(_prior_mg_mean/_prior_mg_shape);

		double logratio = 0;
		logratio += -_prior_b_mu*temp*(1.0/proposed-1.0/_mg_mu[g]);
		logratio += -(double)nrdata*_prior_b_mu*(log(proposed)-log(_mg_mu[g]));
		logratio += log(gsl_ran_gamma_pdf(_mg_mu[g],1.0/pow(_mhstep,2.0),proposed*pow(_mhstep,2.0)))-log(gsl_ran_gamma_pdf(proposed,1.0/pow(_mhstep,2.0),_mg_mu[g]*pow(_mhstep,2.0)));
		logratio += (_prior_mg_shape-1)*(log(proposed)-log(_mg_mu[g]));
		logratio += -(proposed-_mg_mu[g])/(_prior_mg_mean/_prior_mg_shape);

		//cout << "update_mg init " << (proba_prop-proba_old) << " new " << logratio << endl;

		if(log(gsl_ran_flat(rgen,0,1))<logratio)
			_mg_mu[g]=proposed;
	}
}

void datas::update_prior_mg_mean(gsl_rng* rgen) {
	double proposed=gsl_ran_gamma(rgen,1.0/pow(_mhstep,2.0),_prior_mg_mean*pow(_mhstep,2.0));
	double logratio=0;
	logratio += -proposed + _prior_mg_mean; // Exp mean=1 prior
	logratio += log(gsl_ran_gamma_pdf(_prior_mg_mean,1.0/pow(_mhstep,2.0),proposed*pow(_mhstep,2.0)))-log(gsl_ran_gamma_pdf(proposed,1.0/pow(_mhstep,2.0),_prior_mg_mean*pow(_mhstep,2.0))); // proposal
	for (int g=0; g<_G;g++){
		logratio += log(gsl_ran_gamma_pdf(_mg_mu[g],_prior_mg_shape,proposed/_prior_mg_shape));
		logratio -= log(gsl_ran_gamma_pdf(_mg_mu[g],_prior_mg_shape,_prior_mg_mean/_prior_mg_shape));
	}
	if(log(gsl_ran_flat(rgen,0,1))<logratio) {
		_prior_mg_mean=proposed;
	}
}

void datas::update_mg_shape(gsl_rng* rgen){
	double proposed=gsl_ran_gamma(rgen,1.0/pow(_mhstep,2.0),_prior_mg_shape*pow(_mhstep,2.0));

	double templog=0, tempsum=0;
	for (int g=0; g<_G;g++){
		templog+=log(_mg_mu[g]);
		tempsum+= _mg_mu[g];
	}

	// double proba_old = (_prior_mg_shape-1)*templog - (double)_G*gsl_sf_lngamma (_prior_mg_shape)-(double)_G*_prior_mg_shape *log(_prior_mg_mean/_prior_mg_shape)-
	//   _prior_mg_shape/_prior_mg_mean*tempsum+
	//   log(gsl_ran_gamma_pdf(proposed,1.0/pow(_mhstep,2.0),_prior_mg_shape*pow(_mhstep,2.0)))-_prior_mg_shape; //+log(gsl_ran_gamma_pdf(_prior_mg_shape,2.0,0.5));
	// //                (1.1-1) * log(prior_mg_shape)-prior_mg_shape/2;

	// double proba_prop= (proposed-1)*templog - (double)_G*gsl_sf_lngamma (proposed)-(double)_G*proposed *log(_prior_mg_mean/proposed)-
	//   proposed/_prior_mg_mean*tempsum+
	//   log(gsl_ran_gamma_pdf(_prior_mg_shape,1.0/pow(_mhstep,2.0),proposed*pow(_mhstep,2.0)))-proposed; //+log(gsl_ran_gamma_pdf(proposed,2.0,0.5));

	double logratio = 0;
	logratio += (proposed-_prior_mg_shape)*templog;
	logratio += -(double)_G*(gsl_sf_lngamma(proposed)-gsl_sf_lngamma(_prior_mg_shape));
	logratio += -(double)_G*(proposed*log(_prior_mg_mean/proposed)-_prior_mg_shape*log(_prior_mg_mean/_prior_mg_shape));
	logratio += -((proposed-_prior_mg_shape)/_prior_mg_mean)*tempsum;
	logratio += log(gsl_ran_gamma_pdf(_prior_mg_shape,1.0/pow(_mhstep,2.0),proposed*pow(_mhstep,2.0)))-log(gsl_ran_gamma_pdf(proposed,1.0/pow(_mhstep,2.0),_prior_mg_shape*pow(_mhstep,2.0)));
	logratio += -(proposed-_prior_mg_shape);

	//cout << "update_mg_shape init " << (proba_prop-proba_old) << " new " << logratio << endl;

	if(log(gsl_ran_flat(rgen,0,1))<logratio)
		_prior_mg_shape=proposed;
}

void datas::update_prior_mp(gsl_rng* rgen){
	double proposed=gsl_ran_gamma(rgen,1.0/pow(_mhstep,2.0),_prior_mp*pow(_mhstep,2.0));

	double tempsum=0;
	for (int n=0; n<_N;n++)
		tempsum+= _pn[n];

	// double proba_old = -_prior_zetap/_prior_mp*tempsum - (double)_N*_prior_zetap*log(_prior_mp)-_prior_mp+
	//   log(gsl_ran_gamma_pdf(proposed,1.0/pow(_mhstep,2.0),_prior_mp*pow(_mhstep,2.0)));

	// double proba_prop= -_prior_zetap/proposed*tempsum - (double)_N*_prior_zetap*log(proposed)-proposed+
	//   log(gsl_ran_gamma_pdf(_prior_mp,1.0/pow(_mhstep,2.0),proposed*pow(_mhstep,2.0)));

	double logratio = 0;
	logratio += -_prior_zetap*tempsum*(1.0/proposed-1.0/_prior_mp);
	logratio += -(double)_N*_prior_zetap*(log(proposed)-log(_prior_mp));
	logratio += -(proposed-_prior_mp);
	logratio += log(gsl_ran_gamma_pdf(_prior_mp,1.0/pow(_mhstep,2.0),proposed*pow(_mhstep,2.0)))-log(gsl_ran_gamma_pdf(proposed,1.0/pow(_mhstep,2.0),_prior_mp*pow(_mhstep,2.0)));

	//cout << "update_prior_mp init " << (proba_prop-proba_old) << " new " << logratio << endl;

	if(log(gsl_ran_flat(rgen,0,1))<logratio)
		_prior_mp=proposed;
}

void datas::update_prior_zetap(gsl_rng* rgen){
	double proposed=gsl_ran_gamma(rgen,1.0/pow(_mhstep,2.0),_prior_zetap*pow(_mhstep,2.0));

	double templog=0, tempsum=0;
	for (int n=0; n<_N;n++){
		templog+=log(_pn[n]);
		tempsum+= _pn[n];
	}

	double proba_old = (_prior_zetap-1)*templog - _prior_zetap/_prior_mp*tempsum - (double)_N*gsl_sf_lngamma (_prior_zetap)+(double)_N*_prior_zetap *log(_prior_zetap/_prior_mp)+
			log(gsl_ran_gamma_pdf(proposed,1.0/pow(_mhstep,2.0),_prior_zetap*pow(_mhstep,2.0)))-_prior_zetap;
	//                (1.1-1) * log(prior_mg_shape)-prior_mg_shape/2;

	double proba_prop = (proposed-1)*templog - proposed/_prior_mp*tempsum - (double)_N*gsl_sf_lngamma (proposed)+(double)_N*proposed *log(proposed/_prior_mp)+
			log(gsl_ran_gamma_pdf(_prior_zetap,1.0/pow(_mhstep,2.0),proposed*pow(_mhstep,2.0)))-proposed;

	if(log(gsl_ran_flat(rgen,0,1))<proba_prop-proba_old)
		_prior_zetap=proposed;
}

void datas::update_prior_zetaa(gsl_rng* rgen){
	double proposed=gsl_ran_gamma(rgen,1.0/pow(_mhstep,2.0),_prior_zetaa*pow(_mhstep,2.0));

	double templog=0, tempsum=0;
	for (int n=0; n<_N;n++){
		templog+=log(_an[n]);
		tempsum+=_an[n];
	}

	double proba_old = (_prior_zetaa-1)*templog - _prior_zetaa/_prior_ma*tempsum - (double)_N*gsl_sf_lngamma(_prior_zetaa)+(double)_N*_prior_zetaa *log(_prior_zetaa/_prior_ma)+
			log(gsl_ran_gamma_pdf(proposed,1.0/pow(_mhstep,2.0),_prior_zetaa*pow(_mhstep,2.0)))-_prior_zetaa/1;

	double proba_prop= (proposed-1)*templog - proposed/_prior_ma*tempsum - (double)_N*gsl_sf_lngamma (proposed)+(double)_N*proposed *log(proposed/_prior_ma)+
			log(gsl_ran_gamma_pdf(_prior_zetaa,1.0/pow(_mhstep,2.0),proposed*pow(_mhstep,2.0)))-proposed/1;

	if(log(gsl_ran_flat(rgen,0,1))<proba_prop-proba_old)
		_prior_zetaa=proposed;
}

void datas::update_prior_ma(gsl_rng* rgen){
	double tempsum=0;

	double proposed=gsl_ran_gamma(rgen,1.0/pow(_mhstep,2.0),_prior_ma*pow(_mhstep,2.0));
	double proba_prop=0, proba_old=0;

	for (int n=0; n<_N;n++)
		tempsum+=_an[n];

	proba_old = -_prior_zetaa/_prior_ma*tempsum - (double)_N*_prior_zetaa*log(_prior_ma)-_prior_ma+
			log(gsl_ran_gamma_pdf(proposed,1.0/pow(_mhstep,2.0),_prior_ma*pow(_mhstep,2.0)));

	proba_prop= -_prior_zetaa/proposed*tempsum - (double)_N*_prior_zetaa*log(proposed)-proposed+
			log(gsl_ran_gamma_pdf(_prior_ma,1.0/pow(_mhstep,2.0),proposed*pow(_mhstep,2.0)));

	if(log(gsl_ran_flat(rgen,0,1))<proba_prop-proba_old)
		_prior_ma=proposed;
}

void datas::update_s(gsl_rng* rgen){
	if (!_nb_flag) {
		return;
	}
	if(_zetak_flag)
		for (int n=0; n<_N;n++)
			for (int g=0; g<_G;g++)
				_s[g][n]=gsl_ran_gamma(rgen, _prior_zetask[_kg[g][n]]+_z[g][n],1/(_mu[g][_kg[g][n]]*_pn[n]+ _prior_zetask[_kg[g][n]]));
	else
		for (int n=0; n<_N;n++)
			for (int g=0; g<_G;g++)
				_s[g][n]=gsl_ran_gamma(rgen, _prior_zetas+_z[g][n],1/(_mu[g][_kg[g][n]]*_pn[n]+ _prior_zetas));
}

void datas::update_prior_zetas(gsl_rng* rgen){

	double templog=0, tempsum=0;
	double proposed=gsl_ran_gamma(rgen,1.0/pow(_mhstep,2.0),_prior_zetas*pow(_mhstep,2.0));
	double proba_prop=0, proba_old=0;

	for (int g=0; g<_G;g++){
		for (int n=0; n<_N;n++){
			templog+=log(_s[g][n]);
			tempsum+= _s[g][n];
		}
	}

	proba_old = (_prior_zetas-1)*templog - _prior_zetas *tempsum - (double)(_N*_G)*gsl_sf_lngamma (_prior_zetas)+(double)(_N*_G)*_prior_zetas *log(_prior_zetas)+
			log(gsl_ran_gamma_pdf(proposed,1.0/pow(_mhstep,2.0),_prior_zetas*pow(_mhstep,2.0)))+(1-1)*log(_prior_zetas)-_prior_zetas/1.0;
	proba_prop= (proposed-1)*templog - proposed *tempsum - (double)(_N*_G)*gsl_sf_lngamma (proposed)+(double)(_N*_G)*proposed *log(proposed)+
			log(gsl_ran_gamma_pdf(_prior_zetas,1.0/pow(_mhstep,2.0),proposed*pow(_mhstep,2.0)))+(1-1)*log(proposed)-proposed/1.0;

	//PN: pb here - not an exp(1) as prior as written in the text? here a gamma shape=3 scale=3 (mean 1). Of note I see no reason for centering the prior around 1...
	//		cout<<"prior_zetas: "<<_prior_zetas<<" ";
	//		cout<<proposed<<endl;
	//		if(log(gsl_ran_flat(rgen,0,1))<proba_prop-proba_old && proposed>=16)
	if(log(gsl_ran_flat(rgen,0,1))<proba_prop-proba_old && proposed>=0.1)
		_prior_zetas=proposed;

	for (int l=0;l<_K;l++)
		_prior_zetask[l]=_prior_zetas;
}



void datas::run_iteration(gsl_rng* rgen, int iter, int iter_allnequal,int iter_update_kappa,string updates[],int n_updates,bool update_also_an){

	string update;
	for (int i=0;i<n_updates;i++){
		//	for(const auto& update : updates) {   // Range-for!
		update=updates[i];
		//		cout<<update<<endl;


		if (update.compare("update_b")==0) update_b(rgen);
		if (update.compare("update_mg")==0) update_mg(rgen);
		if (update.compare("update_mg_shape")==0) update_mg_shape(rgen);
		if (update.compare("update_rescale_mg_mu_pn")==0) update_rescale_mg_mu_pn(rgen, iter<iter_allnequal);
		if (update.compare("update_pn")==0){
			if(iter>=iter_allnequal)	update_pn(rgen);
			else 						update_pn_allequal(rgen);
		}
		if (update.compare("update_prior_mp")==0) update_prior_mp(rgen);
		if (update.compare("update_prior_zetap")==0)update_prior_zetap(rgen);

		if (update.compare("update_z")==0)
			if (iter%3==0)
				update_z(rgen);
		if (update.compare("update_prior_ma")==0)update_prior_ma(rgen);
		if (update.compare("update_prior_zetaa")==0)update_prior_zetaa(rgen);
		if (update.compare("update_mu_nu")==0) update_mu_nu(rgen);


		if (update.compare("update_kappa_an")==0) {
			if(iter%iter_update_kappa==0) {
				if (iter<iter_allnequal)
					update_kappa_n_comp_an(rgen, (iter<iter_allnequal),true, _G*0.3, 10,update_also_an);
				else if (iter<3000) {
					update_kappa_n_comp_an(rgen, (iter<iter_allnequal), false, _G*0.3, 10,update_also_an);
				} else {
					if (iter%3==0)
						update_kappa_n_comp_an(rgen, (iter<iter_allnequal), false, _G, 10,update_also_an);
				}
				update_kappa_g(rgen);
				update_var_kappa_g(rgen);
			}
//			if (iter<2000)
//				for (int nn=0; nn<_N;nn++)
//					_kappa_n_comp[2][nn]=1;
		}

		if (update.compare("update_s")==0)
			if (iter%3==0)
				update_s(rgen);
		if (update.compare("update_prior_zetas")==0)
			if (iter>=200)
				update_prior_zetas(rgen);
	}

}
