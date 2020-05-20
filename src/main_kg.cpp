# include <iostream>
# include <typeinfo>
# include <fstream>
# include <string>
# include <stdio.h>
# include <sys/stat.h>
# include"datas.h"
# include"additional_functions.h"


using namespace std;
/*???? version 2*/
int main(int argc, char *argv[]){

	bool booltest = false;
	//
	//	string param_file="/Users/mirauta/Google Drive/workspace_gg/SC/data/param_nozeta_crp_kappag.txt";
	//	//	//
	//	//	//	//	string data_file ="/Users/mirauta/Projects/SC/Single_cell/data/test/syntetic_counts_matrix.txt";
	//	//	//	//	string dest_folder="/Users/mirauta/Projects/SC/Single_cell/data/test/";
	//	//	//
	//	//	//
	//	//	//	//	string data_file ="/Users/mirauta/Projects/SC/Single_cell/data/MARSseq/MARSseq_counts_matrix.txt";
	//	//	//	//	string dest_folder="/Users/mirauta/Projects/SC/Single_cell/data/MARSseq/";
	//	//	//
	//	//	//	string data_file ="/Users/mirauta/Google Drive/workspace_gg/SC/data/simpn4an02/simpn4an02_counts_matrix.txt";
	//	//	//	string dest_folder="/Users/mirauta/Google Drive/workspace_gg/SC/data/simpn4an02/";
	//	//	//
	//	//	//	string data_file ="/Users/mirauta/Google Drive/workspace_gg/SC/data/simpn02an6/simpn02an6_counts_matrix.txt";
	//	//	//	string dest_folder="/Users/mirauta/Google Drive/workspace_gg/SC/data/simpn02an6/";
	//	//
	//	//	//	string data_file ="/Users/mirauta/Projects/SC/Single_cell/data/sim/sim_counts_matrix.txt";
	//	//	//	string dest_folder="/Users/mirauta/Projects/SC/Single_cell/data/sim/";
	//	//	//
	//	//	//	string data_file ="/Users/mirauta/Projects/SC/Single_cell/data/10x/10x_counts_matrix.txt";
	//	//	//	string dest_folder="/Users/mirauta/Projects/SC/Single_cell/data/10x/";
	//	string data_file ="/Users/mirauta/Google Drive/workspace_gg/SC/data/SMARTseq2_1000/SMARTseq2_1000_counts_matrix.txt";
	//	string dest_folder="/Users/mirauta/Google Drive/workspace_gg/SC/data/SMARTseq2_1000/";
	//	if (booltest) {	}
	//	else {
	cout<<"Reading arguments"<<argc<<endl;
	string data_file="null";
	string param_file="null";
	string dest_folder="null";
	//	}

	string run_type;
	string detail;
	int K=1;
	int verbose=0;
	int rseed=-1;
	//int K=-1;

	int num_arg = 1;
	while (num_arg < argc)
	{
		if ( (strcmp (argv[num_arg], "-destfolder") == 0)
				&& (num_arg + 1 < argc) ) {
			num_arg++;
			dest_folder=string(argv[num_arg]);
			cout << "dest_folder: " << dest_folder << endl;
			num_arg++;
		} else if ( (strcmp (argv[num_arg], "-datafile") == 0)
				&& (num_arg + 1 < argc) ) {
			num_arg++;
			data_file=string(argv[num_arg]);
			cout << "data_file: " << data_file << endl;
			num_arg++;
		} else if ( (strcmp (argv[num_arg], "-paramfile") == 0)
				&& (num_arg + 1 < argc) ) {
			num_arg++;
			param_file=string(argv[num_arg]);
			cout << "param_file: " << param_file << endl;
			num_arg++;
		} else if ( (strcmp (argv[num_arg], "-rseed") == 0)
				&& (num_arg + 1 < argc) ) {
			num_arg++;
			rseed=atoi(argv[num_arg]);
			cout << "rseed: " << rseed << endl;
			num_arg++;
			if (rseed<0) {
				cerr << "random seed (rseed) needs to be a positive integer" << endl;
				return(1);
			}
		} else if ( (strcmp (argv[num_arg], "-K") == 0)
				&& (num_arg + 1 < argc) ) {
			num_arg++;
			K=atoi(argv[num_arg]);
			cout << "K: " << K << endl;
			num_arg++;
		} else if ( (strcmp (argv[num_arg], "-verbose") == 0)
				&& (num_arg + 1 < argc) ) {
			num_arg++;
			verbose=atoi(argv[num_arg]);
			cout << "verbose: " << K << endl;
			num_arg++;
		} else {
			cerr << " unknown keyword: '" << argv[num_arg] << "' or missing the associated value" << endl;
			return(1);
		}
	}

	if (dest_folder.compare("null")==0 || data_file.compare("null")==0 || param_file.compare("null")==0) {
		cerr << "Usage : " << argv[0] << endl
				<< "   -datafile <filename> -paramfile <filename> -destfolder <foldername>" << endl
				<< "   [-rseed <positive integer> -K <#ofclusters>]" << endl
				<< endl ;
		return (1) ;
	}

	if (verbose) cout<<"verbose\n"<<endl;
	int  iter_out=100, itermax=100000000,
			nb_flag=0, zetak_flag=0, kappan_flag=0,
			kappa_i_flag=1, kappa_w_flag=0 ,kappa_w2_flag=0 , kappa_z_flag=0,
			all_genes_flag=0, iter_allnequal=0,
			iter_start_cooling=-1, iter_cooling_step=10, kappa_g_flag=0, iter_update_kappa=1, predictiveposterior_flag=1,window_adjustment_flag=1;
	double amplitude_cooling_step=0.001;

	read_parameters(param_file, &run_type, &detail, &K,&iter_out, &itermax, &nb_flag, &zetak_flag,&kappa_i_flag,
			&kappa_w_flag,&kappa_w2_flag,&kappa_z_flag,&kappan_flag,&all_genes_flag,&iter_start_cooling,&iter_cooling_step,&amplitude_cooling_step,&kappa_g_flag,&iter_update_kappa,&window_adjustment_flag);

	if (iter_start_cooling<0) {
		iter_start_cooling=itermax+1;
	}

	if (K<=0) {
		cerr << "invalid number of clusters (K): " << K << endl
				<< "K needs to be specified in the command line args or in the parameter files," << endl
				<< "and should be a strictly positive integer" << endl;
		return(1);
	}
	//	num_arg++;

	dest_folder+='/';
	cout<<"Running:"<<run_type<<" with data: "<< data_file <<" to: " << dest_folder<<endl;

	if (rseed<0) {
		rseed=time(0);
	}
	cout << "random generator seed is " << rseed << endl;
	gsl_rng* rgen;     rgen=gsl_rng_alloc(gsl_rng_rand48); gsl_rng_set(rgen,rseed);
	gsl_rng* rgen1=gsl_rng_alloc(gsl_rng_rand48); gsl_rng_set(rgen1,rseed+1);

	ofstream FILE;
	struct stat st;
	if(stat((dest_folder).c_str(),&st) != 0) {
		int retsys=system(("mkdir "+dest_folder).c_str());
		if (retsys!=0) {
			cerr << "Unable to create folder \"" << dest_folder << "\"" << endl;
			return(1);
		}
	}

	datas dat=datas(); /// allocate an empty datas object

	dat.read_metadata(data_file);
	dat.set_metadata(dat._G,dat._N,K);
	dat.allocate_data_kg();
	dat.set_flags((bool)nb_flag, (bool)zetak_flag, (bool)kappan_flag, (bool)kappa_i_flag, (bool)kappa_w_flag,(bool)kappa_w2_flag, (bool)kappa_z_flag, (bool)kappa_g_flag,(bool)window_adjustment_flag);

	dat.read_data(data_file); /// read counts from filexi

	cout<<"Running data: "<<dat._G<<" "<<dat._N<<" "<<dat._K<<endl;

	dat.initialize_data_kg(rgen); /// set initial values to kappa, a, pn, mu, k, z, d
	//	for (int nn=0; nn<dat._N;nn++) printf ("%f \t", dat._an[nn]);

	double* pn_0=new double[dat._N];

	for (int nn=0; nn<dat._N;nn++){
		pn_0[nn]=0;
		for (int g=0; g<dat._G;g++){
//			cout<<dat._counts[g][nn]<<endl;
			if (dat._counts[g][nn]>0)
				pn_0[nn]++;
		}
	}
	for (int nn=0; nn<dat._N;nn++)	cout<<pn_0[nn]<<' ';
	cout<<endl;
	double temp=0;
	for (int nn=0; nn<dat._N;nn++) temp=temp+pn_0[nn]/dat._N;
	for (int nn=0; nn<dat._N;nn++)	pn_0[nn]=pn_0[nn]/double(temp);
	for (int nn=0; nn<dat._N;nn++)	cout<<pn_0[nn]<<' ';

//	exit(EXIT_FAILURE);

	cout<<"Starting iterations"<<endl;
	cout<<"test"<<endl;
	for(int iter=0; iter<itermax;iter++){

		//	if (iter>=iter_start_cooling && (iter-iter_start_cooling)%iter_cooling_step==0) {
		//dat.set_coolingtemp(dat._coolingtemp*(1.0-amplitude_cooling_step));
		//}
		string updates []                          ={"update_b","update_mu_nu","update_mg","update_mg_shape","update_rescale_mg_mu_pn" ,"NOT_update_pn","update_prior_mp","update_prior_zetap",\
				"update_z","update_prior_ma","update_prior_zetaa","update_kappa_an","update_s","update_prior_zetas"};
		bool update_also_an=true;
		string updates0[]={"update_b","update_mu_nu","update_mg","update_mg_shape", "update_rescale_mg_mu_pn", "update_z","update_prior_ma","update_s"};

		//		if (iter<10000000)
		//			for (int n=0; n<dat._N;n++)
		//				for (int icomp=0; icomp<dat.Ncomp; icomp++)
		//					 dat._kappa_n_comp[icomp][n]=0.000001;
		if (iter<500){
			for (int n=0; n<dat._N;n++)
				for (int icomp=0; icomp<dat.Ncomp; icomp++)
					if (dat._kappa_n_comp[icomp][n]>2.0)
						dat._kappa_n_comp[icomp][n]=2.0;
			double temp=0.0;
			for (int nn=0; nn<dat._N;nn++)
				temp+=dat._an[nn]/dat._N;
			//			if (iter>1000)
			//				for (int nn=0; nn<dat._N;nn++)
			//					dat._an[nn]=dat._an[nn]/temp*15.0;
		}
		if (iter <100)
			for (int nn=0; nn<dat._N;nn++)
				dat._pn[nn]=pn_0[nn];
		if (iter==2000)
			updates[5]="update_pn";

		//		for (int g=0; g<dat._G;g++)
		//			dat._kappa_g[g]=4/(iter/20.0+1.0);
		//		if (iter<200)
		//			dat._prior_zetaa=10;
		//
		//		//
		//		if (iter<1000)
		//			dat._prior_ma=1;
		//		if (iter<100){
		//			//			dat._prior_zetap=100;
		//			double temp=0;
		//			for (int nn=0; nn<dat._N;nn++)
		//				temp+=dat._pn[nn];
		//			for (int nn=0; nn<dat._N;nn++)
		//				dat._pn[nn]=dat._pn[nn]/temp*dat._N;
		//			temp=0;
		//			for (int nn=0; nn<dat._N;nn++)
		//				temp+=dat._an[nn];
		//			for (int nn=0; nn<dat._N;nn++)
		//				dat._an[nn]=dat._an[nn]/temp*dat._N;
		//
		//			dat._prior_ma=1;
		//			dat._prior_mp=1;
		//			//			dat._prior_zetaa=100;
		//			//			update_also_an=false;
		//		}
		//string  updates[]={"update_b","update_mu_nu","update_mg","update_mg_shape","update_rescale_mg_mu_pn","update_prior_mp",	"update_z","update_kappa_an","update_s","update_prior_zetas"};

		//		for (int nn=0; nn<dat._N;nn++)  dat._kappa_n_comp[0][nn]=1.0;
		//		dat._prior_zetas=5;
		//
		//		//		dat._prior_ma=0.2;
		//		//		dat._prior_mp=4;
		//		bool update_also_an=true;
		//		string  updates[]={"update_b","update_mg","update_mg_shape","update_rescale_mg_mu_pn","update_pn","update_prior_zetap",\
		//				"update_z","update_prior_zetaa","update_s","update_kappa_an"};
		//
		//
		//		for (int nn=0; nn<dat._G;nn++)  dat._mu[nn][0]=2;
		//		for (int nn=0; nn<dat._N;nn++) dat._an[nn]=.2;
		//		for (int nn=0; nn<dat._N;nn++) dat._pn[nn]=5;
		//		cout<<"Marginal likelihood: "<<dat.compute_likelihood_marginalized_z_s()<<endl;
		//
		//		for (int nn=0; nn<dat._G;nn++)  dat._mu[nn][0]=2;
		//		for (int nn=0; nn<dat._N;nn++) dat._an[nn]=1.0;
		//		for (int nn=0; nn<dat._N;nn++) dat._pn[nn]=1.0;
		//		cout<<"Marginal likelihood: "<<dat.compute_likelihood_marginalized_z_s()<<endl;

		//		exit(EXIT_FAILURE);

		if (iter<100)
			dat.run_iteration(rgen,iter,iter_allnequal,iter_update_kappa,updates0,sizeof(updates)/sizeof(updates[0]),false);
		else
			dat.run_iteration(rgen,iter,iter_allnequal,iter_update_kappa,updates,sizeof(updates)/sizeof(updates[0]),update_also_an);

		if (iter%5000==0){
			cout<<"Likelihood: "<<dat.compute_likelihood()<<endl;
			cout<<"Marginal likelihood: "<<dat.compute_likelihood_marginalized_z_s()<<endl;
		}


		//		for (int k=0; k<1; k++)
		//			for (int nn=0; nn<11;nn++){
		//
		//				cout<<dat._map_mu_nu[nn][k]<<endl;
		//				cout<<dat._mu[nn][k]<<endl;
		//			}
		//				exit(EXIT_FAILURE);
		//		if ((iter%100==99)&(iter>500)){
		//
		////			exit(EXIT_FAILURE);
		//			dat.update_an_pn_ratio_n(rgen, iter_allnequal, iter_update_kappa);
		//			double*xx=dat.compute_likelihood_marginalized_z_s_n();
		//			cout<<"Marginal likelihood: "<<xx[0]<<endl;
		//
		//
		//
		//		}

		//		dat2._an[0]=121;
		//		cout<<dat2._an[0]<<endl;
		//		cout<<dat._an[0]<<endl;
		//		exit(EXIT_FAILURE);

		if(iter%iter_out==0){
			cout<<"at iteration: "<<iter<<endl;
			for (int nn=0; nn<dat._N;nn++) printf ("%f \t", dat._an[nn]);
			cout<<endl<<"zetaa: "<<dat._prior_zetaa<<"ma: "<<dat._prior_ma<< "\n"<<endl;

			for (int nn=0; nn<dat._N;nn++) printf ("%f \t", dat._pn[nn]);
			//			for (int nn=0; nn<11;nn++) printf ("%f \t", dat._mu[nn][0]);
			cout<<endl<<"zetap: "<<dat._prior_zetap<<"mp: "<<dat._prior_mp<<"\n"<<endl;
		}
		/**********************************************************************************************************************************/
		/**********************************************************************************************************************************/
		if(iter%iter_out==0){

			cout<<"At iter "<<iter<<endl;
			for (int icomp=0; icomp<4; icomp++) {
				cout << "_naccept_kappa_n_comp_moves[" << icomp << "] " << dat._naccept_kappa_n_comp_moves[icomp] << "/" << (double)dat._naccept_kappa_n_comp_moves[icomp] + dat._nreject_kappa_n_comp_moves[icomp] << endl;
				dat._naccept_kappa_n_comp_moves[icomp]=0;
				dat._nreject_kappa_n_comp_moves[icomp]=0;
			}
			cout << "_naccept_an_move " << dat._naccept_an_move << "/" << (double)dat._naccept_an_move + dat._nreject_an_move << endl;
			dat._naccept_an_move=0;
			dat._nreject_an_move=0;
			cout << "_naccept_kappag_move " << dat._naccept_kappag_move << "/" << (double)dat._naccept_kappag_move + dat._nreject_kappag_move << endl;
			dat._naccept_kappag_move=0;
			dat._nreject_kappag_move=0;


			if(iter==0)
				write_vec(dat._cells,dat._N, dest_folder,"pn_update",FILE,1);
			write_vec(dat._pn,dat._N, dest_folder,"pn_update",FILE,0);

			if(iter==0)
				write_vec(dat._cells,dat._N, dest_folder,"an_update",FILE,1);
			write_vec(dat._an,dat._N, dest_folder,"an_update",FILE,0);

			if(iter==0) {
				write_vec(dat._cells,dat._N,dest_folder,"kappa_i_update",FILE,1);
				write_vec(dat._cells,dat._N,dest_folder,"kappa_w_update",FILE,1);
				write_vec(dat._cells,dat._N,dest_folder,"kappa_w2_update",FILE,1);
				write_vec(dat._cells,dat._N,dest_folder,"kappa_z_update",FILE,1);
			}
			write_vec(dat._kappa_n_comp[0],dat._N, dest_folder,"kappa_i_update",FILE,0);
			write_vec(dat._kappa_n_comp[3],dat._N, dest_folder,"kappa_w2_update",FILE,0);
			write_vec(dat._kappa_n_comp[1],dat._N, dest_folder,"kappa_w_update",FILE,0);
			write_vec(dat._kappa_n_comp[2],dat._N, dest_folder,"kappa_z_update",FILE,0);

			write_vec(dat._kappa_g,dat._G, dest_folder,"kappa_g_update",FILE,iter==0);

			if(iter==0)
				write_vec(dat._clusters,dat._K, dest_folder,"zetask_update",FILE,1);
			write_vec(dat._prior_zetask,dat._K, dest_folder,"zetask_update",FILE,0);

			//write_matrix_2(dat._map_mu_nu,  dat._genes,dat._genes_frag, dat._clusters, iter, dat._G, dat._K,  dest_folder, "mu_update",  FILE, iter==0);
			write_matrix_2(dat._mu,  dat._genes,dat._genes_frag, dat._clusters, iter, dat._G, dat._K,  dest_folder, "mu_update",  FILE, iter==0);
			write_matrix_2(dat._kg,  dat._genes,dat._genes_frag, dat._cells, iter, dat._G, dat._N,  dest_folder, "kg_update",  FILE, iter==0);
			write_matrix_2(dat._count_kg,  dat._cells,dat._cells, dat._clusters, iter, dat._N, dat._K,  dest_folder, "nk_update",  FILE, iter==0);
			write_matrix_2(dat._gamma,  dat._cells,dat._cells, dat._clusters, iter, dat._N, dat._K,  dest_folder, "gammak_update",  FILE, iter==0);
			write_matrix_2(dat._z, dat._genes, dat._genes_frag, dat._cells, iter, dat._G, dat._N,  dest_folder, "z_update",  FILE, iter==0);
			write_matrix_2(dat._s, dat._genes, dat._genes_frag, dat._cells, iter, dat._G, dat._N,  dest_folder, "s_update",  FILE, iter==0);
			write_one_record_per_gene(dat._mg_mu,  dat._genes,dat._genes_frag, iter, dat._G, dest_folder, "mg_update",  FILE, iter==0);

			//            write_value((float)dat._epsilon,dest_folder,"epsilon_update",FILE,iter==0);
			write_value((float)dat._prior_b_mu,dest_folder,"prior_b_update",FILE,iter==0);
			write_value((float)dat._prior_mp,dest_folder,"prior_mp_update",FILE,iter==0);
			write_value((float)dat._prior_ma,dest_folder,"prior_ma_update",FILE,iter==0);
			write_value((float)dat._prior_zetap,dest_folder,"prior_zetap_update",FILE,iter==0);
			write_value((float)dat._prior_zetaa,dest_folder,"prior_zetaa_update",FILE,iter==0);
			write_value((float)dat._prior_mg_shape,dest_folder,"prior_mg_shape_update",FILE,iter==0);
			write_value((float)dat._prior_mg_mean,dest_folder,"prior_mg_mean_update",FILE,iter==0);
			write_value((float)dat._prior_zetas,dest_folder,"prior_zetas_update",FILE,iter==0);
			write_value((float)dat._var_kappa_g,dest_folder,"var_kappa_g",FILE,iter==0);

			write_value((float)dat._prior_d0,dest_folder,"prior_d0_update",FILE,iter==0);

			write_value(dat._coolingtemp,dest_folder,"coolingtemp",FILE,iter==0);

			write_value(dat.compute_likelihood(),dest_folder,"total_likelihood",FILE,iter==0);
			//write_value(dat.compute_likelihood_marginalized_z(),dest_folder,"total_likelihood_marginalized_z",FILE,iter==0);
			//write_value(dat.compute_likelihood_marginalized_z_k(rgen1,predictiveposterior_flag==1),dest_folder,"total_likelihood_marginalized_z_k",FILE,iter==0);
			if (predictiveposterior_flag==1) {
				write_matrix_2(dat._predictiveposterior_cdf, dat._genes, dat._genes_frag, dat._cells, iter, dat._G, dat._N,  dest_folder, "predictiveposterior_cdf",  FILE, iter==0);
				//write_matrix_2(dat._predictiveposterior_pdf, dat._genes, dat._genes_frag, dat._cells, iter, dat._G, dat._N,  dest_folder, "predictiveposterior_pdf",  FILE, iter==0);
				//write_matrix_2(dat._predictiveposterior_tail, dat._genes, dat._genes_frag, dat._cells, iter, dat._G, dat._N,  dest_folder, "predictiveposterior_tail",  FILE, iter==0);
			}
			if(iter==0){
				write_vec(dat._genes,dat._G, dest_folder,"d_update",FILE,1);
				write_vec(dat._genes_frag,dat._G, dest_folder,"d_update",FILE,0);
			}
			write_vec(dat._d,dat._G, dest_folder,"d_update",FILE,0);
		}
	}
	return 0;
}
