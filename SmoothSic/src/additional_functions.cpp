
# include "additional_functions.h"
# include "datas.h"

using namespace std;

void write_vec (int* x, int size, string dest_folder, string dest_file,   ofstream & file1, bool restart){

	if(restart) file1.open((dest_folder+dest_file).c_str());
	else file1.open((dest_folder+dest_file).c_str(),ios_base::app);
	if(!file1) cerr<<"NOT able to write/append to "<<(dest_folder+dest_file).c_str()<<endl;
	for(int i=0; i<size-1; i++){file1<<x[i]<<"\t";}
	file1<<x[size-1]<<"\n";
	file1.close();
}

void write_vec (float* x, int size, string dest_folder, string dest_file,   ofstream & file1, bool restart){

	if(restart) file1.open((dest_folder+dest_file).c_str());
	else file1.open((dest_folder+dest_file).c_str(),ios_base::app);
	if(!file1) cerr<<"NOT able to write/append to "<<(dest_folder+dest_file).c_str()<<endl;
	for(int i=0; i<size-1; i++){file1<<std::setprecision(4)<<x[i]<<"\t";}
	file1<<x[size-1]<<"\n";
	file1.close();
}

//void write_matrix_header (string* x, int size, string dest_folder, string dest_file,   ofstream & file1){
//
//
//}

void write_matrix_2 (int** x,  string* genes, string* genes_frag,string*cells,int iter, int G, int N, string dest_folder, string dest_file,   ofstream & file1, bool restart){

	if(restart){
		file1.open((dest_folder+dest_file).c_str());
		if(!file1) cerr<<"NOT able to write/append to "<<(dest_folder+dest_file).c_str()<<endl;
		file1<<"gene\tfeature\titeration\t";
		for(int i=0; i<N-1; i++){file1<<cells[i]<<"\t";}
		file1<<cells[N-1]<<"\n";
		file1.close();
	}

	file1.open((dest_folder+dest_file).c_str(),ios_base::app);
	if(!file1) cerr<<"NOT able to write/append to "<<(dest_folder+dest_file).c_str()<<endl;

	for(int g=0; g<G; g++){
		file1<<genes[g]<<"\t"<<genes_frag[g]<<"\t"<<iter<<"\t";
		for(int n=0; n<N-1; n++)
			file1<<x[g][n]<<"\t";
		file1<<x[g][N-1]<<"\n";
	}
	file1.close();
}

void write_one_record_per_gene (double* x,  string* genes, string* genes_frag, int iter, int G,  string dest_folder, string dest_file,   ofstream & file1, bool restart){

	if(restart){
		file1.open((dest_folder+dest_file).c_str());
		if(!file1) cerr<<"NOT able to write/append to "<<(dest_folder+dest_file).c_str()<<endl;
		file1<<"gene\tfeature\titeration\tvalue"<<"\n";
		file1.close();
	}

	file1.open((dest_folder+dest_file).c_str(),ios_base::app);
	if(!file1) cerr<<"NOT able to write/append to "<<(dest_folder+dest_file).c_str()<<endl;

	for(int g=0; g<G; g++)
		file1<<genes[g]<<"\t"<<genes_frag[g]<<"\t"<<iter<<"\t"<<x[g]<<"\n";
	file1.close();
}
void write_matrix_2 (float** x,  string*genes,string* genes_frag,string*cells,int iter, int G, int N, string dest_folder, string dest_file,   ofstream & file1, bool restart){

	if(restart){
		file1.open((dest_folder+dest_file).c_str());
		if(!file1) cerr<<"NOT able to write/append to "<<(dest_folder+dest_file).c_str()<<endl;
		file1<<"gene\tfeature\titeration\t";
		for(int i=0; i<N-1; i++){file1<<cells[i]<<"\t";}
		file1<<cells[N-1]<<"\n";
		file1.close();
	}

	file1.open((dest_folder+dest_file).c_str(),ios_base::app);
	if(!file1) cerr<<"NOT able to write/append to "<<(dest_folder+dest_file).c_str()<<endl;

	for(int g=0; g<G; g++){
		file1<<genes[g]<<"\t"<<genes_frag[g]<<"\t"<<iter<<"\t";
		for(int n=0; n<N-1; n++)
			file1<<x[g][n]<<"\t";
		file1<<x[g][N-1]<<"\n";
	}
	file1.close();
}

void write_matrix_2 (double** x,  string*genes,string* genes_frag,string*cells,int iter, int G, int N, string dest_folder, string dest_file,   ofstream & file1, bool restart){

	if(restart){
		file1.open((dest_folder+dest_file).c_str());
		if(!file1) cerr<<"NOT able to write/append to "<<(dest_folder+dest_file).c_str()<<endl;
		file1<<"gene\tfeature\titeration\t";
		for(int i=0; i<N-1; i++){file1<<cells[i]<<"\t";}
		file1<<cells[N-1]<<"\n";
		file1.close();
	}

	file1.open((dest_folder+dest_file).c_str(),ios_base::app);
	if(!file1) cerr<<"NOT able to write/append to "<<(dest_folder+dest_file).c_str()<<endl;

	for(int g=0; g<G; g++){
		file1<<genes[g]<<"\t"<<genes_frag[g]<<"\t"<<iter<<"\t";
		for(int n=0; n<N-1; n++)
			file1<<x[g][n]<<"\t";
		file1<<x[g][N-1]<<"\n";
	}
	file1.close();
}
void write_matrix (int** x,  string*genes,string*cells,int iter, int G, int N, string dest_folder, string dest_file,   ofstream & file1, bool restart){

	if(restart){
		file1.open((dest_folder+dest_file).c_str());
		if(!file1) cerr<<"NOT able to write/append to "<<(dest_folder+dest_file).c_str()<<endl;
		file1<<"feature\titeration\t";
		for(int i=0; i<N-1; i++){file1<<cells[i]<<"\t";}
		file1<<cells[N-1]<<"\n";
		file1.close();
	}

	file1.open((dest_folder+dest_file).c_str(),ios_base::app);
	if(!file1) cerr<<"NOT able to write/append to "<<(dest_folder+dest_file).c_str()<<endl;

	for(int g=0; g<G; g++){
		file1<<genes[g]<<"\t"<<iter<<"\t";
		for(int n=0; n<N-1; n++)
			file1<<x[g][n]<<"\t";
		file1<<x[g][N-1]<<"\n";
	}
	file1.close();
}

void write_matrix (float** x,  string*genes,string*cells,int iter, int G, int N, string dest_folder, string dest_file,   ofstream & file1, bool restart){

	if(restart){
		file1.open((dest_folder+dest_file).c_str());
		if(!file1) cerr<<"NOT able to write/append to "<<(dest_folder+dest_file).c_str()<<endl;
		file1<<"feature\titeration\t";
		for(int i=0; i<N-1; i++){file1<<cells[i]<<"\t";}
		file1<<cells[N-1]<<"\n";
		file1.close();
	}

	file1.open((dest_folder+dest_file).c_str(),ios_base::app);
	if(!file1) cerr<<"NOT able to write/append to "<<(dest_folder+dest_file).c_str()<<endl;

	for(int g=0; g<G; g++){
		file1<<genes[g]<<"\t"<<iter<<"\t";
		for(int n=0; n<N-1; n++)
			file1<<x[g][n]<<"\t";
		file1<<x[g][N-1]<<"\n";
	}
	file1.close();
}

void write_vec (double* x, int size, string dest_folder, string dest_file,   ofstream & file1, bool restart){

	if(restart) file1.open((dest_folder+dest_file).c_str());
	else file1.open((dest_folder+dest_file).c_str(),ios_base::app);
	if(!file1) cerr<<"NOT able to write/append to "<<(dest_folder+dest_file).c_str()<<endl;
	for(int i=0; i<size; i++){file1<<std::setprecision(4)<<x[i]<<"\t";}
	file1<<"\n";
	file1.close();
}

void write_vec (string* x, int size, string dest_folder, string dest_file,   ofstream & file1, bool restart){

	if(restart) file1.open((dest_folder+dest_file).c_str());
	else file1.open((dest_folder+dest_file).c_str(),ios_base::app);
	if(!file1) cerr<<"NOT able to write/append to "<<(dest_folder+dest_file).c_str()<<endl;
	for(int i=0; i<size; i++){file1<<x[i]<<"\t";}
	file1<<"\n";
	file1.close();
}

void write_value (float x,string dest_folder, string dest_file,   ofstream & file1, bool restart){

	if(restart) file1.open((dest_folder+dest_file).c_str());
	else file1.open((dest_folder+dest_file).c_str(),ios_base::app);
	if(!file1) cerr<<"NOT able to write/append to "<<(dest_folder+dest_file).c_str()<<endl;
	file1<<std::setprecision(4)<<x<<"\n";
	file1.close();
}

void read_parameters(string param_file,string* run_type,string* detail,int* K, int *inter_out, int *itermax, int* nb_flag, int* zetak_flag,\
		int* kappa_flag, int* kappa_w_flag,int* kappa_w2_flag, int* kappa_z_flag, int* kappan_flag, int* all_genes_flag, int* iter_start_cooling,\
		int* iter_cooling_step, double* amplitude_cooling_step, int* kappa_g_flag,int* iter_update_kappa,int* window_adjustment_flag){

	cout<< "Reading algorithm parameters" <<endl;
	ifstream file0;
	file0.open(param_file.c_str());
	if(!file0) {
		cout<<"The parameter file \"" << param_file << "\" is missing"<<endl;exit(2);
	}

	string keyword;
	file0 >> keyword;
	while (!file0.eof()) {
		if (file0.eof()) {
			cerr << "unexpected EOF after \""<< keyword <<"\" in parameter file" << endl;
			exit(1);
		} else if (keyword.compare("Number_clusters")==0) {
			if (*K>0) {
				cerr << "error: number of cluster already specified with -K option in command line args" << endl;
				exit(1);
			}
			file0 >> *K;
			cout << keyword << " " << *K << endl;
		} else if (keyword.compare("Run_type")==0) {
			file0 >> *run_type;
			cout << keyword << " " << *run_type << endl;
		} else if (keyword.compare("Run_detail")==0) {
			file0 >> *detail;
			cout << keyword << " " << *detail << endl;
		} else if (keyword.compare("Tinning_step")==0) {
			file0 >> *inter_out;
			cout << keyword << " " << *inter_out << endl;
		} else if (keyword.compare("Itermax")==0) {
			file0 >> *itermax;
			cout << keyword << " " << *itermax << endl;
		} else if (keyword.compare("iter_start_cooling")==0) {
			file0 >> *iter_start_cooling;
			cout << keyword << " " << *iter_start_cooling << endl;
		} else if (keyword.compare("iter_cooling_step")==0) {
			file0 >> *iter_cooling_step;
			cout << keyword << " " << *iter_cooling_step << endl;
		} else if (keyword.compare("amplitude_cooling_step")==0) {
			file0 >> *amplitude_cooling_step;
			cout << keyword << " " << *amplitude_cooling_step << endl;
		} else if (keyword.compare("Run_with_NB")==0) {
			file0 >> *nb_flag;
			cout << keyword << " " << *nb_flag << endl;
		} else if (keyword.compare("Use_zeta^k")==0) {
			file0 >> *zetak_flag;
			cout << keyword << " " << *zetak_flag << endl;
		} else if (keyword.compare("Use_kappa")==0) {
			file0 >> *kappa_flag;
			cout << keyword << " " << *kappa_flag << endl;
		} else if (keyword.compare("Use_kappa_w")==0) {
			file0 >> *kappa_w_flag;
			cout << keyword << " " << *kappa_w_flag << endl;
		}else if (keyword.compare("Use_kappa_w2")==0) {
			file0 >> *kappa_w2_flag;
			cout << keyword << " " << *kappa_w2_flag << endl;
		} else if (keyword.compare("Use_kappa_z")==0) {
			file0 >> *kappa_z_flag;
			cout << keyword << " " << *kappa_z_flag << endl;
		} else if (keyword.compare("Use_kappa^n")==0) {
			file0 >> *kappan_flag;
			cout << keyword << " " << *kappan_flag << endl;
		} else if (keyword.compare("Use_kappa_g")==0) {
			file0 >> *kappa_g_flag;
			cout << keyword << " " << *kappa_g_flag << endl;
		} else if (keyword.compare("Use_all_genes_d0")==0) {
			file0 >> *all_genes_flag;
			cout << keyword << " " << *all_genes_flag << endl;
		}else if (keyword.compare("Use_gene_length_w")==0) {
			file0 >> *window_adjustment_flag;
			cout << keyword << " " << *window_adjustment_flag << endl;
		} else if (keyword.compare("iter_update_kappa")==0) {
			file0 >> *iter_update_kappa;
			cout << keyword << " " << *iter_update_kappa << endl;
		} else {
			cerr << "unknown keyword \""<< keyword <<"\" in parameter file" << endl;
			exit(1);
		}
		file0 >> keyword;
	}
	file0.close();
	cout<<endl;
}


void simulate_data (gsl_rng* rgen, datas* dat){
	cout<<"Start simulating data"<<endl;
	dat->_GG=100;

	for (int n=0; n<dat->_N;n++)
		dat->_pn [n]=0.4;//(float)(n%20)/50.0+0.05;
	for (int n=0; n<dat->_N;n++)
		dat->_an [n]=0.1;//(float)(n%20)/50.0+0.05;
	for (int g=0; g<dat->_G;g++)
		dat->_bg [g]=0.3;//(float)(n%20)/50.0+0.05;
	for (int g=0; g<dat->_G;g++)
		dat->_window[g]=14;

	for (int g=0; g<dat->_G;g++)
		for (int l=0; l<dat->_K;l++)
			dat->_mu[g][l]=gsl_ran_gamma(rgen,dat->_prior_m_mu*dat->_prior_b_mu,1.0/dat->_prior_b_mu*(l+1));

	int* klusters=new int [dat->_N];
	for (int n=0; n<dat->_N;n++){
		klusters[n]=floor(n/((double)dat->_K+0.00001));
	}

	dat->_a=1;
	char* tempc=new char[1];
	for (int n=0; n<dat->_N;n++){
		tempc[0]='0'+n;
		dat->_cells [n]=string(("cell_"+string(tempc)+"_sim").c_str());//(float)(n%20)/50.0+0.05;
	}
	for (int g=0; g<dat->_G;g++){
		tempc[0]='0'+g;
		dat->_genes [g]=string(("gene_"+string(tempc)+"_sim").c_str());//(float)(n%20)/50.0+0.05;
	}

	for (int g=0; g<dat->_G;g++){
		tempc[0]='0'+g;
		dat->_genes_frag [g]=string(("gene_"+string(tempc)+"_sim").c_str());//(float)(n%20)/50.0+0.05;
	}


	for (int n=0; n<dat->_N;n++){
		for (int g=0; g<dat->_G;g++){

			if(g<20)
				dat->_z[g][n]=gsl_ran_poisson(rgen,(dat->_pn[n])*(dat->_mu[g][klusters[n]]));//  distribution of the sum of nb with mu and overdiseprsion   gsl_ran_negative_binomial_pdf(val,1/(mu*overdispersion+1),window/overdispersion)<<endl; exit(1);
			else
				dat->_z[g][n]=gsl_ran_poisson(rgen,(dat->_pn[n])*(dat->_mu[g][0]));

			if(dat->_z[g][n]==0) {
				dat->_counts[g][n]=0;
			} else {
				double cmean=dat->_z[g][n]*dat->_an[n]*dat->_bg[g]*dat->_window[g];
				double coverdisp=(dat->_kappa_n_comp[0][n] + dat->_kappa_n_comp[1][n]/dat->_window[g] + dat->_kappa_n_comp[2][n]/dat->_z[g][n])*dat->_kappa_g[g];
				dat->_counts[g][n]= gsl_ran_negative_binomial(rgen,1.0/(cmean*coverdisp+1),1.0/coverdisp);
			}
		}
	}

	ofstream FILE;
	int *temp=new int[dat->_N+3];
	temp[0]=dat->_GG; temp[1]=dat->_G; temp[2]=dat->_N;
	write_vec (temp, 3,  "/media/storage/DATA/Results/Single_cell/data/",  "counts_expression_sim",   FILE, 1);
	for( int i=3; i<(dat->_N+3);i++)
		temp[i]=i-2;
	write_vec (temp, dat->_N+3,  "/media/storage/DATA/Results/Single_cell/data/",  "counts_expression_sim",   FILE, 0);

	for (int g=0; g<dat->_G;g++){
		temp[0]=g+1;temp[1]=g+1;  temp[2]=dat->_window[g];for( int i=3; i<(dat->_N+3);i++)temp[i]=dat->_counts[g][i-3];
		write_vec (temp, dat->_N+3,  "/media/storage/DATA/Results/Single_cell/data/",  "counts_expression_sim",   FILE, 0);
	}
	//cout<<"End simulating data"<<endl;
}
