#ifndef AFUNCTIONS_H
#define AFUNCTIONS_H

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

# include "datas.h"

using namespace std;
//void write_matrix_header (string* x, int size, string dest_folder, string dest_file,   ofstream & file1);
void write_matrix(int** x,  string*genes,string*cells,int iter, int G, int N, string dest_folder, string dest_file,   ofstream & file1, bool restart);
void write_matrix(float** x,  string*genes,string*cells,int iter, int G, int N, string dest_folder, string dest_file,   ofstream & file1, bool restart);
void write_matrix_2(int** x,  string*genes,string*genes_frag,string*cells,int iter, int G, int N, string dest_folder, string dest_file,   ofstream & file1, bool restart);
void write_matrix_2(float** x,  string*genes,string*genes_frag,string*cells,int iter, int G, int N, string dest_folder, string dest_file,   ofstream & file1, bool restart);
void write_matrix_2(double** x,  string*genes,string*genes_frag,string*cells,int iter, int G, int N, string dest_folder, string dest_file,   ofstream & file1, bool restart);
void write_vec (float* x, int size, string dest_folder, string dest_file,   ofstream & file1, bool restart);
void write_vec (int* x, int size, string dest_folder, string dest_file,   ofstream & file1, bool restart);
void write_vec (double* x, int size, string dest_folder, string dest_file,   ofstream & file1, bool restart);
void write_vec (string* x, int size, string dest_folder, string dest_file,   ofstream & file1, bool restart);
void write_value (float x, string dest_folder, string dest_file,   ofstream & file1, bool restart);
void write_one_record_per_gene (double* x,  string* genes, string* genes_frag, int iter, int G,  string dest_folder, string dest_file,   ofstream & file1, bool restart);
void read_parameters(string param_file,string* run_type,string* detail,int* K, int *inter_out, int *itermax,int* nb_flag, int* zetak_flag, int* kappa_flag, int* kappa_w_flag,int* kappa_w2_flag, int* kappa_z_flag, int* kappan_flag,int *all_genes_flag, int* iter_start_cooling, int* iter_cooling_step, double* amplitude_cooling_step, int* kappa_g_flag, int* iter_update_kappa,int* window_adjustment_flag);
void simulate_data (gsl_rng* rgen, datas* dat);
#endif
