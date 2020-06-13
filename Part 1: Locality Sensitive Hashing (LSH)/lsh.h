#ifndef LSH_H
#define LSH_H

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string.h>
#include <ctime>
#include "curve.h"
#include "multi_hash_table.h"

using namespace std;


void parsing_saving_curves(string &file_name, vector<curve *>&curves, double &delta, bool q_curve_flag );

double coord_to_grid_coord_transform(double &delta, double &coord);

vector <double> curve_to_grid_curve_transform(vector <double>& curve, unsigned int d, double &delta);

unsigned int get_HT_key(vector<int>& r_vales, curve * curve_ptr, unsigned int &k, double &delta, unsigned int &HT_SIZE);

double get_uniformly_distributed_random_double(double M, double N);

int get_uniformly_distributed_random_int(int M, int N);

void initialize_random_seed();

unsigned long int get_concat_grid_curve_HV(vector<double> &concat_grid_curve, vector<int>& r_values, unsigned long int &M );

vector<double> create_shifted_curve(vector<double> &curve, double &delta);

void print_all_curves(vector<curve *> &curves);

multi_hash_table* create_multi_hash_table(unsigned int &L, unsigned int &HT_SIZE);

void destroy_multi_hash_table(multi_hash_table* MHT_ptr);

void delete_curves(vector<curve * > &curves);

void concatenating_grids_curves(vector<double>& grid_curve, vector<double>& concat_grid_curve);

unsigned int get_curve_dimension(string points_str);

vector <double>  remove_consecutive_duplicate_points(vector <double> dubl_curve, unsigned int d);

void LSH(vector< int >& r_values, multi_hash_table* MHT_ptr, curve * curve_ptr, unsigned int &L,
         unsigned int &k, unsigned int &HT_SIZE, double &delta);

void get_lsh_nearest_neighbor(vector<int>& r_values, multi_hash_table* MHT_ptr, curve * query_ptr,
                          unsigned int &k, double &delta, unsigned int &HT_SIZE);

void get_real_nearest_neighbor(vector<curve *> &curves, curve * query_ptr);

double frechet_distance(vector<double> &curve1, vector<double> &curve2, unsigned int d);

double eucl_distance(vector<double> &point1, vector<double> &point2);

vector<double> get_a_point_from_curve(vector<double> &curve, unsigned int point_index, unsigned int d);

void remove_dublicate_curves(vector<curve *> &matched_curves);

double get_min_DFD(vector<curve * > &curves, curve * query_ptr, unsigned int &min_index, double &elapsed_time);

void write_results_in_file(ofstream &outfile, curve * query_ptr);

void write_stats_in_file(ofstream &outfile, curve * query_ptr);

void run_without_stats(vector<curve *> &curves, string &q_data_set_file_name, string &output_file_name,
                       unsigned int &L, unsigned int &k, unsigned int &HT_SIZE, double &delta);

void run_with_stats(vector<curve *> &curves, string &q_data_set_file_name, string &output_file_name,
                    unsigned int &L, unsigned int &k, unsigned int &HT_SIZE, unsigned int &repeats, double &delta);

bool repeat(string &q_data_set_file_name, string &output_file_name);

void check_argv_string(string &argv, string proper_argv, string flag_option);


class InputParser{  // a class for parsing command line arguements and their options
    public:
        InputParser (int &argc, char **argv){
            for (int i=1; i < argc; ++i)
                this->tokens.push_back(string(argv[i]));
        }

        const string& getCmdOption(const string &option) const{
            vector<string>::const_iterator itr;
            itr = find(this->tokens.begin(), this->tokens.end(), option);
            if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                return *itr;
            }
            static const string empty_string("");
            return empty_string;
        }

        bool cmdOptionExists(const string &option) const{
            return find(this->tokens.begin(), this->tokens.end(), option) != this->tokens.end();
        }
    private:
        vector <string> tokens;
};

//bool compare_curves(vector<double>& curve1, vector<double>& curve2);
//hash_table* create_hash_table(unsigned int &HT_SIZE);
//void destroy_hash_table(hash_table* HT_ptr);


#endif // LSH_H
