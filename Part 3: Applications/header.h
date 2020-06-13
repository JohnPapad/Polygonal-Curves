
#ifndef HEADER_H
#define HEADER_H
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include "cluster.h"

int get_uniformly_distributed_random_int(int M, int N);
double get_uniformly_distributed_random_double2(double M, double N);
void initialize_random_seed(); //for rand function
void parse_proteins(string input_fn, vector<curve *> &proteins, unsigned int d);
void print_all_curves(vector<curve *> &curves, bool print_coords);
double calc_dist(curve * x_ptr, curve * y_ptr, bool dist_choice);
void print_points(vector<vector<double>> &original_curve_points, unsigned int dimension);
vector<vector<double>> select_random_points(vector<vector<double>> points, unsigned int n_points);
vector<curve *> random_k_points_selection(vector<curve *> curves, unsigned int k);
void simple_clustering(vector<curve *>& curves, vector<curve *>& centroids, vector<vector<double> > &dists, bool dist_choice);
vector<curve *> k_means_pp(vector<curve *> curves, vector<vector<double> > &dists, unsigned int k, bool dist_choice);
void lloyed_assignment(vector<curve *> &centroids, curve * curve_ptr, vector<vector<double> > &dists, bool dist_choice);
bool PAM(vector<curve *> &curves, vector<curve *> &centroids, vector<vector<double> > &dists, double bound, bool dist_choice);
double get_objective_fun_cost(vector<curve *> &curves, vector<curve *> &temp_centroids, vector<vector<double> > &dists, bool dist_choice);
double get_cluster_cost(vector<curve *> &cluster_members, curve * checking_curve, vector<vector<double> > &dists, int checking_member_index);
double get_silhouette(vector<curve *> &centroids, vector<vector<double> > &dists, vector<double> &clusters_silhts, bool dist_choice);
double get_all_clusters_objective_cost(vector<curve *> &centroids);
double eucl_distance(vector<double> &point1, vector<double> &point2);
double get_dist(curve * curve1_ptr, curve * curve2_ptr, bool dist_choice, vector<vector<double> > &dists);
double get_min_dist(vector<curve * > &curves, curve * query_ptr, vector<vector<double> > &dists, double &sec_min_dist, unsigned int &min_index, unsigned int &sec_min_index, bool dist_choice);
void reset_cluster(cluster * cluster_ptr);
void reset_all_curves(vector<curve *> &curves);
void reset_curves(vector<curve *> &curves);
void reset_all_clusters(vector<curve *> &centroids);
void reset_curve(curve * curve_ptr);
void set_new_centroid(curve * new_centroid, cluster * cluster_ptr);
void parse_roads(string &file_name, vector<curve *> &roads, bool skip_first_line, unsigned int d);
double frechet_distance(vector<vector<double>> &curve1, vector<vector<double>> &curve2);
bool read_flag() ;
void create_clusters(vector<curve *>& centroids);
void destroy_clusters(vector<curve *>& centroids, bool mean_flag);
void write_results_in_file(vector<curve *> &centroids, double overall_silht, string output_fn, double clustering_time, bool segments_flag);
void print_counters();
void delete_curves(vector<curve * > &curves);
vector<vector<double>> cut_off_points(vector<vector<double>> points, unsigned int n_points);
bool compare_nodes(vector<double> &node1, vector<double> &node2);
bool check_range(vector<double> &first_node, vector<double> &last_node, vector<double> &checking_node);
bool check_for_junction(vector<curve *> &roads, vector<double> &checking_node, string checking_node_type, unsigned int checking_node_i);
void create_segments_file(vector<curve *> &segs, string output_fn);
double get_circumcircle_radius(vector<double> &point_a, vector<double> &point_b, vector<double> &point_c);
void seperate_roads_to_segments(vector<curve *> &roads, vector<curve *> &segs, double radius_bound, unsigned int min_seg_size, unsigned int max_seg_size);
vector<double> calc_point_c( vector<vector<double>> &points);


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



#endif // HEADER_H
