#ifndef CURVE_H
#define CURVE_H

#include <iostream>
#include <vector>
#include <string>
#include <limits>

using namespace std;

struct q_curve_stats;

class curve
{
    public:
        ~curve();
        curve(vector<double> &original_curve, vector<double> &grid_curve, string &curve_name,
              unsigned int &number_of_points, unsigned int &dimension);
        void print_curve();
        vector<double>& get_original_curve();
        vector<double>& get_grid_curve();
        string get_curve_name();
        unsigned int get_number_of_points();
        unsigned int get_curve_dimension();
        bool compare_grid_curves(curve * query_ptr);

        void create_q_stats();
        void set_stat_found_grid_curve(int &found_grid_curve);
        void set_stat_lsh_nn(curve * curve_ptr);
        void set_stat_real_nn(curve * curve_ptr);
        void set_stat_lsh_dist(double &lsh_dist);
        void set_stat_true_dist(double &true_dist);
        void set_stat_t_lsh(double &t_lsh);
        void set_stat_t_true(double &t_true);
        void set_stat_t_lsh_min(double & t_lsh_min);
        void set_stat_t_lsh_max(double & t_lsh_max);
        void set_stat_t_lsh_avg(double & t_lsh_avg);
        void set_stat_min_dist_lsh(double &min_dist_lsh);
        void set_stat_max_dist_lsh(double &max_dist_lsh);
        void set_stat_avg_dist_lsh(double &avg_dist_lsh);

        int    get_stat_found_grid_curve();
        curve* get_stat_lsh_nn();
        curve* get_stat_real_nn();
        double get_stat_lsh_dist();
        double get_stat_true_dist();
        double get_stat_t_lsh();
        double get_stat_t_true();
        double get_stat_t_lsh_min();
        double get_stat_t_lsh_max();
        double get_stat_t_lsh_avg();
        double get_stat_min_dist_lsh();
        double get_stat_max_dist_lsh();
        double get_stat_avg_dist_lsh();


    private:
        vector<double> original_curve ;
        vector<double> grid_curve;
        string curve_name;
        unsigned int number_of_points;
        unsigned int dimension;
        q_curve_stats * q_stats_ptr ;
};


struct q_curve_stats
{
    int found_grid_curve = -1;
    curve * lsh_nn = NULL; //lsh_nearest_neighbor
    curve * real_nn = NULL; //real_nearest_neighbor
    double lsh_dist = -1.0 ; //distance between query and lsh nearest neighbor
    double true_dist = -1.0 ; //distance between query and true nearest neighbor
    double t_lsh = -1.0; //time spent for finding the lsh nearest neighbor
    double t_true = -1.0; //time spent for finding the true nearest neighbor
    double t_lsh_min = -1.0; // the minimum time spent for finding the lsh nearest neighbor (when the program has been repeated 100 times)
    double t_lsh_max = -1.0; // the maximum time spent for finding the lsh nearest neighbor (when the program has been repeated 100 times)
    double t_lsh_avg = -1.0; // the average time spent for finding the lsh nearest neighbor (when the program has been repeated 100 times)
    double min_dist_lsh = -1.0; //the minimum distance between the query and the lsh nearest neighbor (when the program has been repeated 100 times)
    double max_dist_lsh = -1.0; //the maximum distance between the query and the lsh nearest neighbor (when the program has been repeated 100 times)
    double avg_dist_lsh = -1.0; //the average distance between the query and the lsh nearest neighbor (when the program has been repeated 100 times)
};

#endif // CURVE_H
