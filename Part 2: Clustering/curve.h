#ifndef CURVE_H
#define CURVE_H

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>

using namespace std;

class cluster;

class curve
{
    public:
        ~curve();
        curve(vector<double> &original_curve, vector<vector<double>> &original_curve_points, vector<double> &grid_curve, string curve_name,
              unsigned int number_of_points, unsigned int dimension);
        curve(curve * curve_ptr, unsigned int index);
        curve(vector<vector<double>> &original_curve_points, unsigned int number_of_points, unsigned int dimension);

        void print_curve(bool print_coords);
        vector<double>& get_original_curve();
        vector<vector<double>>& get_original_curve_points();
        vector<double>& get_grid_curve();
        string get_curve_name();
        unsigned int get_number_of_points();
        unsigned int get_curve_dimension();
        bool compare_grid_curves(curve * query_ptr);
        void set_original_curve(vector<double> &original_curve);
        void set_grid_curve(vector<double> &grid_curve);
        void set_curve_name(string name);

        void set_cluster(cluster * cl);
        void set_is_centroid(bool flag);
        void set_is_assigned(bool flag);
        void set_index(unsigned int index);
        void set_curve_centroid(curve * centroid);
        void set_curve_poss_centroid(curve * centroid);
        void set_dist_to_centroid(double dist);
        void set_dist_to_poss_centroid(double dist);
        double get_dist_to_centroid();
        double get_dist_to_poss_centroid();
        curve* get_curve_centroid();
        curve* get_curve_poss_centroid();
        cluster * get_cluster();
        bool is_centroid();
        unsigned int get_index();
        bool is_assigned();


    private:
        vector<vector<double>> original_curve_points;

        vector<double> original_curve ;
        vector<double> grid_curve;
        string curve_name;
        unsigned int number_of_points;
        unsigned int dimension;

        bool centroid_flag; //whether it is a cluster's centroid or not
        unsigned int index ; //index for accessing the array that contains all the dfds
        bool assigned_flag; //whether it has been assigned or not
        cluster * is_centroid_of_cluster ; //a pointer to its cluster (only if the curve is a centroid)

        curve * assigned_centroid ; //the curve's centroid (if it is a cluster's member)
        curve * poss_assigned_centroid; //the second closest centroid ( if it is a cluster's member)
        double dist_to_centroid; //the dfd between itself and its centroid
        double dist_to_poss_centroid; //the dfd between itself and its second closest centroid

};


#endif // CURVE_H
