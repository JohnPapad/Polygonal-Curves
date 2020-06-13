#ifndef CURVE_H
#define CURVE_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <limits>

using namespace std;

class cluster;

class curve
{
    public:
        curve(vector<vector<double>> &points, vector<double> &c_point, string type,
        unsigned int id, unsigned int number_of_points, unsigned int dimension);

        curve(vector<vector<double>> &points, unsigned int way_id, unsigned int number_of_points, unsigned int dimension);

        void print(bool print_coords);
        vector<double>& get_c_point();
        vector<vector<double>>& get_original_curve_points();
        string get_type();
        unsigned int get_number_of_points();
        unsigned int get_curve_dimension();
        unsigned int get_id();


        void set_cluster(cluster * cl);
        void set_is_centroid(bool flag);
        void set_is_assigned(bool flag);
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
        vector<double> c_point;


        string type;
        unsigned int id;
        unsigned int number_of_points;
        unsigned int dimension;

        bool centroid_flag; //whether it is a cluster's centroid or not
      //  unsigned int index ; //index for accessing the array that contains all the dfds
        bool assigned_flag; //whether it has been assigned or not
        cluster * is_centroid_of_cluster ; //a pointer to its cluster (only if the curve is a centroid)

        curve * assigned_centroid ; //the curve's centroid (if it is a cluster's member)
        curve * poss_assigned_centroid; //the second closest centroid ( if it is a cluster's member)
        double dist_to_centroid; //the dfd between itself and its centroid
        double dist_to_poss_centroid; //the dfd between itself and its second closest centroid

};


#endif // CURVE_H
