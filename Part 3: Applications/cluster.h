#ifndef CLUSTER_H
#define CLUSTER_H

#include "curve.h"

class cluster
{
    public:
        cluster(string id);
        virtual ~cluster();

        void print(string centroid_name);
        string get_id();
        unsigned int get_cluster_size();

        void update_cost(double cost);

        double get_total_cost();
        void set_total_cost(double cost);

        vector<curve *> get_members();
        void add_member(curve * curve_ptr);
        void delete_members();

    private:
        string id;
        double total_cost; // the dists' sum from all cluster's members to their centroid
        vector<curve *> members;
};

#endif // CLUSTER_H
