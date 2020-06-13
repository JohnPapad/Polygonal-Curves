
#include "header.h"

vector<curve *> random_k_points_selection(vector<curve *> curves, unsigned int k)
{
    vector<curve *> centroids;

    for(unsigned int i = 0 ; i < k ; i++)
    {
        unsigned int centroid_index = get_uniformly_distributed_random_int(0, curves.size() -1 ) ;
        curve * centroid_ptr = curves[centroid_index];
        centroids.push_back(centroid_ptr) ;
        curves.erase(curves.begin() + centroid_index) ;
    }

    return centroids ;
}


vector<curve *> k_means_pp(vector<curve *> curves, vector<vector<double> > &dists, unsigned int k, bool dist_choice)
{
    double double_null = 0.0 ;
    unsigned int int_null = 0 ;
    vector<curve *> centroids;
    curve* centroid_ptr ;

    //select the first one randomly
    unsigned int first_centroid_index = get_uniformly_distributed_random_int(0, curves.size() -1 ) ;
    centroid_ptr = curves[first_centroid_index];
    centroids.push_back(centroid_ptr);
    curves.erase(curves.begin() + first_centroid_index) ;

    if (k == 1)
    {
        return centroids;
    }

    unsigned int k_counter = 1;
    do
    {

        vector<pair<double, int>> Probs_CurveIndex ; //use pairs for sorting and using upper_bound(sorting will be done according to the dfd)
        double PartialSum = 0.0 ;
        for(unsigned int r = 0 ; r < curves.size(); r++)
        { //get the partial sum of the first r curves
            double min_dist = get_min_dist(centroids, curves[r], dists, double_null, int_null, int_null, dist_choice);
            PartialSum += min_dist * min_dist ;
            Probs_CurveIndex.emplace_back(PartialSum, r) ;
        }

        sort(Probs_CurveIndex.begin(), Probs_CurveIndex.end()) ;
        double x = get_uniformly_distributed_random_double2(0, Probs_CurveIndex[Probs_CurveIndex.size() -1 ].first) ;
        unsigned int upper_bound_index = upper_bound(Probs_CurveIndex.begin(), Probs_CurveIndex.end(), make_pair(x, 0)) - Probs_CurveIndex.begin() ;
        unsigned int r_index = Probs_CurveIndex[upper_bound_index].second ;
        centroid_ptr = curves[r_index];
        centroids.push_back(centroid_ptr);
        curves.erase(curves.begin() + r_index) ;
        k_counter++ ;
    }while(k_counter < k );

    return centroids ;
}


void simple_clustering(vector<curve *>& curves, vector<curve *>& centroids, vector<vector<double> > &dists, bool dist_choice)
{
    for(unsigned int i = 0; i < curves.size(); i++)
    { //for each curve that has not been marked as a centroid
        if ( curves[i]->is_centroid() == 0)
        {
            lloyed_assignment(centroids, curves[i], dists, dist_choice);
        }
    }
}

void lloyed_assignment(vector<curve *> &centroids, curve * curve_ptr, vector<vector<double> > &dists, bool dist_choice)
{
    double min_dist;
    double sec_min_dist;
    unsigned int centroid_index ;
    unsigned int poss_centroid_index;
    curve* centroid_ptr ;
    curve* poss_centroid_ptr;

    min_dist = get_min_dist(centroids, curve_ptr, dists, sec_min_dist, centroid_index, poss_centroid_index, dist_choice);
    centroid_ptr = centroids[centroid_index];
    poss_centroid_ptr = centroids[poss_centroid_index];

    curve_ptr->set_curve_centroid(centroid_ptr); //set its closest centroid
    curve_ptr->set_dist_to_centroid(min_dist); //set its dfd with its closest centroid

    centroid_ptr->get_cluster()->add_member(curve_ptr); //add curve to the centroid's cluster
    centroid_ptr->get_cluster()->update_cost(min_dist); // update the cluster's sum of dfds

    curve_ptr->set_curve_poss_centroid(poss_centroid_ptr); //set its second closest centroid
    curve_ptr->set_dist_to_poss_centroid(sec_min_dist); // set its dfd with the second closest centroid

}

bool PAM(vector<curve *> &curves, vector<curve *> &centroids, vector<vector<double> > &dists, double bound, bool dist_choice)
{
    vector<curve *> new_centroids = centroids;
    double total_cost = get_all_clusters_objective_cost(centroids);
    double min_objective_cost = numeric_limits<double>::infinity();

    for(unsigned int i = 0 ; i < centroids.size(); i++)
    {
        vector<curve *> cluster_members = centroids[i]->get_cluster()->get_members();
        vector<curve *> temp_centroids = centroids;

        for(unsigned int j = 0; j < cluster_members.size(); j++)
        {
            cluster_members[j]->set_is_centroid(1);
            temp_centroids[i] = cluster_members[j]; //set a cluster's member as a possible cluster's centroid

            double cost = get_objective_fun_cost(curves, temp_centroids, dists, dist_choice);
            cluster_members[j]->set_is_centroid(0);
            unsigned int int_null;
            double double_null;
            double dist = get_min_dist(temp_centroids, centroids[i], dists, double_null, int_null, int_null, dist_choice);

            cost += dist;
            if (cost < min_objective_cost) //same logic as in Mean Discrete Frechet curve
            {
                min_objective_cost = cost;
                new_centroids.clear();
                new_centroids = temp_centroids;
            }
        }
    }

    double cost_diff = total_cost - min_objective_cost;
    if (cost_diff > bound)
    { //exact same logic as in Mean Discrete Frechet curve
        for(unsigned int i = 0; i < new_centroids.size(); i++)
        {
            if (new_centroids[i]->get_id() != centroids[i]->get_id())
            {
                set_new_centroid(new_centroids[i], centroids[i]->get_cluster());
                reset_curve(centroids[i]);
                centroids[i] = new_centroids[i];
            }
        }
        return 0;
    }
    else
    {
        return 1;
    }
}



void create_clusters(vector<curve *>& centroids)
{ //for each centroid create a class cluster to assign its pointer to centroid
    for(unsigned int i=0; i < centroids.size(); i++)
    {
        stringstream ss;
        unsigned int int_id = i + 1;
        ss << int_id;
        string cluster_id = ss.str();
        centroids[i]->set_is_centroid(1); //mark curve as centroid
        cluster * cluster_ptr = new cluster(cluster_id);
        centroids[i]->set_cluster(cluster_ptr); //assign cluster's pointer to curve
    }
}


void destroy_clusters(vector<curve *>& centroids, bool mean_flag)
{
    for(unsigned int i=0; i < centroids.size(); i++)
    {
        delete centroids[i]->get_cluster();
        if (mean_flag == 1) //Mean Discrete Curve has been selected as an update method
        {
            delete centroids[i]; //so we must delete centroids too because they are members of the dataset and have been created afterwards
        }

        centroids[i] = NULL;
    }
    centroids.clear();

}



void reset_cluster(cluster * cluster_ptr)
{
    cluster_ptr->set_total_cost(0.0);
    cluster_ptr->delete_members();
}

void reset_all_curves(vector<curve *> &curves)
{
    for(unsigned int i=0; i < curves.size(); i++)
    {
        if (curves[i]->is_centroid() == 0) //if a curve is not marked as centroid reset it
        {
            reset_curve(curves[i]);
        }
    }
}

void reset_curves(vector<curve *> &curves)
{
    for(unsigned int i=0; i < curves.size(); i++)
    {
        reset_curve(curves[i]);
    }
}

void reset_all_clusters(vector<curve *> &centroids)
{
    for(unsigned int i=0; i < centroids.size(); i++)
    {
        reset_cluster(centroids[i]->get_cluster());
    }
}


void reset_curve(curve * curve_ptr)
{ //set original default values
    curve_ptr->set_is_centroid(0);
    curve_ptr->set_is_assigned(0);
    curve_ptr->set_cluster(NULL);
    curve_ptr->set_curve_centroid(NULL);
    curve_ptr->set_curve_poss_centroid(NULL);
    curve_ptr->set_dist_to_centroid(numeric_limits<double>::infinity());
    curve_ptr->set_dist_to_poss_centroid(numeric_limits<double>::infinity());
}


void set_new_centroid(curve * new_centroid, cluster * cluster_ptr)
{
    new_centroid->set_is_centroid(1);
    new_centroid->set_cluster(cluster_ptr);
}


double get_objective_fun_cost(vector<curve *> &curves, vector<curve *> &temp_centroids, vector<vector<double> > &dists, bool dist_choice)
{
    double total_cost = 0.0 ;
    for(unsigned int i= 0; i < curves.size(); i++)
    {
        if (curves[i]->is_centroid()== 0)
        {
            double sec_min_dist;
            unsigned int centroid_index ;
            unsigned int poss_centroid_index;

            double min_dist = get_min_dist(temp_centroids, curves[i], dists, sec_min_dist, centroid_index, poss_centroid_index, dist_choice);
            total_cost += min_dist ;
        }
    }

    return total_cost;
}

double get_cluster_cost(vector<curve *> &cluster_members, curve * checking_curve, vector<vector<double> > &dists, int checking_member_index, bool dist_choice)
{ // get the sum of dfds between a curve and a cluster
    double cost = 0.0 ;
    int cluster_members_size = cluster_members.size();
    for(int i=0; i < cluster_members_size; i++)
    {
        if (i == checking_member_index) // if the curve belongs to the cluster to calculate the dfd with itself
        {
            continue;
        }

        cost += get_dist(cluster_members[i], checking_curve, dist_choice, dists);
    }

    return cost;
}


double get_silhouette(vector<curve *> &centroids, vector<vector<double> > &dists, vector<double> &clusters_silhts, bool dist_choice)
{
    double overall_silht = 0;
    for(unsigned int i =0; i < centroids.size(); i++)
    {
        vector<curve *> cluster_members = centroids[i]->get_cluster()->get_members();
        double s_i  = 0; //cluster's silhouette

        for(unsigned int j=0; j < cluster_members.size(); j++)
        {
            vector<curve *> poss_cluster_members;

            double cluster_cost;
            double poss_cluster_cost ;

            if (cluster_members[j]->get_curve_poss_centroid() != NULL)
            { //if second closest centroid has already been found take its cluster's member
                poss_cluster_members = cluster_members[j]->get_curve_poss_centroid()->get_cluster()->get_members();
            }
            else //if not, find it
            {
                double double_null;
                unsigned int real_centroid_index;
                unsigned int poss_centroid_index;
                get_min_dist(centroids, cluster_members[j], dists, double_null, real_centroid_index, poss_centroid_index, dist_choice);

                curve * real_centroid = centroids[real_centroid_index];
                curve * poss_centroid = centroids[poss_centroid_index];
                curve * cluster_member_centroid = cluster_members[j]->get_curve_centroid();

                if (cluster_member_centroid->get_id() == poss_centroid->get_id())
                { //if curve has been mistakenly assigned to the second closest centroid (this may happen when LSH is selected)
                    poss_cluster_members = real_centroid->get_cluster()->get_members();
                }
                else
                {
                    poss_cluster_members = poss_centroid->get_cluster()->get_members();
                }
            }

            double b_j;
            if ( poss_cluster_members.size() == 0)
            {
                b_j = 0;
            }
            else
            {
                poss_cluster_cost = get_cluster_cost(poss_cluster_members, cluster_members[j], dists, -1, dist_choice);
                b_j = poss_cluster_cost / (poss_cluster_members.size());
            }

            cluster_cost = get_cluster_cost(cluster_members, cluster_members[j], dists, j, dist_choice);
            double a_j = cluster_cost / (cluster_members.size());

            double max_ab = a_j > b_j ? a_j : b_j ;
            double s_j = (b_j - a_j) / max_ab;

            s_i += s_j;
        }

        if ( cluster_members.size() != 0)
        {
            s_i = s_i / cluster_members.size() ;
        }
        else
        {
            s_i = 0;
        }
        clusters_silhts.push_back(s_i);
        overall_silht += s_i ;
    }

    overall_silht = overall_silht / centroids.size();
    return overall_silht;
}

double get_all_clusters_objective_cost(vector<curve *> &centroids)
{ //sum all objective costs from each cluster
    double total_cost = 0.0 ;
    for(unsigned int i=0 ; i < centroids.size() ; i++)
    {
        total_cost += centroids[i]->get_cluster()->get_total_cost();
    }

    return total_cost;
}
