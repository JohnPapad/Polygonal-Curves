#include "curve.h"

curve::curve(vector<double> &original_curve, vector<double> &grid_curve, string &curve_name,
             unsigned int &number_of_points, unsigned int &dimension)
{
    this->original_curve = original_curve;
    this->grid_curve = grid_curve;
    this->dimension = dimension;
    this->curve_name = curve_name;
    this->number_of_points = number_of_points;
    this->q_stats_ptr = NULL;
}


curve::~curve() // if q_curve_stats instance has been dynamically created for query curve we have to delete it manually
{
    if (q_stats_ptr != NULL)
    {
        delete q_stats_ptr;
    }
}


void curve::print_curve()
{
    typedef numeric_limits< double > dbl;
    cout<<"----------------------------------------"<<endl;
    cout<<endl;
    cout<<"PRINTING ORIGINAL CURVE"<<endl;
    cout<<"original curve size="<<original_curve.size()<<endl;
    cout<<"(";
    cout.precision(dbl::max_digits10);
    for(unsigned int i=0; i < original_curve.size();i++)
    {
        cout<<fixed<<original_curve[i];
        if (original_curve.size()-1 != i)
            cout<<",";
    }
    cout<<")"<<endl;
    cout<<endl;

    cout<<"PRINTING GRID CURVE"<<endl;
    cout<<"grid curve size="<<grid_curve.size()<<endl;
    cout<<"(";
    cout.precision(4);
    for(unsigned int i=0; i < grid_curve.size();i++)
    {
        cout<<fixed<<grid_curve[i];
        if (grid_curve.size()-1 != i)
            cout<<",";
    }
    cout<<")"<<endl;

    cout<<"curve id:"<<curve_name<<endl;
    cout<<"curve's number of points:"<<number_of_points<<endl;
    cout<<"curve's dimension="<<dimension<<endl;
    cout<<"----------------------------------------"<<endl;
    cout<<endl;
}


vector<double>& curve::get_grid_curve()
{
    return grid_curve;
}


vector<double>& curve::get_original_curve()
{
    return original_curve;
}


unsigned int curve::get_curve_dimension()
{
    return dimension;
}


string curve::get_curve_name()
{
    return curve_name;
}


bool curve::compare_grid_curves(curve * query_ptr) //compare this curve with another by their grid curves
{
    vector<double> grid_curve = query_ptr->get_grid_curve();

    if (this->grid_curve.size() != grid_curve.size()) //the 2 grid curves must have the same number of points
    {
        return 0;
    }

    for(unsigned int i=0; i < grid_curve.size(); i++)
    {
        if ( this->grid_curve[i] != grid_curve[i])  //all coordinates must be the same
        {
            return 0;
        }
    }

    return 1;
}


void curve::create_q_stats()
{
    if ( q_stats_ptr == NULL)
    {
        q_stats_ptr = new q_curve_stats;
    }
}


void curve::set_stat_found_grid_curve(int &found_grid_curve)
{
    if ( q_stats_ptr != NULL)
    {
        q_stats_ptr->found_grid_curve = found_grid_curve;
    }
}


void curve::set_stat_lsh_nn(curve * curve_ptr)
{
    if ( q_stats_ptr != NULL)
    {
        q_stats_ptr->lsh_nn = curve_ptr;
    }
}


void curve::set_stat_real_nn(curve * curve_ptr)
{
    if ( q_stats_ptr != NULL)
    {
        q_stats_ptr->real_nn = curve_ptr;
    }
}


void curve::set_stat_lsh_dist(double &lsh_dist)
{
    if ( q_stats_ptr != NULL)
    {
        q_stats_ptr->lsh_dist = lsh_dist;
    }
}


void curve::set_stat_true_dist(double &true_dist)
{
    if ( q_stats_ptr != NULL)
    {
        q_stats_ptr->true_dist = true_dist;
    }
}


void curve::set_stat_t_lsh(double &t_lsh)
{
    if ( q_stats_ptr != NULL)
    {
        q_stats_ptr->t_lsh = t_lsh;
    }
}


void curve::set_stat_t_true(double &t_true)
{
    if ( q_stats_ptr != NULL)
    {
        q_stats_ptr->t_true = t_true;
    }
}


void curve::set_stat_t_lsh_min(double & t_lsh_min)
{
    if ( q_stats_ptr != NULL)
    {
        q_stats_ptr->t_lsh_min = t_lsh_min;
    }
}


void curve::set_stat_t_lsh_max(double & t_lsh_max)
{
    if ( q_stats_ptr != NULL)
    {
        q_stats_ptr->t_lsh_max = t_lsh_max;
    }
}


void curve::set_stat_t_lsh_avg(double & t_lsh_avg)
{
    if ( q_stats_ptr != NULL)
    {
        q_stats_ptr->t_lsh_avg = t_lsh_avg;
    }
}


void curve::set_stat_min_dist_lsh(double &min_dist_lsh)
{
    if ( q_stats_ptr != NULL)
    {
        q_stats_ptr->min_dist_lsh = min_dist_lsh;
    }
}


void curve::set_stat_max_dist_lsh(double &max_dist_lsh)
{
    if ( q_stats_ptr != NULL)
    {
        q_stats_ptr->max_dist_lsh = max_dist_lsh;
    }
}


void curve::set_stat_avg_dist_lsh(double &avg_dist_lsh)
{
    if ( q_stats_ptr != NULL)
    {
        q_stats_ptr->avg_dist_lsh = avg_dist_lsh;
    }
}


int curve::get_stat_found_grid_curve()
{
    return q_stats_ptr->found_grid_curve;
}


curve* curve::get_stat_lsh_nn()
{
    return q_stats_ptr->lsh_nn ;
}


curve* curve::get_stat_real_nn()
{
    return q_stats_ptr->real_nn ;
}


double curve::get_stat_lsh_dist()
{
    return q_stats_ptr->lsh_dist ;
}


double curve::get_stat_true_dist()
{
    return q_stats_ptr->true_dist ;
}


double curve::get_stat_t_lsh()
{
    return q_stats_ptr->t_lsh ;
}


double curve::get_stat_t_true()
{
    return q_stats_ptr->t_true ;
}


double curve::get_stat_t_lsh_min()
{
    return q_stats_ptr->t_lsh_min ;
}


double curve::get_stat_t_lsh_max()
{
    return q_stats_ptr->t_lsh_max ;
}


double curve::get_stat_t_lsh_avg()
{
    return q_stats_ptr->t_lsh_avg ;
}


double curve::get_stat_min_dist_lsh()
{
    return q_stats_ptr->min_dist_lsh ;
}


double curve::get_stat_max_dist_lsh()
{

    return q_stats_ptr->max_dist_lsh ;
}


double curve::get_stat_avg_dist_lsh()
{
    return q_stats_ptr->avg_dist_lsh ;
}
