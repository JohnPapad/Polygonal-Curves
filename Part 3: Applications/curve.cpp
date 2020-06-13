#include "curve.h"
#include "cluster.h"

static unsigned int static_index = 0 ;


curve::curve(vector<vector<double>> &points, vector<double> &c_point, string type,
                unsigned int id, unsigned int number_of_points, unsigned int dimension)
{
    this->original_curve_points = points;
    this->c_point = c_point;
    this->type = type;
    this->id = id;
    this->number_of_points = number_of_points;
    this->dimension = dimension;

    this->is_centroid_of_cluster = NULL;
    this->assigned_centroid = NULL;
    this->poss_assigned_centroid = NULL;
    this->centroid_flag = 0;
    this->assigned_flag = 0;
    this->dist_to_centroid = numeric_limits<double>::infinity();
    this->dist_to_poss_centroid =numeric_limits<double>::infinity();
}


curve::curve(vector<vector<double>> &points, unsigned int way_id, unsigned int number_of_points, unsigned int dimension)
{
    this->original_curve_points = points;
    this->type = to_string(way_id);
    this->id = static_index;
    static_index++;

    this->number_of_points = number_of_points;
    this->dimension = dimension;

    this->is_centroid_of_cluster = NULL;
    this->assigned_centroid = NULL;
    this->poss_assigned_centroid = NULL;
    this->centroid_flag = 0;
    this->assigned_flag = 0;
    this->dist_to_centroid = numeric_limits<double>::infinity();
    this->dist_to_poss_centroid =numeric_limits<double>::infinity();
}



void curve::print(bool print_coords)
{
    cout.precision(5);
    cout<<"----------------------------------------"<<endl;
    cout<<endl;
    cout<<"PRINTING "<<type<<endl;
    if (print_coords == 1)
    {
        cout<<"\n PRINTING POINTS"<<endl;
        cout<<"{";
        for(unsigned int point_i=0; point_i < original_curve_points.size(); point_i++)
        {
            cout<<"(";
            for(unsigned int coord_i =0; coord_i < dimension; coord_i++)
            {
                cout.precision(5);
                cout<<fixed<<original_curve_points[point_i][coord_i];
                if (dimension - 1 != coord_i)
                    cout<<", ";
            }
            if (point_i != original_curve_points.size() -1 )
            {
                cout<<"), ";
            }
            else
            {
                cout<<")";
            }
        }
        cout<<"}"<<endl;
    }

    cout<<"PRINTING C POINT"<<endl;
    cout<<"c point size: "<<c_point.size()<<endl;
    cout<<"(";
    for(unsigned int i = 0; i < dimension ; i++)
    {
        cout<<fixed<<c_point[i];
        if (dimension - 1 != i)
            cout<<", ";
    }
    cout<<")"<<endl;


    cout<<"type: "<<type<<endl;
    cout<<"id: "<<id<<endl;
    cout<<"number of points: "<<number_of_points<<endl;
    cout<<"dimension: "<<dimension<<endl;
    cout<<"----------------------------------------"<<endl;
    cout<<endl;

}



vector<vector<double>>& curve::get_original_curve_points()
{
    return original_curve_points;
}

vector<double>& curve::get_c_point()
{
    return c_point;
}



unsigned int curve::get_curve_dimension()
{
    return dimension;
}

unsigned int curve::get_number_of_points()
{
    return number_of_points;
}


string curve::get_type()
{
    return type;
}

unsigned int curve::get_id()
{
    return id;
}



bool curve::is_centroid()
{
    return centroid_flag;
}

unsigned int curve::get_index()
{
    return id;
}

bool curve::is_assigned()
{
    return assigned_flag;
}


curve* curve::get_curve_centroid()
{
    return assigned_centroid;
}


curve* curve::get_curve_poss_centroid()
{
    return poss_assigned_centroid;
}


cluster * curve::get_cluster()
{
    return is_centroid_of_cluster;
}


double curve::get_dist_to_centroid()
{
    return dist_to_centroid;
}


double curve::get_dist_to_poss_centroid()
{
    return dist_to_poss_centroid;
}


void curve::set_cluster(cluster * cl)
{
    is_centroid_of_cluster = cl ;
}


void curve::set_is_centroid(bool flag)
{
    centroid_flag = flag;
}


void curve::set_is_assigned(bool flag)
{
    assigned_flag = flag;
}


void curve::set_curve_centroid(curve * centroid)
{
    assigned_centroid = centroid;
}


void curve::set_curve_poss_centroid(curve * centroid)
{
    poss_assigned_centroid = centroid;
}


void curve::set_dist_to_centroid(double dist)
{
    dist_to_centroid = dist ;
}


void curve::set_dist_to_poss_centroid(double dist)
{
    dist_to_poss_centroid = dist;
}

