#include "curve.h"
#include "cluster.h"

static unsigned int static_index = 0 ;

curve::curve(vector<double> &original_curve, vector<vector<double>> &original_curve_points, vector<double> &grid_curve, string curve_name,
             unsigned int number_of_points, unsigned int dimension)
{
    this->original_curve = original_curve;
    this->original_curve_points = original_curve_points;
    this->grid_curve = grid_curve;
    this->dimension = dimension;
    this->curve_name = curve_name;
    this->number_of_points = number_of_points;

    this->is_centroid_of_cluster = NULL;
    this->assigned_centroid = NULL;
    this->poss_assigned_centroid = NULL;
    this->centroid_flag = 0;
    this->index = static_index; //unique index
    static_index++;
    this->assigned_flag = 0;
    this->dist_to_centroid = numeric_limits<double>::infinity();
    this->dist_to_poss_centroid =numeric_limits<double>::infinity();
}


curve::~curve() // if q_curve_stats instance has been dynamically created for query curve we have to delete it manually
{

}


curve::curve(curve * curve_ptr, unsigned int index)
{
    stringstream ss;
    ss << static_index;
    static_index++;
    string name = ss.str();
    this->curve_name = "CURVE_ID:" + name; //unique name. the <!> is used to indicate that this curve is a mean curve (a possible centroid)

    this->index = index;
    this->original_curve = curve_ptr->get_original_curve();
    this->original_curve_points = curve_ptr->get_original_curve_points();
    this->grid_curve = curve_ptr->get_grid_curve();
    this->number_of_points = curve_ptr->get_number_of_points();
    this->dimension = curve_ptr->get_curve_dimension();
    this->centroid_flag = curve_ptr->is_centroid();
    this->assigned_flag = curve_ptr->is_assigned();
    this->is_centroid_of_cluster = curve_ptr->get_cluster();
    this->assigned_centroid = curve_ptr->get_curve_centroid();
    this->poss_assigned_centroid = curve_ptr->get_curve_poss_centroid();
    this->dist_to_centroid = curve_ptr->get_dist_to_centroid();
    this->dist_to_poss_centroid = curve_ptr->get_dist_to_poss_centroid();

}

curve::curve(vector<vector<double>> &original_curve_points, unsigned int number_of_points, unsigned int dimension)
{
    stringstream ss;
    ss << static_index;
    string name = ss.str();
    this->curve_name = "!CURVE_ID:" + name; //unique name. the <!> is used to indicate that this curve is a mean curve (a possible centroid)
    this->index = 999999;  //the index will be set correctly later
    static_index++;

    this->original_curve_points = original_curve_points;
    this->number_of_points = number_of_points;
    this->dimension = dimension;
}


void curve::print_curve(bool print_coords)
{
    typedef numeric_limits< double > dbl;
    cout.precision(dbl::max_digits10);
    cout<<"----------------------------------------"<<endl;
    cout<<endl;
    cout<<"PRINTING ORIGINAL CURVE"<<endl;
    cout<<"original curve size="<<original_curve.size()<<endl;
    if (print_coords == 1)
    {
        cout<<"(";
        for(unsigned int i=0; i < original_curve.size();i++)
        {
            cout<<fixed<<original_curve[i];
            if (original_curve.size()-1 != i)
                cout<<", ";
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
                cout<<", ";
        }
        cout<<")"<<endl;

        cout<<"\n PRINTING CURVE'S POINTS"<<endl;
        cout<<"{";
        for(unsigned int point_i=0; point_i < original_curve_points.size(); point_i++)
        {
            cout<<"(";
            for(unsigned int coord_i =0; coord_i < dimension; coord_i++)
            {
                cout.precision(dbl::max_digits10);
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

    cout<<"curve id:"<<curve_name<<endl;
    cout<<"curve index:"<<index<<endl;
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

vector<vector<double>>& curve::get_original_curve_points()
{
    return original_curve_points;
}


unsigned int curve::get_curve_dimension()
{
    return dimension;
}

unsigned int curve::get_number_of_points()
{
    return number_of_points;
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


bool curve::is_centroid()
{
    return centroid_flag;
}

unsigned int curve::get_index()
{
    return index;
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


void curve::set_curve_name(string name)
{
    curve_name = name;
}

void curve::set_original_curve(vector<double> &original_curve)
{
    this->original_curve = original_curve;
}

void curve::set_grid_curve(vector<double> &grid_curve)
{
    this->grid_curve = grid_curve;
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

void curve::set_index(unsigned int index)
{
    this->index = index;
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
