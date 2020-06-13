#include "lsh.h"

void LSH(vector< int >& r_values, multi_hash_table* MHT_ptr, curve * curve_ptr, unsigned int &L,
         unsigned int &k, unsigned int &HT_SIZE, double &delta)
{
    for(unsigned i=0; i < L ; i++) //create L hash tables and store curve to a different one each time
    {
        unsigned int HT_key = get_HT_key(r_values, curve_ptr, k, delta, HT_SIZE);
        MHT_ptr->hash_table_id[i]->store_curve_to_HT(curve_ptr, HT_key);
    }
}

vector<double> create_shifted_curve(vector<double>& curve, double &delta)
{  //shift a grid curve by adding to each of its coordianates a uniformly created random double t
    vector<double> shifted_curve = curve; //initialize shifted grid curve as the original grid curve so that a random double t to be added to all its coordinates
    double M = 0.00;
    double t = get_uniformly_distributed_random_double(M, delta); // at [0.0 , delta]
    transform(shifted_curve.begin(), shifted_curve.end(), shifted_curve.begin(),bind2nd(plus<double>(), t)); //add double number t to all grid curve's coordinates
    return shifted_curve;
}


unsigned long int get_concat_grid_curve_HV(vector<double> &concat_grid_curve, vector<int>& r_values, unsigned long int &M )
{ // hashing the concatenated grid curve( which has been created by concatenating all k grid curves) and getting an int

    double hash_value = 0.00;

    if (concat_grid_curve.size() > r_values.size()) //if there are not enough random r's for the linear combination
    {
        unsigned int limit = concat_grid_curve.size() - r_values.size();
        int r_lower_limit = 2 * concat_grid_curve.size() * -1;
        int r_upper_limit = 2 * concat_grid_curve.size();

        for(unsigned int i=0; i < limit ; i++) //create some more
        {
            int r = get_uniformly_distributed_random_int(r_lower_limit, r_upper_limit);
            r_values.push_back(r);  // and store them to the r_vector for future use
        }
    }

    for(unsigned int i=0; i < concat_grid_curve.size() ; i++)   //calculate the linear combination of the concatenated grid curve's coordinates and r_values
    {
        hash_value += concat_grid_curve[i] * r_values[i];
    }

    if (hash_value < 0.00 ) //converting  it to a positive double
    {
        hash_value *= -1;
    }

    unsigned long int rounded_hash_value= round(hash_value); //converting it to an int by rounding
    unsigned long int hash_value_mod = rounded_hash_value % M;
    return hash_value_mod;
}


unsigned int get_HT_key(vector<int>& r_values, curve * curve_ptr, unsigned int &k, double &delta, unsigned int &HT_SIZE)
{   // implementing lsh of curves , concatenating all k grid curves and getting the concatenated grid curve's hash value
    unsigned long int M = pow(2, 32) - 5; //a large first number
    vector<double> concat_grid_curve;

    for(unsigned int i = 0; i < k ;i++) //lsh of curves - creating k shifted grid curves
    {                                   //by adding a uniformly distributed random t each time to the original curve
        vector<double> shifted_curve = create_shifted_curve(curve_ptr->get_grid_curve() , delta);
        concatenating_grids_curves(shifted_curve, concat_grid_curve); // concatenating each shifted grid curve to the others
        //concat_grid_curve = remove_consecutive_duplicate_points(concat_grid_curve, curve_ptr->get_curve_dimension()); // removing consecutive duplicate points from the concatenated grid curve
    }

    unsigned long int hash_value = get_concat_grid_curve_HV(concat_grid_curve, r_values, M); // get concatenated grid curve's hash value
    unsigned int hash_table_key = hash_value % HT_SIZE; //and fit it to the hash table size
    return hash_table_key;
}


void concatenating_grids_curves(vector<double>& grid_curve, vector<double>& concat_grid_curve)
{ //storing grid_curve at the end of concat_grid_curve
    concat_grid_curve.insert( concat_grid_curve.end(), grid_curve.begin(), grid_curve.end() );
}

double coord_to_grid_coord_transform(double &delta, double &coord)
{ // match each original coordinate to a grid one
    return round(coord / delta) * delta;
}


vector<double> curve_to_grid_curve_transform(vector <double>& curve, unsigned int d, double &delta)
{ //creating a grid curve from an original curve and also removing consecutive duplicate points
    vector<double> grid_curve;
    vector<double> temp_point1;
    vector<double> temp_point2;
    double grid_coord ;

    for(unsigned int i=0; i < d; i++) //for each fist point's coordinate
    {
        grid_coord = coord_to_grid_coord_transform(delta, curve[i]);
        grid_curve.push_back(grid_coord); //add first point's grid coordinate to grid curve
        temp_point1.push_back(grid_coord); //and also keep it for future comparison for identical points
    }

    for(unsigned int i=d; i < curve.size(); i++) //for all the rest coordinates
    {
        grid_coord = coord_to_grid_coord_transform(delta, curve[i]);
        temp_point2.push_back(grid_coord); //forming next grid curve point so that we can compare it with the previously saved one

        if(temp_point2.size() == d) //if next grid curve point has been formed
        {
            if (temp_point1 != temp_point2) // and it is not identical to the previous one
            {
                for(unsigned int i=0; i < temp_point2.size(); i++)
                {
                    grid_curve.push_back(temp_point2[i]); //add it to grid curve
                }
                temp_point1.clear();
                temp_point1 = temp_point2;
                temp_point2.clear();
            }
            else //there are two consecutive duplicate grid curve points
            {
                temp_point2.clear(); // so discard the second one and don't add it to the grid curve
            }
        }
    }
    return grid_curve;
}
