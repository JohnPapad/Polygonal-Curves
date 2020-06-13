#include "lsh.h"


int main(int argc, char *argv[])
{
      double delta ;
      unsigned int HT_SIZE;
      unsigned int repeats;
      unsigned int w;

      bool flag_stats ;
      unsigned int L;
      unsigned int k;
      string input_file;
      string output_file;
      string query_file;
      string function;
      string hash;

      InputParser input(argc, argv);

      if(input.cmdOptionExists("-d"))
      {
        const string &arg_str = input.getCmdOption("-d");
        input_file = arg_str;
      }
      else
      {
          cout<<"->Please enter input file name:"<<endl;
          cin >> input_file;
      }

      if(input.cmdOptionExists("-o"))
      {
          const string &arg_str = input.getCmdOption("-o");
          output_file = arg_str;
      }
      else
      {
          cout<<"->Please enter output file name:"<<endl;
          cin >> output_file;
      }

      if(input.cmdOptionExists("-q"))
      {
          const string &arg_str = input.getCmdOption("-q");
          query_file = arg_str;
      }
      else
      {
          cout<<"->Please enter query file name:"<<endl;
          cin >> query_file;
      }

      if(input.cmdOptionExists("-k"))
      {
          const string &arg_str = input.getCmdOption("-k");
          string temp_str = arg_str;
          k = stoi(temp_str);
      }
      else
      {
          k = 2 ; //default value
      }

      if(input.cmdOptionExists("-L"))
      {
          const string &arg_str = input.getCmdOption("-L");
          string temp_str = arg_str;
          L = stoi(temp_str);
      }
      else
      {
          L = 3;  // default value
      }

      if(input.cmdOptionExists("-stats"))
      {
          flag_stats = 1;
          cout<<"->Please enter the number of repetitions"<<endl;
          cin >> repeats ;
      }
      else
      {
          flag_stats = 0;
      }

      if(input.cmdOptionExists("-function"))
      {
          const string &arg_str = input.getCmdOption("-function");
          function = arg_str;
          check_argv_string(function, "DFT", "-function");
      }
      else
      {
          cout<<"->Please enter distance function:"<<endl;
          cin>> function;
          check_argv_string(function, "DFT", "-function");
      }

      if(input.cmdOptionExists("-hash"))
      {
          const string &arg_str = input.getCmdOption("-hash");
          hash = arg_str;
          check_argv_string(hash, "classic", "-hash");
      }
      else
      {
          cout<<"->Please enter hash function:"<<endl;
          cin>> hash;
          check_argv_string(hash, "classic", "-hash");
      }

      cout<<"->Please enter delta:"<<endl;
      cin >> delta;

      cout<<"->Hash Table Size = number of curves / w .Please enter w:"<<endl;
      cin >> w;

      cout<<"----program's variables----"<<endl;
      cout<<"-input file: "<<input_file<<endl;
      cout<<"-query file: "<<query_file<<endl;
      cout<<"-k: "<<k<<endl;
      cout<<"-L: "<<L<<endl;
      cout<<"-output file: "<<output_file<<endl;
      if (flag_stats == 1)
      {
          cout<<"-stats: Yes"<<endl;
      }
      else
      {
          cout<<"-stats: No"<<endl;
      }
      cout<<"-function: DFT"<<endl;
      cout<<"-hash: classic"<<endl;
      cout<<"-delta: "<<delta<<endl;
      if (flag_stats == 1)
      {
          cout<<"-repetitions: "<<repeats<<endl;
      }
      cout<<"-w: "<<w<<endl;

    initialize_random_seed();

    vector<curve * > curves;
    parsing_saving_curves(input_file, curves, delta, 1);

    HT_SIZE = curves.size() / w ;
    if (HT_SIZE == 0)
    {
        HT_SIZE = 1;
    }
    cout<<"-hash table size: "<<HT_SIZE<<endl;;
    cout<<"==================================="<<endl;

    if ( flag_stats == 0)
    {
        run_without_stats(curves, query_file, output_file, L, k, HT_SIZE, delta);
    }
    else
    {
        run_with_stats(curves, query_file, output_file, L, k, HT_SIZE, repeats, delta);
    }

    return 0;
}


void check_argv_string(string &argv, string proper_argv, string flag_option)
{   //checking for option's -hash and -function correctness
    if ((flag_option == "-hash") && (argv == "probabilistic"))
    {
        cout<<"You chose 'probabilistic' as hash function but this is not a team project, " <<endl;
        cout<<"and therefore it is not supported.Hash function set to 'classic' by default "<<endl;
        cout<<endl;
        return;
    }

    if ((flag_option == "-function") && (argv == "DTW"))
    {
        cout<<"You chose 'DTW' as distance function but this is not a team project, " <<endl;
        cout<<"and therefore it is not supported.Distance function set to 'DFT' by default "<<endl;
        cout<<endl;
        return;
    }

    if (argv == " ")
    {
        cout<<"You forgot to enter a value for option:"<<flag_option<<endl;
        cout<<"Option: "<<flag_option<<" set to '"<<proper_argv<<"' by default"<<endl;
        cout<<endl;
        argv = proper_argv;
    }
    else if (argv != proper_argv)
    {
        cout<<"You entered false input for option:"<<flag_option<<endl;
        cout<<"Option: "<<flag_option<<" set to '"<<proper_argv<<"' by default"<<endl;
        cout<<endl;
        argv = proper_argv ;
    }

    return;
}


void run_without_stats(vector<curve *> &curves, string &q_data_set_file_name, string &output_file_name,
                       unsigned int &L, unsigned int &k, unsigned int &HT_SIZE, double &delta)
{   //when -stats option was not given
    vector<int> r_values; //a vector that stores uniformly created numbers
    multi_hash_table* MHT_ptr = create_multi_hash_table(L, HT_SIZE); //creating the L hash tables

    for (unsigned int i=0; i<curves.size() ; i++) //for each dataset's curve
    {
        LSH(r_values, MHT_ptr, curves[i], L, k, HT_SIZE, delta);  // hash curve and store it to L hash tables
    }
    //hash tables and curves' pointers stored to them will not be destroyed to each possible repetition
    //the previous steps will only be executed once

    do
    {
        vector<curve * > queries;  // queries will be destroyed each time the process is repeated
        ofstream outfile (output_file_name, ios_base::app | ios_base::out);  //so that we can append info to already created files

        parsing_saving_curves(q_data_set_file_name, queries, delta, 1); //read and save queries

        for (unsigned int i=0; i < queries.size() ; i++)  // for each query
        {
            queries[i]->create_q_stats(); // create a class member struct called q_curve_stats to store all statistics
            get_lsh_nearest_neighbor(r_values, MHT_ptr, queries[i], k, delta, HT_SIZE);
            get_real_nearest_neighbor(curves, queries[i]);
            write_results_in_file(outfile, queries[i]);
        }

        delete_curves(queries);  //taking care of memory leaks (for each repetition)

        if (repeat(q_data_set_file_name, output_file_name) == 0)  //exit
        {
            outfile.close();
            break;
        }

    } while (1);

    //deleting all dynamically created structures
    destroy_multi_hash_table(MHT_ptr);
    delete_curves(curves);
}


bool repeat(string &q_data_set_file_name, string &output_file_name)
{
    string answer;
    do
    {
        cout<<"->Would you like to repeat the process with another searching data set ?"<<endl;
        cout<<"->If you would like please press 'yes' , if not press 'no' to exit."<<endl;
        cin >> answer ;
    } while ((answer != "yes")&&(answer != "no")); // checking for correctness

    if(answer == "no")  //exit
    {
        return 0;
    }
    else
    {
        cout<<"->Please type searching data set's file name:"<<endl;
        cin >> q_data_set_file_name;
        cout<<"->Please type output file's name:"<<endl;
        cin >> output_file_name;
        return 1;
    }
}


void run_with_stats(vector<curve *> &curves, string &q_data_set_file_name, string &output_file_name,
                       unsigned int &L, unsigned int &k, unsigned int &HT_SIZE, unsigned int &repeats, double &delta)
{   //when -stats option was given
    do
    {
        ofstream outfile (output_file_name, ios_base::app | ios_base::out);
        vector<curve * > queries;
        parsing_saving_curves(q_data_set_file_name, queries, delta, 1); //read and save queries

        for (unsigned int i=0; i < queries.size() ; i++)  // for each query
        {

            double minDistanceLSH = numeric_limits<double>::infinity();
            double maxDistanceLSH = -1.0 ;
            double sumDistanceLSH =  0.0 ;
            double tLSHmin = numeric_limits<double>::infinity();
            double tLSHmax = -1.0 ;
            double tLSHsum =  0.0 ;

            queries[i]->create_q_stats(); // create a class member struct called q_curve_stats to store all statistics
            get_real_nearest_neighbor(curves, queries[i]); // calculate min real dist and real Nearest neighbor only once for each curve

            for(unsigned int j =0 ; j < repeats ; j++)  //repeat process 100 times for each curve to get statistics
            {
                vector<int> r_values; //a vector that stores uniformly created numbers
                multi_hash_table* MHT_ptr = create_multi_hash_table(L, HT_SIZE); //creating the L hash tables
                                                                                 //each time will be destroyed and created again
                for (unsigned int k=0; k<curves.size() ; k++) //for each dataset's curve
                {
                    LSH(r_values, MHT_ptr, curves[k], L, k, HT_SIZE, delta);  // hash curve and store it to L hash tables
                }

                get_lsh_nearest_neighbor(r_values, MHT_ptr, queries[i], k, delta, HT_SIZE);

                if (queries[i]->get_stat_lsh_dist() < minDistanceLSH)
                {
                    minDistanceLSH = queries[i]->get_stat_lsh_dist();
                }

                if (queries[i]->get_stat_lsh_dist() > maxDistanceLSH)
                {
                    maxDistanceLSH = queries[i]->get_stat_lsh_dist();
                }

                sumDistanceLSH += queries[i]->get_stat_lsh_dist();

                if (queries[i]->get_stat_t_lsh() < tLSHmin)
                {
                    tLSHmin = queries[i]->get_stat_t_lsh();
                }

                if (queries[i]->get_stat_t_lsh() > tLSHmax)
                {
                    tLSHmax = queries[i]->get_stat_t_lsh();
                }

                tLSHsum += queries[i]->get_stat_t_lsh();

                destroy_multi_hash_table(MHT_ptr); //destroy the L hash tables in order to be recreated at the next repetition
            }

            //save all statistics to query
            queries[i]->set_stat_min_dist_lsh(minDistanceLSH);
            queries[i]->set_stat_max_dist_lsh(maxDistanceLSH);
            double avgDistanceLSH = sumDistanceLSH / repeats ;
            queries[i]->set_stat_avg_dist_lsh(avgDistanceLSH);

            queries[i]->set_stat_t_lsh_min(tLSHmin);
            queries[i]->set_stat_t_lsh_max(tLSHmax);
            double tLSHavg = tLSHsum / repeats ;
            queries[i]->set_stat_t_lsh_avg(tLSHavg);

            write_stats_in_file(outfile, queries[i]);
        }//finished with all queries

        delete_curves(queries); //taking care of memory leaks

        if (repeat(q_data_set_file_name, output_file_name) == 0)  //exit
        {
            outfile.close();
            break;
        }

    }while(1);

    delete_curves(curves);
}


void write_results_in_file(ofstream &outfile, curve * query_ptr)
{
    outfile.precision(8);
    outfile << "Query: " << query_ptr->get_curve_name()<< endl;
    outfile << "DistanceFunction: DFT"<< endl;
    outfile << "HashFunction: Classic" << endl;
    outfile << "FoundGridCurve: ";
    if (query_ptr->get_stat_found_grid_curve() == 1)
    {
        outfile<<"True"<< endl;
    }
    else
    {
        outfile<<"False"<< endl;
    }
    outfile<< "LSH Nearest neighbor: "<< query_ptr->get_stat_lsh_nn()->get_curve_name()<<endl;
    outfile<< "True Nearest neighbor: "<< query_ptr->get_stat_real_nn()->get_curve_name()<<endl;
    outfile<< "distanceLSH: "<< fixed<<query_ptr->get_stat_lsh_dist()<<endl;
    outfile<< "distanceTrue: "<< fixed<<query_ptr->get_stat_true_dist()<<endl;
    outfile<< endl;
    outfile<<"========================================================="<<endl;
}


void write_stats_in_file(ofstream &outfile, curve * query_ptr) //when option stats has been choosen
{
    outfile.precision(8);
    outfile << "Query: " << query_ptr->get_curve_name()<< endl;
    outfile << "DistanceFunction: DFT"<< endl;
    outfile << "HashFunction: Classic" << endl;
    outfile << "|minDistanceLSH – distanceTrue|: "<<fixed<< fabs(query_ptr->get_stat_min_dist_lsh() - query_ptr->get_stat_true_dist()) <<endl;
    outfile << "|maxDistanceLSH – distanceTrue|: "<<fixed<< fabs(query_ptr->get_stat_max_dist_lsh() - query_ptr->get_stat_true_dist()) <<endl;
    outfile << "|avgDistanceLSH – distanceTrue|: "<<fixed<< fabs(query_ptr->get_stat_avg_dist_lsh() - query_ptr->get_stat_true_dist()) <<endl;
    outfile << "tLSHmin: "<<fixed<<query_ptr->get_stat_t_lsh_min()<<endl;
    outfile << "tLSHmax: "<<fixed<<query_ptr->get_stat_t_lsh_max()<<endl;
    outfile << "tLSHavg: "<<fixed<<query_ptr->get_stat_t_lsh_avg()<<endl;
    outfile << "tTrue: "<<fixed<<query_ptr->get_stat_t_true()<<endl;
    outfile<< endl;
    outfile<<"========================================================="<<endl;
}


void parsing_saving_curves(string &file_name, vector<curve *> &curves, double & delta, bool q_curve_flag)
{
    ifstream file(file_name);

    if (!file) //checking if file exists
    {
        cout<<"-ERROR-There is no input file: '"<<file_name<<"'"<<endl;
        exit(1);
    }

    string line;

    while(getline(file, line)) //get a dataset line
    {
        if (q_curve_flag == 1) //if function being used for parsing queries bypass first line (which defines R , as this is not a team project, so the R is always 0)
        {
            q_curve_flag = 0;
            continue;
        }

        vector <double> original_curve;
        istringstream iss(line);
        string curve_name; //curve's id
        unsigned int points; //number of curve's points
        iss>>curve_name;
        iss>>points;
        string curve_info; //remove the first 2 numbers of the line (curve_id and number of curve's points)
        getline(iss, curve_info, '\t');
        string points_str;
        getline(iss, points_str, '\t'); //get all curve's points
        unsigned int curve_dimension = get_curve_dimension(points_str);
        istringstream iss2(points_str);
        string coord_str;

        while(getline(iss2, coord_str, ' ')) //get all point's coordinates
        {
            //seperating the actual double coordinate from the useless characters
            if ((coord_str.at(0) == '(' ) || (coord_str.at(0) == ',' ))
            {
                coord_str.erase(coord_str.begin());
            }

            if ((coord_str.back() == ')') || (coord_str.back() == ','))
            {
                coord_str.pop_back();
            }

            istringstream iss3(coord_str);
            double coord;
            iss3 >> coord;
            original_curve.push_back(coord); // save each coord to a vector of doubles
        }
        vector<double> grid_curve = curve_to_grid_curve_transform(original_curve, curve_dimension, delta); //creating the grid curve
        curve * curve_ptr = new curve(original_curve, grid_curve, curve_name, points, curve_dimension); // save all curve's info to a class
        curves.push_back(curve_ptr); //add new class curve's pointer to a vector
    }
}


unsigned int get_curve_dimension(string points_str)
{
    unsigned int d = 1;
    for (unsigned int i=0; i < points_str.size(); i++) // basically counting how many commas are between one open and one close parenthesis
    {
        if (points_str.at(i) == ')')
        {
            break;
        }
        else if (points_str.at(i) == ',')
        {
            d++;
        }
    }
    return d;
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


double get_uniformly_distributed_random_double(double M, double N)
{  //pseudocode taken from eclass notes
    double result;
    do
    {
        result = M + (rand() / (RAND_MAX + 1.0))*(N-M);
    }while(result == N);
    return result;
}


int get_uniformly_distributed_random_int(int M, int N)
{  ///pseudocode taken from eclass notes
    int result = M + (rand() / (RAND_MAX + 1.0))*(N-M+1);
    return result;
}


vector<double> create_shifted_curve(vector<double>& curve, double &delta)
{  //shift a grid curve by adding to each of its coordinates a uniformly created random double t
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
        concat_grid_curve = remove_consecutive_duplicate_points(concat_grid_curve, curve_ptr->get_curve_dimension()); // removing consecutive duplicate points from the concatenated grid curve
    }

    unsigned long int hash_value = get_concat_grid_curve_HV(concat_grid_curve, r_values, M); // get concatenated grid curve's hash value
    unsigned int hash_table_key = hash_value % HT_SIZE; //and fit it to the hash table size
    return hash_table_key;
}


void concatenating_grids_curves(vector<double>& grid_curve, vector<double>& concat_grid_curve)
{ //storing grid_curve at the end of concat_grid_curve
    concat_grid_curve.insert( concat_grid_curve.end(), grid_curve.begin(), grid_curve.end() );
}


vector <double>  remove_consecutive_duplicate_points(vector <double> dubl_curve, unsigned int d)
{
    vector<double> unq_curve;
    vector<double> temp_point1;
    vector<double> temp_point2;

    for(unsigned int i=0; i < d; i++)
    {
        unq_curve.push_back(dubl_curve[i]); //add first point's coordinates to new curve
        temp_point1.push_back(dubl_curve[i]); //and also keep them for future comparison for identical points
    }

    for(unsigned int i=d; i < dubl_curve.size(); i++)
    {
        temp_point2.push_back(dubl_curve[i]); //forming next curve point
        if(temp_point2.size() == d) //if next curve point has been formed
        {
            if (temp_point1 != temp_point2) // and it is not identical to the previous one
            {
                for(unsigned int i=0; i < temp_point2.size(); i++)
                {
                    unq_curve.push_back(temp_point2[i]); //add it to the new curve
                }
                temp_point1.clear();
                temp_point1 = temp_point2;
                temp_point2.clear();
            }
            else //there are two consecutive duplicate curve points
            {
                temp_point2.clear(); // so discard the second one and don't add it to the new curve
            }
        }
    }
    return unq_curve;
}


void remove_dublicate_curves(vector<curve *> &matched_curves)
{
    for(unsigned int i=0; i < matched_curves.size(); i++)
    {
        for(unsigned int j=i+1; j < matched_curves.size(); j++)
        {
            if (matched_curves[i]->get_curve_name() == matched_curves[j]->get_curve_name())  //searching by curves' names (ids) because the vector contains pointers to curves
            {
                matched_curves.erase(matched_curves.begin() + i);  //removing duplicate
                j--;
            }
        }
    }
}


void LSH(vector< int >& r_values, multi_hash_table* MHT_ptr, curve * curve_ptr, unsigned int &L,
         unsigned int &k, unsigned int &HT_SIZE, double &delta)
{
    for(unsigned i=0; i < L ; i++) //create L hash tables and store curve to a different one each time
    {
        unsigned int HT_key = get_HT_key(r_values, curve_ptr, k, delta, HT_SIZE);
        MHT_ptr->hash_table_id[i]->store_curve_to_HT(curve_ptr, HT_key);
    }
}


void get_real_nearest_neighbor(vector<curve *> &curves, curve * query_ptr)
{
    double min_DFD;
    double elapsed_time;
    unsigned int min_index ; // an index of matched_curves vector that indicates the curve(basically its pointer)
                             // that has the minimum frechet distance from the query curve
    min_DFD = get_min_DFD(curves, query_ptr, min_index, elapsed_time);

    // all stats are saved to a query curve's struct member called q_curve_stats
    query_ptr->set_stat_real_nn(curves[min_index]);
    query_ptr->set_stat_true_dist(min_DFD);
    query_ptr->set_stat_t_true(elapsed_time);
}

void get_lsh_nearest_neighbor(vector<int>& r_values, multi_hash_table* MHT_ptr,
                          curve * query_ptr, unsigned int &k, double &delta, unsigned int &HT_SIZE)
{
    int found_grid_curve;
    unsigned int HT_key = get_HT_key(r_values, query_ptr, k, delta, HT_SIZE); //get a hash value for query curve
    vector<curve * > matched_curves; // a vector containing pointers to curves that match with query curve(have the same grid curve)
    MHT_ptr->searching_MHT(query_ptr, matched_curves, HT_key, found_grid_curve); //searching all hash tables for matches
    query_ptr->set_stat_found_grid_curve(found_grid_curve);
    remove_dublicate_curves(matched_curves);

    double min_DFD;
    double elapsed_time;
    unsigned int min_index ; // an index of matched_curves vector that indicates the curve(basically its pointer)
                             // that has the minimum frechet distance from the query curve
    if (matched_curves.size() != 0)
    {
        min_DFD = get_min_DFD(matched_curves, query_ptr, min_index, elapsed_time);
        // all stats are saved to a query curve's struct member called q_curve_stats
        query_ptr->set_stat_lsh_nn(matched_curves[min_index]);
        query_ptr->set_stat_lsh_dist(min_DFD);
        query_ptr->set_stat_t_lsh(elapsed_time);
    }
}


double get_min_DFD(vector<curve * > &curves, curve * query_ptr, unsigned int &min_index, double &elapsed_time)
{   //get the minimum frechet distance between a query curve and a vector of curves
    double min_DFD = numeric_limits<double>::infinity();
    clock_t begin = clock();

    for(unsigned int i = 0 ; i < curves.size() ; i++)
    {
        double DFD = frechet_distance(query_ptr->get_original_curve(), curves[i]->get_original_curve(), query_ptr->get_curve_dimension());

        if (DFD < min_DFD)
        {
            min_DFD = DFD;
            min_index = i ;
        }

        if (min_DFD == 0.0)  // we have found the closest neighbor
        {
            break;
        }
     }

     clock_t end = clock();
     elapsed_time = double(end - begin) / CLOCKS_PER_SEC;
     return min_DFD;
}


double frechet_distance(vector<double> &curve1, vector<double> &curve2, unsigned int d)
{  // implementing the pseudocode
    unsigned int n = curve1.size() / d; // n,m = number of points
    unsigned int m = curve2.size() / d;

    double L[n][m] ;

    for(unsigned int i=0; i < n ; i++)
    {
        for(unsigned int j=0; j < m ; j++)
        {
            if ((i == 0) && (j == 0))
            {
                vector<double> p = get_a_point_from_curve(curve1, 0, d);
                vector<double> q = get_a_point_from_curve(curve2, 0, d);
                L[i][j] = eucl_distance(p, q);
            }
            else if (i == 0)
            {
                vector<double> p = get_a_point_from_curve(curve1, 0, d);
                vector<double> q = get_a_point_from_curve(curve2, j, d);
                L[i][j] = max( eucl_distance(p, q) , L[0][j - 1] );
            }
            else if (j == 0)
            {
                vector<double> p = get_a_point_from_curve(curve1, i, d);
                vector<double> q = get_a_point_from_curve(curve2, 0, d);
                L[i][j] = max( eucl_distance(p, q) , L[i - 1][0] );
            }
            else
            {
                vector<double> p = get_a_point_from_curve(curve1, i, d);
                vector<double> q = get_a_point_from_curve(curve2, j, d);
                double min_L = min( L[i-1][j], min( L[i][j-1], L[i-1][j-1] ));
                L[i][j] = max( eucl_distance(p, q) , min_L);
                //L[i][j] = max( eucl_distance(p, q) , min( L[i-1][j], L[i][j-1], L[i-1][j-1] ));
            }
        }
    }
    return L[n-1][m-1] ;
}


double eucl_distance(vector<double> &point1, vector<double> &point2)
{  //calculate the euclidean distance between 2 points of any dimension
    double sum = 0.0;

    for(unsigned int i=0 ; i < point1.size(); i++) //for each points' coordinate (point1 and point2 have the same dimension)
    {
        double d_coor = point1[i] - point2[i] ;
        sum += pow(d_coor, 2.0);
    }

    double result = sqrt(sum);
    return result;
}


vector<double> get_a_point_from_curve(vector<double> &curve, unsigned int point_index, unsigned int d)
{  // forming a point out of curve's coordinates according to point's index and curve's dimension
    vector<double> point;
    unsigned int l_limit = point_index * d;
    unsigned int u_limit = l_limit + d ;

    for(unsigned int i= l_limit; i < u_limit ; i++)
    {
        point.push_back(curve[i]);
    }

    return point;
}


multi_hash_table* create_multi_hash_table(unsigned int &L, unsigned int &HT_SIZE)
{  // creating L hash tables (basically a multi hash table) and returning a pointer to them
    multi_hash_table* MHT_ptr = new multi_hash_table(L, HT_SIZE);
    return MHT_ptr;
}


void destroy_multi_hash_table(multi_hash_table* MHT_ptr)
{
    delete MHT_ptr;
    MHT_ptr = NULL;
}


void delete_curves(vector<curve * > &curves)
{
    for (unsigned int i=0; i < curves.size(); i++)
    {
        delete curves[i];
    }
}


void initialize_random_seed() //for rand function
{
    srand (time(NULL));
}


void print_all_curves(vector<curve *> &curves)
{
    for (unsigned int i=0; i < curves.size() ; i++)
    {
        curves[i]->print_curve();
    }
}


/*
hash_table* create_hash_table(unsigned int &HT_SIZE)
{
    hash_table* HT_ptr = new hash_table(HT_SIZE);
    return HT_ptr;
}

void destroy_hash_table(hash_table * HT_ptr)
{
    delete HT_ptr;
}

bool compare_curves(vector<double>& curve1, vector<double>& curve2)
{
    if (curve1.size() != curve2.size() )
    {
        return 0;
    }

    for(unsigned int i=0; i < curve1.size(); i++)
    {
        if ( curve1[i] != curve2[i])
        {
            return 0;
        }
    }
    return 1;
}
*/
