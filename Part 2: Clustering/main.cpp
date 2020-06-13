#include "lsh.h"
#include "bst.h"

int main(int argc, char *argv[])
{
    InputParser input(argc, argv);

    string input_file;
    string config_file;
    string output_file;
    bool complete_flag;

    if(input.cmdOptionExists("-i"))
    {
        const string &arg_str = input.getCmdOption("-i");
        input_file = arg_str;
    }
    else
    {
        cout<<"->Please enter input file name:"<<endl;
        cin >> input_file;
    }

    if(input.cmdOptionExists("-c"))
    {
          const string &arg_str = input.getCmdOption("-c");
          config_file = arg_str;
    }
    else
    {
        cout<<"->Please enter configuration file name:"<<endl;
        cin >> config_file;
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

    if(input.cmdOptionExists("-d"))
    {
        const string &arg_str = input.getCmdOption("-d");
        string dist_name = arg_str;
        if (dist_name != "DFD")
        {
            if(dist_name == "DTW")
            {
                cout<<"->You chose DTW as the distance metric but this is not a team project, so it is not supported."<<endl;
            }
            else
            {
                cout<<"-WRONG- input as a distance metric!"<<endl;
            }
            cout<<"->Distance metric set to DFD by default"<<endl;
        }
    }
    else
    {
        cout<<"->You forgot to choose distance metric. DFD was set by default."<<endl;
    }

    if(input.cmdOptionExists("-complete"))
    {
        complete_flag = 1;
    }
    else
    {
        complete_flag = 0;
    }

    unsigned int L ;
    unsigned int kappa ;
    unsigned int k ;

    vector<pair<string, unsigned int >> vars;
    vars.push_back(make_pair("number_of_clusters:", 0));
    vars.push_back(make_pair("number_of_grid_curves:", 2));
    vars.push_back(make_pair("number_of_hash_tables:", 3));

    read_config_file(config_file, vars);
    k = vars[0].second;
    kappa= vars[1].second;
    L = vars[2].second;

    cout<<"---->Configs<----"<<endl;
    cout<<"Number of clusters (k): "<<k<<endl;
    cout<<"Number of grid curves (kappa): "<<kappa<<endl;
    cout<<"Number of hash tables (L): "<<L<<endl;
    cout<<endl;
    cout<<"Input file: "<<input_file<<endl;
    cout<<"Configuration file: "<<config_file<<endl;
    cout<<"Output file: "<<output_file<<endl;
    cout<<"Distance metric: DFD"<<endl;
    cout<<"-----------------"<<endl;

    //exit(1);

    double delta = 0.02 ;
    unsigned int HT_SIZE;
    unsigned int w = 4;
    initialize_random_seed();

    vector<curve * > curves;
    parsing_saving_curves(input_file, curves, delta, 1);
    HT_SIZE = curves.size() / w ;
    cout<<"HT_SIZE: "<<HT_SIZE<<endl;
    cout<<"Number of curves: "<<curves.size()<<endl;

    vector<int> r_values;
    multi_hash_table* MHT_ptr = initialize_MHT(curves, r_values, L, kappa, HT_SIZE, delta);
    cout<<"--->Multi hash table created<---"<<endl;

    bool mean_curve_flag ;
    bool exit_flag;
    bool lsh_flag;
    bool init_flag;

    do
    {
        cout<<"->Please choose among program's different variations"<<endl;
        cout<<"->For initialization with random selection of k-points please press '0'."<<endl;
        cout<<"->For initialization with k-means++ please press '1'."<<endl;
        init_flag = read_flag();

        cout<<"->For Lloyd's assignment please press '0'. For Assignment by Range Search (LSH) please press '1;."<<endl;
        lsh_flag = read_flag();

        cout<<"->For Partitioning Around Medoids (PAM) please press '0'. For computing the Mean Discrete Frechet curve please press '1'."<<endl;
        mean_curve_flag = read_flag();

        print_program_variation(init_flag, lsh_flag, mean_curve_flag);

        vector<vector<double> > DFDs;
        if (mean_curve_flag == 0)
        {
            vector<vector<double> > temp_DFDs( curves.size(), vector<double>(curves.size(), -1));
            DFDs = temp_DFDs;
        }
        else
        {
            vector<vector<double> > temp_DFDs( curves.size() + 2*k, vector<double>(curves.size() + 2*k, -1));
            DFDs = temp_DFDs;
        }

        vector<curve *> centroids;
        if (init_flag == 0)
        {
            centroids = random_k_points_selection(curves, k);
        }
        else
        {
            centroids = k_means_pp(curves, DFDs, k);
        }
        cout<<"--->Chose k centroids<---"<<endl;

        create_clusters(centroids);
        cout<<"--->Clusters created<---"<<endl;

        if (mean_curve_flag == 1)
        {
            transform_starting_centroids(centroids, curves.size());
        }

        double bound;
        if (( lsh_flag == 1) && (mean_curve_flag == 1))
        {
            bound = 0.01;
        }
        else
        {
            bound = 0.001;
        }
        unsigned int iterations = 0;
        clock_t begin = clock();

        do
        {
            if (iterations == 4)
            {
                iterations = 0;
                bound *= 4;
            }

            if (lsh_flag == 0)
            {
                simple_clustering(curves, centroids, DFDs);
            }
            else
            {
                clustering_with_LSH(curves, centroids, DFDs, r_values, MHT_ptr, kappa, delta, HT_SIZE, mean_curve_flag);
            }

            bool end_program;
            if (mean_curve_flag == 0)
            {
                end_program = PAM(curves, centroids, DFDs, bound);
            }
            else
            {
                end_program = find_mean_centroids(curves, centroids, DFDs, bound, delta, lsh_flag);
            }
            //cout<<"--->Done updating<---"<<endl;

            if ( end_program == 1)
            {
                break;
            }
            else
            {
                reset_all_curves(curves);
                reset_all_clusters(centroids);
            }
            iterations++;
        }while(1);

        clock_t end = clock();
        double elapsed_time = double(end - begin) / CLOCKS_PER_SEC;

        vector<double> clusters_silhts;
        double overall_silht = get_silhouette(centroids, DFDs, clusters_silhts);

        write_results_in_file(centroids, clusters_silhts, overall_silht, output_file, complete_flag, elapsed_time, init_flag, lsh_flag, mean_curve_flag);
        reset_all_clusters(centroids);
        destroy_clusters(centroids, mean_curve_flag);
        reset_curves(curves);

        cout<<"->Please press '0' if you would like to run an other variation, or press 1 to exit."<<endl;
        exit_flag = read_flag();
        if (exit_flag == 1)
        {
            break;
        }

    }while(1);

    destroy_multi_hash_table(MHT_ptr);
    delete_curves(curves);

    return 0;
}


bool read_flag() //checking for false input
{
    bool flag;
    while (!(std::cin >> flag))
    {
        cout<<"-WRONG INPUT- Please press 0 or 1."<<endl;
        cin.clear();
        cin.ignore(numeric_limits<std::streamsize>::max(), '\n');
    }
    return flag;
}

void print_program_variation(bool init_flag, bool lsh_flag, bool mean_curve_flag)
{
    cout<<"->Initialization method: ";
    if (init_flag == 0)
    {
        cout<<"random selection of k-points";
    }
    else
    {
        cout<<"k-means++";
    }
    cout<<endl;

    cout<<"->Assignment method: ";
    if (lsh_flag == 0)
    {
        cout<<"Lloyd's";
    }
    else
    {
        cout<<"LSH";
    }
    cout<<endl;

    cout<<"->Updating method: ";
    if (mean_curve_flag == 0)
    {
        cout<<"PAM";
    }
    else
    {
        cout<<"Mean Discrete Frechet curve";
    }
    cout<<endl;
}

void read_config_file(string config_file, vector<pair<string, unsigned int >> &vars)//unsigned int &k, unsigned int &kappa, unsigned int &L)
{
    //string file_name = "cluster.conf" ;
    vector<pair<string, unsigned int>> configs;
    ifstream file(config_file);

    if (!file) //checking if file exists
    {
        cout<<"-ERROR 404-There is no config file: "<<endl;
        exit(1);
    }

    string line;
    while(getline(file, line)) //read line
    {
        istringstream iss(line);
        string var_name;
        unsigned int var_val;
        iss>>var_name;
        iss>>var_val;
        configs.push_back(make_pair(var_name, var_val));
    }

    if (configs.size() == 0)
    {
        cout<<"-ERROR-Config file is empty"<<endl;
        exit(-1);
    }

    if (configs[0].first != vars[0].first)
    {
        cout<<"-ERROR-Config does not contain number of clusters which is mandatory"<<endl;
        exit(-2);
    }
    else
    {
        vars[0].second = configs[0].second ;
    }


    if (configs.size() == 1)
    {
        cout<<"Variables will be set to their default values"<<endl;
        return;
    }
    else
    {
        for(unsigned int i=1; i < vars.size(); i++)
        {
            for(unsigned j=1; j < configs.size(); j++)
            {
                if (vars[i].first == configs[j].first)
                {
                    vars[i].second = configs[j].second;
                }
            }
        }
    }
}

void write_results_in_file(vector<curve *> &centroids, vector<double> &clusters_silhts, double overall_silht, string output_file_name, bool complete_flag, double clustering_time, bool init_flag, bool lsh_flag, bool mean_curve_flag)
{
    ofstream outfile (output_file_name, ios_base::app | ios_base::out);
    outfile.precision(8);
    unsigned int  i=init_flag+1;
    unsigned int  j=lsh_flag+1;

    unsigned int k;
    if (mean_curve_flag ==1)
    {
        k = 1;
    }
    else
    {
        k = 2;
    }

    outfile << "Algorithm: "  << i<<"x"<<j<<"x"<<k<<endl;
    outfile << "Metric: Frechet"<< endl<<endl;
    for(unsigned int i=0; i < centroids.size(); i++)
    {
        outfile << "CLUSTER-" << centroids[i]->get_cluster()->get_id()<<" {size: "<< centroids[i]->get_cluster()->get_cluster_size()<<", centroid: ";
        if (mean_curve_flag == 0)
        {
            outfile<<centroids[i]->get_curve_name()<<" }"<<endl;
        }
        else
        {
            outfile << endl;
            vector<vector<double>> original_curve_points = centroids[i]->get_original_curve_points();
            unsigned int dimension  = centroids[i]->get_curve_dimension();
            outfile<<"[";
            for(unsigned int point_i=0; point_i < original_curve_points.size(); point_i++)
            {
                outfile<<"(";
                for(unsigned int coord_i =0; coord_i < dimension; coord_i++)
                {
                    outfile<<fixed<<original_curve_points[point_i][coord_i];
                    if (dimension - 1 != coord_i)
                        outfile<<", ";
                }
                if (point_i != original_curve_points.size() -1 )
                {
                    outfile<<"), ";
                }
                else
                {
                    outfile<<")";
                }
            }
            outfile<<"] }"<<endl;
        }
        outfile<<"-------------------------------------------------------------------"<<endl;
    }

    outfile << "clustering_time: "<< clustering_time<< endl;
    outfile << "Silhouette: [";
    for(unsigned int i=0; i < clusters_silhts.size(); i++)
    {
        outfile <<fixed<<clusters_silhts[i]<<", ";
    }
    outfile<<fixed<<overall_silht<<"]"<<endl;

    outfile<<endl;
    outfile<<endl;
    if (complete_flag == 1)
    {
        for(unsigned int i=0; i < centroids.size(); i++)
        {
            outfile << "CLUSTER-" << centroids[i]->get_cluster()->get_id()<<" {";
            vector<curve *> cluster_members = centroids[i]->get_cluster()->get_members();
            for(unsigned int j=0; j < cluster_members.size(); j++)
            {
                outfile<<cluster_members[j]->get_curve_name();
                if ( j != cluster_members.size() -1)
                {
                    outfile<<", ";
                }
            }
            outfile<<"}"<<endl;

            outfile<<"-------------------------------------------------------------------"<<endl;
        }
    }

    outfile<<"==================================================================="<<endl;
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

void transform_starting_centroids(vector<curve *> &centroids, unsigned int curves_size)
{ // create a new curve (centroid) for each originally selected dataset's curve (centroid)
    for(unsigned int i=0; i < centroids.size(); i++)
    {
        curve * centroid = new curve(centroids[i], curves_size + i);
        reset_curve(centroids[i]); //mark originally selected centroid as a simple curve
        centroids[i] = centroid;
    }
}



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

vector<curve *> k_means_pp(vector<curve *> curves, vector<vector<double> > &DFDs, unsigned int k)
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
            double min_dfd = get_min_DFD(centroids, curves[r], DFDs, double_null, int_null, int_null);
            PartialSum += min_dfd * min_dfd ;
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

void simple_clustering(vector<curve *>& curves, vector<curve *>& centroids, vector<vector<double> > &DFDs)
{
    for(unsigned int i = 0; i < curves.size(); i++)
    { //for each curve that has not been marked as a centroid
        if ( curves[i]->is_centroid() == 0)
        {
            lloyed_assignment(centroids, curves[i], DFDs);
        }
    }
}

void lloyed_assignment(vector<curve *> &centroids, curve * curve_ptr, vector<vector<double> > &DFDs)
{
    double min_DFD;
    double sec_min_DFD;
    unsigned int centroid_index ;
    unsigned int poss_centroid_index;
    curve* centroid_ptr ;
    curve* poss_centroid_ptr;

    min_DFD = get_min_DFD(centroids, curve_ptr, DFDs, sec_min_DFD, centroid_index, poss_centroid_index);
    centroid_ptr = centroids[centroid_index];
    poss_centroid_ptr = centroids[poss_centroid_index];

    curve_ptr->set_curve_centroid(centroid_ptr); //set its closest centroid
    curve_ptr->set_dist_to_centroid(min_DFD); //set its dfd with its closest centroid

    centroid_ptr->get_cluster()->add_member(curve_ptr); //add curve to the centroid's cluster
    centroid_ptr->get_cluster()->update_cost(min_DFD); // update the cluster's sum of dfds

    curve_ptr->set_curve_poss_centroid(poss_centroid_ptr); //set its second closest centroid
    curve_ptr->set_dist_to_poss_centroid(sec_min_DFD); // set its dfd with the second closest centroid

}


void clustering_with_LSH(vector<curve *>& curves, vector<curve *>& centroids,  vector<vector<double> > &DFDs, vector<int>& r_values,
                         multi_hash_table* MHT_ptr, unsigned int &k, double &delta, unsigned int &HT_SIZE, bool mean_flag)
{
    double centroid_dfd ;
    double min_centroid_dfd = numeric_limits<double>::infinity();

    unsigned int int_null = 0;
    double double_null = 0.0 ;

    vector<curve *> temp_centroids = centroids;

    for(unsigned int i = 0; i < centroids.size() - 1; i++)
    { //calculate the min dfd difference among all centroids
        temp_centroids.erase(temp_centroids.begin());
        centroid_dfd = get_min_DFD(temp_centroids, centroids[i], DFDs, double_null, int_null, int_null);

        if (centroid_dfd < min_centroid_dfd)
        {
            min_centroid_dfd = centroid_dfd ;
        }
    }

    unsigned int no_assignment_iterations = 0;
    double range_limit = min_centroid_dfd;
    vector< pair<string, curve *>> arranged_curves ;
    do
    {
        vector< pair<string, curve *>> checked_curves ; //all curves that have been checked during a range iteration
        for(unsigned int i = 0; i < centroids.size(); i++)
        {
            assignment_range_search(checked_curves, DFDs, r_values, MHT_ptr, centroids[i], range_limit, k, delta, HT_SIZE);
        }

        for(unsigned int i = 0; i < checked_curves.size(); i++)
        {
            no_assignment_iterations = 0 ;
            checked_curves[i].second->set_is_assigned(1); //mark curve as assigned
            arranged_curves.push_back(checked_curves[i]); //keep all assigned curves to a vector
            checked_curves[i].second->get_curve_centroid()->get_cluster()->add_member(checked_curves[i].second); //add curve to its centroid's cluster
            checked_curves[i].second->get_curve_centroid()->get_cluster()->update_cost(checked_curves[i].second->get_dist_to_centroid()); //update cluster's cost
        }

        if ( mean_flag == 0) //PAM has been selected as the update method which means that the centroids are part of the dataset
        {
            if ( arranged_curves.size() == (curves.size() - centroids.size()))
            {
                break;
            }
        }
        else //Mean Discrete Curve has been selected as the update method which means that the centroid are not part of the dataset
        {
            if ( arranged_curves.size() == curves.size())
            {
                break;
            }
        }

        if (checked_curves.size() == 0)
        { //if no curves found to be assigned
            no_assignment_iterations++ ;

            if ( no_assignment_iterations >= 1)
            {
                vector<curve *> remained_curves;
                for(unsigned int i=0; i < curves.size(); i++)
                { //add all remaining curves to a vector
                    if ((curves[i]->is_centroid() == 0)&&(curves[i]->is_assigned()==0))
                    { //must have not been marked as centroid or assigned
                        remained_curves.push_back(curves[i]);
                    }
                }

                simple_clustering(remained_curves, centroids, DFDs); //assign the remaining curves with lloyd's algorithm
                break;
            }
        }
        range_limit *= 2 ; //at each iteration double the range
    }while(1);
}

void assignment_range_search(vector< pair<string, curve *>> &checked_curves,  vector<vector<double> > &DFDs, vector<int>& r_values, multi_hash_table* MHT_ptr,
                          curve * centroid, double &range_limit, unsigned int &k, double &delta, unsigned int &HT_SIZE)
{
    unsigned int HT_key = get_HT_key(r_values, centroid, k, delta, HT_SIZE); //get a hash value for centroid
    vector< pair<string, curve *> > matched_curves;
    MHT_ptr->searching_MHT(matched_curves, HT_key); //searching all hash tables for matches

        double dfd ;
        for(unsigned int i =0 ; i < matched_curves.size(); i++)
        {
            dfd = get_DFD(matched_curves[i].second, centroid, DFDs);
            //get the dfd between centroid and matched curve
            if (dfd < range_limit)
            {
                if (dfd < matched_curves[i].second->get_dist_to_centroid())
                {
                    if (matched_curves[i].second->get_dist_to_centroid() == numeric_limits<double>::infinity())
                    { //add each curve only one time at each range_limit iteration to avoid duplicates
                        checked_curves.push_back(matched_curves[i]);
                    }
                    matched_curves[i].second->set_dist_to_poss_centroid(matched_curves[i].second->get_dist_to_centroid());
                    matched_curves[i].second->set_curve_poss_centroid(matched_curves[i].second->get_curve_centroid());

                    matched_curves[i].second->set_dist_to_centroid(dfd);
                    matched_curves[i].second->set_curve_centroid(centroid);
                }
                else if (dfd < matched_curves[i].second->get_dist_to_poss_centroid())
                {
                    matched_curves[i].second->set_dist_to_poss_centroid(dfd);
                    matched_curves[i].second->set_curve_poss_centroid(centroid);
                }
            }
        }
}


bool find_mean_centroids(vector<curve *> &curves, vector<curve *> &centroids, vector<vector<double> > &DFDs, double bound, double delta, bool lsh_flag)
{
    vector<curve *> new_centroids = centroids;
    double total_cost = get_all_clusters_objective_cost(centroids);
    double min_objective_cost = numeric_limits<double>::infinity();
    unsigned int base_index = centroids[centroids.size() -1]->get_index() + 1; // index in the dfds table

    for(unsigned int i = 0; i < centroids.size(); i++)
    {
        curve * new_centroid;
        vector<curve *> cluster_members = centroids[i]->get_cluster()->get_members(); //get all centroid's cluster members
        vector<curve *> temp_centroids = centroids;

        if (cluster_members.size() == 0)
        {
            continue;
        }
        else if (cluster_members.size() == 1)
        {//if the cluster contains only one curve copy this curve to a new one (possible centroid)
            new_centroid = new curve(cluster_members[0], base_index + i);
        }
        else
        { //if the cluster contains 2 or more curves
            tree * my_tree = new tree(cluster_members); //create a complete binary tree and store all curves at its leafs
            new_centroid = my_tree->get_overall_mean_DF_curve(); //use the tree to get the overall mean curve recursively
            delete my_tree; //the tree is no longer useful so delete it to avoid memory leaks
            new_centroid->set_index(base_index + i); //set the correct index to access the dfds array
        }

        temp_centroids[i] = new_centroid;

        for(unsigned int i=0; i < DFDs.size(); i++) //reset possible centroid's values at the dfds array
        {
            DFDs[new_centroid->get_index()][i] = -1;
            DFDs[i][new_centroid->get_index()] = -1;
        }

        double cost = get_objective_fun_cost(curves, temp_centroids, DFDs); //calculate the objective function cost
        if (cost < min_objective_cost)
        { // if cost is less than the already best
            min_objective_cost = cost;
			for(unsigned int i = 0; i < new_centroids.size(); i++)
			{
				if (new_centroids[i]->get_curve_name() != centroids[i]->get_curve_name())
				{
					delete new_centroids[i];
				}
			}
			
            new_centroids.clear();
            new_centroids = temp_centroids; //keep current configuration of centroids
        }
		else
		{
			delete new_centroid;
		}
    }

    double cost_diff = total_cost - min_objective_cost;
    if (cost_diff > bound)
    { // if new cost is way better than the previous one we need to change some or all centroids
        for(unsigned int i = 0; i < new_centroids.size(); i++)
        {
            if (new_centroids[i]->get_curve_name() != centroids[i]->get_curve_name())
            { // change only the centroids that are different from the previous ones
                set_new_centroid(new_centroids[i], centroids[i]->get_cluster());
                unsigned int old_centroid_index = centroids[i]->get_index();
                unsigned int new_centroid_index = new_centroids[i]->get_index();
                delete centroids[i];

                if ((lsh_flag == 1) && (new_centroids[i]->get_curve_name().at(0) == '!' ))
                {  //if lsh has been selected as assignment method we need to compute the grid curve
                    vector<double> original_curve;
                    transform_curve_points_to_original_curve(new_centroids[i]->get_original_curve_points(), original_curve);
                    new_centroids[i]->set_original_curve(original_curve);

                    vector<double> grid_curve = curve_to_grid_curve_transform(original_curve, new_centroids[i]->get_curve_dimension(), delta); //creating the grid curve
                    new_centroids[i]->set_grid_curve(grid_curve);
                }

                centroids[i] = new_centroids[i];
                centroids[i]->set_index(old_centroid_index);

                for(unsigned int i=0; i < DFDs.size(); i++)
                {
                    DFDs[old_centroid_index][i] = DFDs[new_centroid_index][i] ;
                    DFDs[i][old_centroid_index] = DFDs[i][new_centroid_index];
                }
            }
        }
        return 0; //changes have been made . so send a signal not to finish the program
    }
    else //no changes .keep previous centroid configuration
    {
		for(unsigned int i = 0; i < new_centroids.size(); i++)
        {
            if (new_centroids[i]->get_curve_name() != centroids[i]->get_curve_name())
			{
				delete new_centroids[i];
			}
		}
		
        return 1; //and send a signal to finish the program
    }
}

curve * get_mean_DF_curve(curve * curve1_ptr, curve * curve2_ptr)
{ //based on the pseudocode from eclass notes
    if ( curve2_ptr == NULL)
    {
        return curve1_ptr ;
    }

    vector<vector<double>> &curve1 = curve1_ptr->get_original_curve_points();
    vector<vector<double>> &curve2 = curve2_ptr->get_original_curve_points();

    unsigned int n = curve1.size();
    unsigned int m = curve2.size();

    vector<vector<double> > L(n , vector<double>(m, 0));
    dfd_array(curve1, curve2, L);

    vector<vector<double>> mean_curve_points;
    mean_curve_points.push_back(get_mean_point(curve1[n-1], curve2[m-1]));

    unsigned int i = n - 1 ;
    unsigned int j = m - 1;
    unsigned int index ;
    unsigned int index2;
    double temp_dfd ;

    while (1)
    {
        if  ((i == 0) &&( j == 0))
        {
            break;
        }

        if ( i == 0)
        {
            j--;
        }
        else if ( j == 0)
        {
            i--;
        }
        else
        {
            if( L[i-1][j] < L[i][j-1] )
            {
                index2 = 0;
                temp_dfd = L[i-1][j] ;
            }
            else
            {
                index2 = 1;
                temp_dfd = L[i][j-1] ;
            }

            index = temp_dfd < L[i-1][j-1] ? index2 : 2 ;

            if (index == 0)
            {
                i--;
            }
            else if (index == 1)
            {
                j--;
            }
            else if ( index == 2)
            {
                i--;
                j--;
            }
        }

        mean_curve_points.push_back(get_mean_point(curve1[i], curve2[j]));
    }

    /*
    if (mean_curve_points.size() > 500)
    {
        mean_curve_points.resize(500);
    }*/

    vector<double> vector_double_null;
    curve * mean_curve = new curve(mean_curve_points, mean_curve_points.size(), mean_curve_points[0].size());
    return mean_curve; //new mean curve is not initialized with a corresponding grid curve. it will be set later if needed
}

vector<double> get_mean_point(vector<double> &point1, vector<double> &point2)
{ // get the mean between two curve points of any dimension
    vector<double> mean_point;

    for(unsigned int i = 0 ; i < point1.size(); i++)
    {
        mean_point.push_back( (point1[i] + point2[i]) / 2.0 );
    }

    return mean_point;
}

void transform_curve_points_to_original_curve(vector<vector<double>> &curve_points, vector<double> &original_curve)
{ //transform a vector of points to a vector of coordinates
    for(unsigned int i = 0; i < curve_points.size(); i++)
    {
        for(unsigned int j = 0; j < curve_points[i].size(); j++)
        {
            original_curve.push_back(curve_points[i][j]);
        }
    }
}

bool PAM(vector<curve *> &curves, vector<curve *> &centroids, vector<vector<double> > &DFDs, double bound)
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

            double cost = get_objective_fun_cost(curves, temp_centroids, DFDs);
            cluster_members[j]->set_is_centroid(0);
            unsigned int int_null;
            double double_null;
            double DFD = get_min_DFD(temp_centroids, centroids[i], DFDs, double_null, int_null, int_null);

            cost += DFD;
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
            if (new_centroids[i]->get_curve_name() != centroids[i]->get_curve_name())
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

double get_objective_fun_cost(vector<curve *> &curves, vector<curve *> &temp_centroids, vector<vector<double> > &DFDs)
{
    double total_cost = 0.0 ;
    for(unsigned int i= 0; i < curves.size(); i++)
    {
        if (curves[i]->is_centroid()== 0)
        {
            double sec_min_DFD;
            unsigned int centroid_index ;
            unsigned int poss_centroid_index;

            double min_DFD = get_min_DFD(temp_centroids, curves[i], DFDs, sec_min_DFD, centroid_index, poss_centroid_index);
            total_cost += min_DFD ;
        }
    }

    return total_cost;
}

double get_cluster_cost(vector<curve *> &cluster_members, curve * checking_curve, vector<vector<double> > &DFDs, int checking_member_index)
{ // get the sum of dfds between a curve and a cluster
    double cost = 0.0 ;
    int cluster_members_size = cluster_members.size();
    for(int i=0; i < cluster_members_size; i++)
    {
        if (i == checking_member_index) // if the curve belongs to the cluster to calculate the dfd with itself
        {
            continue;
        }

        cost += get_DFD(cluster_members[i], checking_curve, DFDs);
    }

    return cost;
}


double get_silhouette(vector<curve *> &centroids, vector<vector<double> > &DFDs, vector<double> &clusters_silhts)
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
                get_min_DFD(centroids, cluster_members[j], DFDs, double_null, real_centroid_index, poss_centroid_index);

                curve * real_centroid = centroids[real_centroid_index];
                curve * poss_centroid = centroids[poss_centroid_index];
                curve * cluster_member_centroid = cluster_members[j]->get_curve_centroid();

                if (cluster_member_centroid->get_curve_name() == poss_centroid->get_curve_name())
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
                poss_cluster_cost = get_cluster_cost(poss_cluster_members, cluster_members[j], DFDs, -1);
                b_j = poss_cluster_cost / (poss_cluster_members.size());
            }

            cluster_cost = get_cluster_cost(cluster_members, cluster_members[j], DFDs, j);
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


//============================================================



void parsing_saving_curves(string &file_name, vector<curve *> &curves, double & delta, bool skip_first_line)
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
        if (skip_first_line == 1) //if function being used for parsing queries bypass first line (which defines R , as this is not a team project, so the R is always 0)
        {
            skip_first_line = 0;
            continue;
        }

        vector<vector<double>> original_curve_points ;
        vector<double> point ;

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

        unsigned int coord_counter = 0;
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

            point.push_back(coord);
            coord_counter++;
            if (coord_counter == curve_dimension)
            {
                coord_counter = 0;
                original_curve_points.push_back(point);
                point.clear();
            }
        }
        vector<double> grid_curve = curve_to_grid_curve_transform(original_curve, curve_dimension, delta); //creating the grid curve
        curve * curve_ptr = new curve(original_curve, original_curve_points, grid_curve, curve_name, points, curve_dimension); // save all curve's info to a class
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

int get_uniformly_distributed_random_int(int M, int N)
{  //pseudocode taken from eclass notes
    int result = M + (rand() / (RAND_MAX + 1.0))*(N-M+1);
    return result;
}

double get_uniformly_distributed_random_double2(double M, double N)
{  //pseudocode taken from eclass notes
    double result;
    result = M + (rand() / (RAND_MAX + 1.0))*(N-M);
    return result;
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

double get_DFD(curve * curve1_ptr, curve * curve2_ptr, vector<vector<double> > &DFDs)
{ //checking if the dfd between curve1 and curve2 has already been calculated and stored to the dfds array
            double DFD;
            if ( (DFDs[curve1_ptr->get_index()][curve2_ptr->get_index()] == -1)
                &&(DFDs[curve2_ptr->get_index()][curve1_ptr->get_index()] == -1) )
            {
                DFD = frechet_distance(curve1_ptr->get_original_curve_points(), curve2_ptr->get_original_curve_points());
                DFDs[curve1_ptr->get_index()][curve2_ptr->get_index()] = DFD;
                DFDs[curve2_ptr->get_index()][curve1_ptr->get_index()] = DFD;
            }
            else
            {
                DFD = DFDs[curve1_ptr->get_index()][curve2_ptr->get_index()];
            }

        return DFD;
}

double get_min_DFD(vector<curve * > &curves, curve * query_ptr, vector<vector<double> > &DFDs, double &sec_min_DFD, unsigned int &min_index, unsigned int &sec_min_index)
{   //get the minimum frechet distance between a query curve and a vector of curves
    double min_DFD = numeric_limits<double>::infinity();
    sec_min_DFD = numeric_limits<double>::infinity();

    double DFD;
    for(unsigned int i = 0 ; i < curves.size() ; i++)
    {
        DFD = get_DFD(query_ptr, curves[i], DFDs);

        if (DFD < min_DFD)
        {
            sec_min_DFD = min_DFD ;
            sec_min_index = min_index;

            min_DFD = DFD;
            min_index = i ;
        }
        else
        {
            if (DFD < sec_min_DFD)
            {
                sec_min_DFD = DFD;
                sec_min_index = i ;
            }
        }
/*
        if (min_DFD == 0.0)  // we have found the closest neighbor
        {
            break;
        }*/
     }

     return min_DFD;
}

double frechet_distance(vector<vector<double>> &curve1, vector<vector<double>> &curve2)
{  // implementing the pseudocode
    unsigned int n = curve1.size();
    unsigned int m = curve2.size();

    double L[n][m] ;

    for(unsigned int i=0; i < n ; i++)
    {
        for(unsigned int j=0; j < m ; j++)
        {
            if ((i == 0) && (j == 0))
            {
                L[i][j] = eucl_distance(curve1[0], curve2[0]);
            }
            else if (i == 0)
            {
                double eucl_dist = eucl_distance(curve1[0], curve2[j]);
                L[i][j] = eucl_dist > L[0][j - 1] ? eucl_dist : L[0][j - 1] ;
            }
            else if (j == 0)
            {
                double eucl_dist = eucl_distance(curve1[i], curve2[0]);
                L[i][j] = eucl_dist > L[i-1][0] ? eucl_dist : L[i-1][0] ;
            }
            else
            {
                double min_L1 = L[i][j-1] < L[i-1][j-1] ? L[i][j-1] : L[i-1][j-1] ;
                double min_L2 = min_L1 < L[i-1][j] ? min_L1 : L[i-1][j] ;
                double eucl_dist = eucl_distance(curve1[i], curve2[j]);
                L[i][j] = min_L2 > eucl_dist ? min_L2 : eucl_dist ;
            }
        }
    }

    return L[n-1][m-1] ;
}


void dfd_array(vector<vector<double>> &curve1, vector<vector<double>> &curve2, vector<vector<double> > &L)
{  // getting the array with all dfds between curve1 and curve2
    for(unsigned int i=0; i < L.size() ; i++)
    {
        for(unsigned int j=0; j < L[i].size() ; j++)
        {
            if ((i == 0) && (j == 0))
            {
                L[i][j] = eucl_distance(curve1[0], curve2[0]);
            }
            else if (i == 0)
            {
                double eucl_dist = eucl_distance(curve1[0], curve2[j]);
                L[i][j] = eucl_dist > L[0][j - 1] ? eucl_dist : L[0][j - 1] ;
            }
            else if (j == 0)
            {
                double eucl_dist = eucl_distance(curve1[i], curve2[0]);
                L[i][j] = eucl_dist > L[i-1][0] ? eucl_dist : L[i-1][0] ;
            }
            else
            {
                double min_L1 = L[i][j-1] < L[i-1][j-1] ? L[i][j-1] : L[i-1][j-1] ;
                double min_L2 = min_L1 < L[i-1][j] ? min_L1 : L[i-1][j] ;
                double eucl_dist = eucl_distance(curve1[i], curve2[j]);
                L[i][j] = min_L2 > eucl_dist ? min_L2 : eucl_dist ;
            }
        }
    }
}


double eucl_distance(vector<double> &point1, vector<double> &point2)
{  //calculate the euclidean distance between 2 points of any dimension
    double sum = 0.0;

    for(unsigned int i=0 ; i < point1.size(); i++) //for each points' coordinate (point1 and point2 have the same dimension)
    {
        double d_coor = point1[i] - point2[i] ;
        sum += d_coor * d_coor ;
    }

    double result = sqrt(sum);
    return result;
}

multi_hash_table* initialize_MHT(vector<curve *> &curves, vector<int> &r_values, unsigned int &L, unsigned int &k, unsigned int &HT_SIZE, double &delta)
{
    multi_hash_table* MHT_ptr = new multi_hash_table(L, HT_SIZE);

    for (unsigned int i=0; i<curves.size() ; i++) //for each dataset's curve
    {
        LSH(r_values, MHT_ptr, curves[i], L, k, HT_SIZE, delta);  // hash curve and store it to L hash tables
    }

    return MHT_ptr;
}



void print_all_curves(vector<curve *> &curves, bool print_coords)
{
    for (unsigned int i=0; i < curves.size() ; i++)
    {
        curves[i]->print_curve(print_coords);
    }
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


void print_clusters(vector<curve *>& centroids)
{
    cout<<"-----PRINTING CLUSTERS-----"<<endl;
    for(unsigned int i =0 ; i < centroids.size(); i++)
    {
        centroids[i]->get_cluster()->print(centroids[i]->get_curve_name());
        cout<<"---------------------------------------"<<endl;
    }
}


void initialize_random_seed() //for rand function
{
    srand (time(NULL));
}
