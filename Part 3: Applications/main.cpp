#include "header.h"
#include "cluster.h"

int main(int argc, char *argv[])
{

    InputParser input(argc, argv);

    string input_file;
    string input_file2;
    string dist_name;

    if(input.cmdOptionExists("-i1"))
    {
        const string &arg_str = input.getCmdOption("-i1");
        input_file = arg_str;
    }
    else
    {
        cout<<"->Please enter input file name for proteins:"<<endl;
        cin >> input_file;
    }

    if(input.cmdOptionExists("-i2"))
    {
          const string &arg_str = input.getCmdOption("-i2");
          input_file2 = arg_str;
    }
    else
    {
        cout<<"->Please enter input file name for roads:"<<endl;
        cin >> input_file2;
    }

    if(input.cmdOptionExists("-d"))
    {
        const string &arg_str = input.getCmdOption("-d");
        string dist_name = arg_str;
        if ((dist_name != "DFD") && (dist_name != "CRMSD"))
        {
            if(dist_name == "DTW")
            {
                cout<<"->You chose DTW as the distance metric but this is not a team project, so it is not supported."<<endl;
            }
            else
            {
                cout<<"-WRONG- input as a distance metric!"<<endl;
            }
            cout<<"->Distance metric set by default"<<endl;
            dist_name = "CRMSD";
        }
    }
    else
    {
        dist_name = "CRMSD";
    }

    initialize_random_seed();

    vector<curve *> proteins;
    vector<curve *> roads;
    vector<curve *> segs;


    for(unsigned int i = 0; i < 2; i++) // i = 0 for first part and i = 1 for second part
    {
        string output_fn ;
        bool dist_choice;

        cout<<"-----------------------------------------"<<endl;
        cout<<"->i: "<<i<<endl;

        if (i == 0)
        {
            cout<<"--->Processing and clustering proteins<---"<<endl;
            cout<<"->Input file name: "<<input_file<<endl;

            parse_proteins(input_file, proteins, 3);

            if ( dist_name == "CRMSD")
            {
                cout<<"->Distance metric: c-RMSD"<<endl;
                output_fn = "crmsd.dat" ;
                dist_choice = 0;
            }
            else
            {
                cout<<"->Distance metric: DFD"<<endl;
                output_fn = "frechet.dat" ;
                dist_choice = 1;
            }
        }
        else //i=1
        {
            delete_curves(proteins);

            dist_choice = 1;

            cout<<"--->Processing and clustering roads' segments<---"<<endl;
            cout<<"->Distance metric: DFD"<<endl;
            cout<<"->Input file name: "<<input_file2<<endl;

            output_fn = "kmeans_ways_frechet.dat" ;

            parse_roads(input_file2, roads, 0, 2);
            seperate_roads_to_segments(roads, segs, 0.03, 10, 20);
            delete_curves(roads);

            cout<<"->Roads have been divided in segments"<<endl;

            create_segments_file(segs, "segment.csv");
        }

        vector<vector<double>> dists;
        unsigned int k_ub;

        if (i == 0)
        {
            vector<vector<double>> temp_dists( proteins.size(), vector<double>(proteins.size(), -1));
            dists = temp_dists;
            k_ub = floor(proteins.size() * 0.10) ;
        }
        else
        {
            vector<vector<double>> temp_dists( segs.size(), vector<double>(segs.size(), -1));
            dists = temp_dists;
            k_ub = floor(segs.size() * 0.10) ;
        }

        unsigned int k_lb = 4;
        unsigned int k = k_lb;
        double k_lb_silht;
        double k_ub_silht;
        double mid_k;

        vector<pair<double, unsigned int>> silht_k ;
        unsigned int j = 0;

        clock_t begin = clock();
        do
        {
            vector<curve *> centroids;

            if (i == 0)
            {
                centroids = k_means_pp(proteins, dists, k, dist_choice);
            }
            else
            {
                centroids = k_means_pp(segs, dists, k, dist_choice);
            }

            cout<<"--->Chose "<< k << " centroids<---"<<endl;
            create_clusters(centroids);

            double bound = 0.002;
            unsigned int iterations = 0;
            do
            {
                if (iterations == 4)
                {
                    iterations = 0;
                    bound *= 4;
                }

                bool end_program;
                if (i == 0)
                {
                    simple_clustering(proteins, centroids, dists, dist_choice);
                    end_program = PAM(proteins, centroids, dists, bound, dist_choice);
                }
                else
                {
                    simple_clustering(segs, centroids, dists, dist_choice);
                    end_program = PAM(segs, centroids, dists, bound, dist_choice);
                }

                //cout<<"--->Done updating<---"<<endl;

                if ( end_program == 1)
                {
                    break;
                }
                else
                {
                    if ( i == 0)
                    {
                        reset_all_curves(proteins);
                    }
                    else
                    {
                        reset_all_curves(segs);
                    }

                    reset_all_clusters(centroids);
                }

                iterations++;
            }while(1);

            vector<double> clusters_silhts;
            double overall_silht = get_silhouette(centroids, dists, clusters_silhts, dist_choice);

            if (j == 0)
            {
                silht_k.emplace_back(overall_silht, k);
                k = k_ub;
                k_lb_silht = overall_silht;
            }
            else if (j == 1)
            {
                silht_k.emplace_back(overall_silht, k);
                k_ub_silht = overall_silht ;
                k = k_lb + (k_ub - k_lb) / 2 + (k_ub - k_lb) % 2;
            }
            else
            {
                bool exit = 0;

                if ( (k == k_lb) || (k == k_ub) )
                {
                    exit = 1;
                }
                else if( (k_ub - k_lb <= 4) && (abs(k_ub_silht - k_lb_silht) <= 0.025) )
                {
                    exit = 1;
                }

                if (exit == 1)
                {
                    clock_t end = clock();
                    double elapsed_time = double(end - begin) / CLOCKS_PER_SEC;

                    write_results_in_file(centroids, overall_silht, output_fn, elapsed_time, i);

                    reset_all_clusters(centroids);
                    destroy_clusters(centroids, 0);

                    if(i == 0)
                    {
                        reset_curves(proteins);
                    }
                    else
                    {
                        reset_curves(segs);
                    }

                    break;
                }


                silht_k.emplace_back(overall_silht, k);
                sort(silht_k.begin(), silht_k.end()) ;

                if (silht_k[1].second < silht_k[2].second)
                {
                    k_lb = silht_k[1].second;
                    k_lb_silht = silht_k[1].first;

                    k_ub = silht_k[2].second;
                    k_ub_silht = silht_k[2].first;
                }
                else
                {
                    k_lb = silht_k[2].second;
                    k_lb_silht = silht_k[2].first;

                    k_ub = silht_k[1].second;
                    k_ub_silht = silht_k[1].first;
                }

                mid_k = k_lb + (k_ub - k_lb) / 2 + (k_ub - k_lb) % 2;


                if (mid_k == k)
                {
                    if(k_lb_silht < k_ub_silht)
                    {
                        k_lb = mid_k ;
                        k_lb_silht = overall_silht;
                    }
                    else
                    {
                        k_ub = mid_k ;
                        k_ub_silht = overall_silht;
                    }

                    mid_k = k_lb + (k_ub - k_lb) / 2 + (k_ub - k_lb) % 2;
                }

                k = mid_k ;
                silht_k.clear();
                silht_k.emplace_back(k_lb_silht, k_lb);
                silht_k.emplace_back(k_ub_silht, k_ub);

            }

            reset_all_clusters(centroids);
            destroy_clusters(centroids, 0);

            if(i == 0)
            {
                reset_curves(proteins);
            }
            else
            {
                reset_curves(segs);
            }

            j++;
        }while(1);
    }

    delete_curves(segs);
    return 0;
}


void write_results_in_file(vector<curve *> &centroids, double overall_silht, string output_fn, double clustering_time, bool segments_flag)
{
    ofstream outfile (output_fn, ios_base::app | ios_base::out);
    outfile.precision(8);

    outfile << "k: " << centroids.size()<<endl;
    outfile << "s: " << overall_silht <<endl;

    if(segments_flag == 1)
    {
        outfile << "clustering_time: " << clustering_time <<endl;
    }

    outfile << "------------------------------------------------"<<endl;


        for(unsigned int i=0; i < centroids.size(); i++)
        {
            vector<curve *> cluster_members = centroids[i]->get_cluster()->get_members();
            vector<unsigned int > cluster_members_ids;

            for(unsigned int j=0; j < cluster_members.size(); j++)
            {
                cluster_members_ids.push_back(cluster_members[j]->get_id());
            }

            sort(cluster_members_ids.begin(), cluster_members_ids.end()) ;

            for(unsigned int j=0; j < cluster_members_ids.size(); j++)
            {
                outfile<<cluster_members_ids[j];
                outfile<<"    ";
            }

            outfile<<endl;
            outfile << "------------------------------------------------"<<endl;

        }

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


void print_all_curves(vector<curve *> &curves, bool print_coords)
{
    for (unsigned int i=0; i < curves.size() ; i++)
    {
        curves[i]->print(print_coords);
    }
}


void print_points(vector<vector<double>> &original_curve_points,unsigned int dimension)
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

vector<vector<double>> cut_off_points(vector<vector<double>> points, unsigned int n_points)
{
    points.resize(n_points);
    return points;
}

vector<vector<double>> select_random_points(vector<vector<double>> points, unsigned int n_points)
{
    vector<vector<double>> random_points;

    for(unsigned int i = 0 ; i < n_points ; i++)
    {
        unsigned int point_i = get_uniformly_distributed_random_int(0, points.size() -1 ) ;
        random_points.push_back(points[point_i]);
        points.erase(points.begin() + point_i) ;
    }

    return random_points ;
}


void delete_curves(vector<curve * > &curves)
{
    for (unsigned int i=0; i < curves.size(); i++)
    {
        delete curves[i];
    }
}



void print_clusters(vector<curve *>& centroids)
{
    cout<<"-----PRINTING CLUSTERS-----"<<endl;
    for(unsigned int i =0 ; i < centroids.size(); i++)
    {
        //centroids[i]->get_cluster()->print(centroids[i]->get_curve_name());
        cout<<"---------------------------------------"<<endl;
    }
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

void initialize_random_seed() //for rand function
{
    srand (time(NULL));
}
