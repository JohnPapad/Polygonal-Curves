#include "header.h"


void parse_roads(string &file_name, vector<curve *> &roads, bool skip_first_line, unsigned int d)
{
    ifstream file(file_name);

    if (!file) //checking if file exists
    {
        cout<<"-ERROR-There is no input file: '"<<file_name<<"'"<<endl;
        exit(1);
    }

    string line;

    while(getline(file, line))
    {
        if (skip_first_line == 1)
        {
            skip_first_line = 0;
            continue;
        }

        vector<vector<double>> original_curve_points ;
        vector<double> point ;

        vector<double> c_point(d, 0.0);


        istringstream iss(line);
        string type;
        unsigned int id;

        string id_str;
        getline(iss, id_str, ' ');
        id_str.pop_back();
        id = stoi(id_str);

        string type_str;
        getline(iss, type_str, ' ');
        type_str.pop_back();
        type = type_str;

        string coord_str;

        unsigned int coord_counter = 0;
        while(getline(iss, coord_str, ' '))
        {
            double coord;
            if (coord_str.back() == ',')
            {
                coord_str.pop_back();
            }
            coord = stod(coord_str);

            point.push_back(coord);
            coord_counter++;
            if (coord_counter == d)
            {
                coord_counter = 0;
                original_curve_points.push_back(point);
                point.clear();
            }
        }
        curve * road_ptr = new curve(original_curve_points, c_point, type, id, original_curve_points.size(), d); // save all curve's info to a class
        roads.push_back(road_ptr);
    }
}


void seperate_roads_to_segments(vector<curve *> &roads, vector<curve *> &segs, double radius_bound, unsigned int min_seg_size, unsigned int max_seg_size)
{
    for(unsigned int i=0; i < roads.size(); i++)
    {
        vector<vector<double>> &road_nodes = roads[i]->get_original_curve_points();
        vector<vector<double>> seg_nodes = road_nodes;

        if (seg_nodes.size() <= min_seg_size)
        {
            curve * seg_ptr = new curve(seg_nodes, roads[i]->get_id(), seg_nodes.size(), roads[i]->get_curve_dimension());
            segs.push_back(seg_ptr);
            continue;
        }
        else
        {
            seg_nodes.clear();
        }

        for(unsigned int j = 0; j < road_nodes.size(); j ++ )
        {
            bool cut_off = 0;
            if ( (seg_nodes.size() <= (min_seg_size / 2)) || ( j >= (road_nodes.size() - min_seg_size / 2) ) )
            {
                seg_nodes.push_back(road_nodes[j]);
                continue;
            }
            else
            {
                if( check_for_junction(roads, road_nodes[j], roads[i]->get_type(), i) == 1 )
                {
                    cut_off = 1;
                }
                else if( get_circumcircle_radius(road_nodes[j-1], road_nodes[j], road_nodes[j+1]) > radius_bound )
                {
                    cut_off = 1;
                }
                else
                {
                    seg_nodes.push_back(road_nodes[j]);
                }
            }

            if ( (seg_nodes.size() == max_seg_size) || (cut_off == 1) )
            {
                curve * seg_ptr = new curve(seg_nodes, roads[i]->get_id(), seg_nodes.size(), roads[i]->get_curve_dimension());
                segs.push_back(seg_ptr);
                seg_nodes.clear();
            }
        }

        if (seg_nodes.size() != 0)
        {
            curve * seg_ptr = new curve(seg_nodes, roads[i]->get_id(), seg_nodes.size(), roads[i]->get_curve_dimension());
            segs.push_back(seg_ptr);
            seg_nodes.clear();
        }

    }
}


double get_circumcircle_radius(vector<double> &point_a, vector<double> &point_b, vector<double> &point_c)
{
    double a;
    double b;
    double c;

    a = eucl_distance(point_b, point_c);
    b = eucl_distance(point_a, point_c);
    c = eucl_distance(point_a, point_b);

    double sqrt_res = (a+b+c) * (b + c -a) * (c + a -b) * (a+b-c) ;
    if (sqrt_res == 0)
    {
        return 0;
    }

    double res = (a * b * c) / sqrt( sqrt_res );
    return res;

}

bool check_for_junction(vector<curve *> &roads, vector<double> &checking_node, string checking_node_type, unsigned int checking_node_i)
{
    for(unsigned int i = 0; i < roads.size(); i++)
    {
        if (i  ==  checking_node_i)
        {
            continue;
        }

        if (roads[i]->get_type() == checking_node_type)
        {
            vector<vector<double>> &road_i_nodes = roads[i]->get_original_curve_points();

            if ( check_range(road_i_nodes[0], road_i_nodes[road_i_nodes.size()-1], checking_node) == 1 )
            {
                for(unsigned int j=0; j < road_i_nodes.size(); j++)
                {
                    if (compare_nodes(road_i_nodes[j], checking_node) == 1)
                    {
                        return 1;
                    }
                }
            }
        }

    }

    return 0;
}


bool check_range(vector<double> &first_node, vector<double> &last_node, vector<double> &checking_node)
{
    for(unsigned int i=0; i < first_node.size(); i++)
    {
        double min_coor;
        double max_coor;
        if (first_node[i] < last_node[i])
        {
            min_coor = first_node[i];
            max_coor = last_node[i];
        }
        else
        {
            max_coor = first_node[i];
            min_coor = last_node[i];
        }

        if ( ! ((min_coor <= checking_node[i]) && (checking_node[i] <= max_coor)) )
        {
            return 0;
        }
    }

    return 1;
}

bool compare_nodes(vector<double> &node1, vector<double> &node2)
{
    for(unsigned int i=0; i < node1.size(); i++)
    {
        if (node1[i] != node2[i])
        {
            return 0;
        }
    }

    return 1;
}


void create_segments_file(vector<curve *> &segs, string output_fn)
{
    ofstream outfile (output_fn, ios_base::app | ios_base::out);
    outfile.precision(8);

        for(unsigned int i=0; i < segs.size(); i++)
        {
            outfile<<segs[i]->get_id()<<", "<<segs[i]->get_type()<<", "<<segs[i]->get_number_of_points()<<", ";

            vector<vector<double>> &seg_nodes = segs[i]->get_original_curve_points();
            for(unsigned int j=0; j < seg_nodes.size(); j++)
            {
                for(unsigned int k=0; k < seg_nodes[j].size(); k++)
                {
                    outfile<<seg_nodes[j][k]<<", " ;
                }
            }

            outfile<<endl;
        }

}

