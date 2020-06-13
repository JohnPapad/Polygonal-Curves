
#include "header.h"


void parse_proteins(string input_fn, vector<curve *> &proteins, unsigned int d)
{
    ifstream file(input_fn);
    if (!file) //checking if file exists
    {
        cout<<"-ERROR 404-No input file with name: '"<<input_fn<<"' found"<<endl;
        exit(1);
    }

    unsigned int num_conform;
    unsigned int number_of_points;

    string line;

    getline(file, line);
    istringstream iss(line);
    iss>>num_conform;

    getline(file, line);
    istringstream iss2(line);
    iss2>>number_of_points;

    for(unsigned int i=0; i < num_conform ; i++)
    {
        vector<vector<double>> protein_points;

        vector<double> c_point(d, 0.0);

        for(unsigned int j = 0; j < number_of_points; j++)
        {
            vector<double> point;
            getline(file, line);
            istringstream iss(line);

            for(unsigned int k=0 ; k < d; k++)
            {
                double coor;
                iss>>coor;
                point.push_back(coor);

                c_point[k] += coor;
            }

            protein_points.push_back(point);
        }


        for(unsigned int i=0 ; i < d; i++)
        {

            c_point[i] = c_point[i] / number_of_points;
        }

        curve * protein_ptr = new curve(protein_points, c_point, "protein", i, number_of_points, d);
        proteins.push_back(protein_ptr);

    }
}



