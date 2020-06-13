#include "header.h"
#include "distances.h"

double calc_dist(curve * x_ptr, curve * y_ptr, bool dist_choice)
{

	vector<vector<double>> &x = x_ptr->get_original_curve_points();
	vector<double> &x_c = x_ptr->get_c_point();

	vector<vector<double>> &y = y_ptr->get_original_curve_points();
	vector<double> &y_c = y_ptr->get_c_point();

    if (dist_choice == 1)
    {
        return frechet_distance(x, y);
    }

	unsigned int dimension = x_ptr->get_curve_dimension();
	unsigned int x_number_of_points = x_ptr->get_number_of_points();
	unsigned int y_number_of_points = y_ptr->get_number_of_points();

	unsigned int number_of_points = x_number_of_points;

    if(x_number_of_points != y_number_of_points)
    {
        if ( x_number_of_points < y_number_of_points)
        {
            number_of_points = x_number_of_points;
            y = cut_off_points(y, number_of_points);

        }
        else
        {
            number_of_points = y_number_of_points;
            x = cut_off_points(x, number_of_points);
        }

        y_c = calc_point_c(y);
        x_c = calc_point_c(x);
    }

    MatrixXd X(number_of_points , dimension), Y(number_of_points , dimension);

    for(unsigned int i=0; i < number_of_points; i++)
    {
    	for(unsigned int j=0 ; j < dimension ; j++)
    	{
    		X(i, j) = x[i][j] - x_c[j];
    		Y(i, j) = y[i][j] - y_c[j];
    	}
    }

    JacobiSVD<MatrixXd> svd( X.transpose() * Y, ComputeFullU | ComputeFullV);

	JacobiSVD<MatrixXd>::SingularValuesType singular = svd.singularValues();

	if ((singular(dimension-1) <= 0) && (dist_choice == 0))
	{
        cout<<"Singular value 3 <= 0"<<endl;
		exit(-3);
	}

	MatrixXd Q(dimension , dimension);
	Q = svd.matrixU() * svd.matrixV().transpose();

	if ( Q.determinant() < 0)
	{
		MatrixXd U(dimension, dimension);
		U = svd.matrixU();

		for(unsigned int i = 0 ; i < dimension ;i++)
		{
			U(i, dimension - 1) *= -1;
		}

		Q = U * svd.matrixV().transpose();
	}

    MatrixXd XQ(number_of_points, dimension) ;
	XQ = X * Q ;

    if ( dist_choice == 0)
    {
        XQ = XQ - Y ;
        return XQ.norm() / sqrt(number_of_points) ;
	}
	else
	{
        vector<vector<double>> Y_vec = convert_matrix_to_vector(Y);
        vector<vector<double>> XQ_vec = convert_matrix_to_vector(XQ);

        return frechet_distance(XQ_vec, Y_vec);
	}

}


vector<double> calc_point_c(vector<vector<double>> &points)
{
    vector<double> c_point(points[0].size(), 0.0);

    for(unsigned int i=0; i < points.size(); i++)
    {
        for(unsigned int j=0; j < points[i].size(); j++)
        {
            c_point[j] += points[i][j];
        }
    }

    for(unsigned int i=0 ; i < points[0].size(); i++)
    {
        c_point[i] = c_point[i] / points.size();
    }

    return c_point;

}

vector<vector<double>> convert_matrix_to_vector(MatrixXd &mtrx)
{
    vector<vector<double>> points( mtrx.rows(), vector<double>(mtrx.cols(), 0.0));

    for(unsigned int i=0 ; i < mtrx.rows(); i++)
    {
        for(unsigned int j=0; j < mtrx.cols(); j++)
        {
            points[i][j] = mtrx(i, j);
        }
    }

    return points;
 }


double get_dist(curve * curve1_ptr, curve * curve2_ptr, bool dist_choice, vector<vector<double> > &dists)
{ //checking if the dfd between curve1 and curve2 has already been calculated and stored to the dfds array
            double dist;
            if ( (dists[curve1_ptr->get_index()][curve2_ptr->get_index()] == -1)
                &&(dists[curve2_ptr->get_index()][curve1_ptr->get_index()] == -1) )
            {
                dist = calc_dist(curve1_ptr, curve2_ptr, dist_choice);

                dists[curve1_ptr->get_index()][curve2_ptr->get_index()] = dist;
                dists[curve2_ptr->get_index()][curve1_ptr->get_index()] = dist;
            }
            else
            {
                dist = dists[curve1_ptr->get_index()][curve2_ptr->get_index()];
            }

        return dist;
}




double get_min_dist(vector<curve * > &curves, curve * query_ptr, vector<vector<double> > &dists, double &sec_min_dist, unsigned int &min_index, unsigned int &sec_min_index, bool dist_choice)
{   //get the minimum frechet distance between a query curve and a vector of curves
    double min_dist = numeric_limits<double>::infinity();
    sec_min_dist = numeric_limits<double>::infinity();

    double dist;
    for(unsigned int i = 0 ; i < curves.size() ; i++)
    {
        dist = get_dist(query_ptr, curves[i], dist_choice, dists);

        if (dist < min_dist)
        {
            sec_min_dist = min_dist ;
            sec_min_index = min_index;

            min_dist = dist;
            min_index = i ;
        }
        else
        {
            if (dist < sec_min_dist)
            {
                sec_min_dist = dist;
                sec_min_index = i ;
            }
        }
/*
        if (min_DFD == 0.0)  // we have found the closest neighbor
        {
            break;
        }*/
     }

     return min_dist;
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
