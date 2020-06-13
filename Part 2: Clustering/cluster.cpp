#include "cluster.h"

cluster::cluster(string id)
{
    this->id = id;
    this->total_cost = 0.0 ;
}

cluster::~cluster()
{
    //dtor
}


void cluster::print(string centroid_name)
{
    cout.precision(15);

    cout<<"Cluster's id:"<<id<<endl;
    cout<<"Cluster's centroid:"<<centroid_name<<endl;
    cout<<"Cluster's total cost:"<<fixed<<total_cost<<endl;
    cout<<"Cluster's members:"<<endl;
    cout<<"[ ";
    for(unsigned int i=0; i < members.size(); i++)
    {
        cout<<members[i]->get_curve_name();
        if (i != members.size() - 1)
        {
            cout<<", ";
        }
    }
    cout<<" ]"<<endl;
}

string cluster::get_id()
{
    return id;
}

unsigned int cluster::get_cluster_size()
{
    return members.size();
}

void cluster::update_cost(double cost)
{
    total_cost += cost ;
}

double cluster::get_total_cost()
{
    return total_cost;
}

vector<curve *> cluster::get_members()
{
    return members;
}

void cluster::set_total_cost(double cost)
{
    total_cost = cost;
}

void cluster::add_member(curve * curve_ptr)
{
    members.push_back(curve_ptr);
}

void cluster::delete_members()
{
    members.clear();
}

