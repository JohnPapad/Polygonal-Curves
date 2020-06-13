#include "multi_hash_table.h"

multi_hash_table::multi_hash_table(unsigned int &L, unsigned int &HT_SIZE)
{
    this->L = L;
    for (unsigned int i=0; i < L ; i++) //dynamically create all L hash tables and store their pointers to a vector
    {
        hash_table * HT_ptr = new hash_table(HT_SIZE);
        hash_table_id.push_back(HT_ptr);
    }
}


multi_hash_table::~multi_hash_table()
{  //taking care of memory leaks
    for (unsigned int i=0; i < L ; i++)
    {
        delete hash_table_id[i];
    }
}


hash_table::hash_table(unsigned int &HT_SIZE)
{
    this->HT_SIZE = HT_SIZE;
    for (unsigned int i=0; i < HT_SIZE ; i++) //dynamically create all HT_SIZE buckets and store their pointers to a vector
    {
        bucket * bucket_ptr = new bucket;
        bucket_index.push_back(bucket_ptr);
    }
}


hash_table::~hash_table()
{   //taking care of memory leaks
    for (unsigned int i=0; i < HT_SIZE ; i++)
    {
        delete bucket_index[i];
    }
}


void multi_hash_table::print_MHT()
{
    cout<<"-----PRINTING MULTI HASH TABLE-----"<<endl;
    for(unsigned int i=0 ; i < L ; i++)
    {
        cout<<"---Printing hash table #"<<i<<"---"<<endl;
        hash_table_id[i]->print_HT();
    }
    cout<<"------------------------------------"<<endl;
}


void hash_table::print_HT()
{
    for(unsigned int i=0 ; i < HT_SIZE ; i++)
    {
        cout<<"->Printing bucket #"<<i<<endl;
        bucket_index[i]->print_bucket();
    }
    cout<<"------------------------------------"<<endl;
}


void bucket::print_bucket()
{
    cout<<"[";
    for(unsigned int i=0 ; i < bucket_curve.size() ; i++)
    {
        cout<<bucket_curve[i]->get_curve_name();
        if (bucket_curve.size()-1 != i)
            cout<<",";
    }
    cout<<"]"<<endl;
}


void bucket::add_curve_to_bucket(curve * curve_ptr) //add a curve's pointer to a bucket
{
    bucket_curve.push_back(curve_ptr);
}


void hash_table::store_curve_to_HT(curve * curve_ptr, unsigned int &HT_key)
{
    bucket_index[HT_key]->add_curve_to_bucket(curve_ptr);
}


void bucket::searching_bucket(curve * query_ptr, vector<curve *> &matched_curves)
{   //searching the entire bucket in order to find curves that match with the query curve (the comparison is conducted by comparing grid curves)
    for(unsigned int i=0; i < bucket_curve.size() ; i++)
    {
        if (bucket_curve[i]->compare_grid_curves(query_ptr) == 1) // if there is a match
        {
            matched_curves.push_back(bucket_curve[i]) ;  // add it to the according vector
        }
    }
}


void bucket::get_all_curves(vector<curve *> &matched_curves) // add all bucket's curve pointers to the matched vector
{
    for(unsigned int i=0; i < bucket_curve.size() ; i++)
    {
        matched_curves.push_back(bucket_curve[i]) ;
    }
}


void hash_table::searching_HT(curve * query_ptr, vector<curve *> &matched_curves, unsigned int& HT_key)
{
    bucket_index[HT_key]->searching_bucket(query_ptr, matched_curves);
}


void hash_table::get_all_curves_from_bucket(vector<curve *> &matched_curves, unsigned int& HT_key)
{
    bucket_index[HT_key]->get_all_curves(matched_curves);
}


void multi_hash_table::searching_MHT(curve * query_ptr, vector<curve *> &matched_curves, unsigned int &HT_key, int &found_grid_curve)
{   //searching all L hash tables in order to find curves that match with the query curve.All matches are saved to the matched vector
    for(unsigned int i=0; i < L ; i++)
    {
        hash_table_id[i]->searching_HT(query_ptr, matched_curves, HT_key);
    }

    if (matched_curves.size() == 0) // if no match found get all curves from all buckets, with index = HT_KEY, from all hash tables
    {
        found_grid_curve = 0;

        for(unsigned int i=0; i < L ; i++)
        {
            hash_table_id[i]->get_all_curves_from_bucket(matched_curves, HT_key);
        }
    }
    else
    {
        found_grid_curve = 1;
    }
}
