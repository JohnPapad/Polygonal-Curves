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


void bucket::get_all_curves(vector< pair<string, curve *> > &matched_curves) // add all bucket's curve pointers to the matched vector
{
    for(unsigned int i=0; i < bucket_curve.size() ; i++)
    {
        if ((bucket_curve[i]->is_centroid() == 0)&&(bucket_curve[i]->is_assigned()==0))
        {
            if (!(find(matched_curves.begin(), matched_curves.end(), make_pair(bucket_curve[i]->get_curve_name(), bucket_curve[i])) != matched_curves.end()))
            { //checking for duplicates among the buckets . Also curves must not be a centroid. Adding only the ones what have not been assigned yet
                matched_curves.emplace_back(bucket_curve[i]->get_curve_name(), bucket_curve[i]) ;
            }
        }
    }
}


void hash_table::get_all_curves_from_bucket(vector< pair<string, curve *> >  &matched_curves, unsigned int& HT_key)
{
    bucket_index[HT_key]->get_all_curves(matched_curves);
}


void multi_hash_table::searching_MHT(vector< pair<string, curve *> >  &matched_curves, unsigned int &HT_key)
{   //searching all L hash tables' buckets that match with the HT_key.All matches are saved to the matched vector.No duplicates allowed
    for(unsigned int i=0; i < L ; i++)
    {
        hash_table_id[i]->get_all_curves_from_bucket(matched_curves, HT_key);
    }
}
