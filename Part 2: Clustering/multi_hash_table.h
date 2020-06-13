#ifndef MULTI_HASH_TABLE_H
#define MULTI_HASH_TABLE_H

#include <vector>
#include <iostream>
#include <algorithm>
#include "curve.h"

using namespace std;

struct bucket  //hash table's bucket
{
    vector<curve * > bucket_curve;  // basically a bucket is a vector that contains pointers to already saved curves
    void print_bucket();
    void add_curve_to_bucket(curve * curve_ptr);
    void get_all_curves(vector< pair<string, curve *> >  &matched_curves);
};

struct hash_table
{
    unsigned int HT_SIZE;
    vector<bucket *> bucket_index; //basically a hash table is a vector that contains pointers to its buckets
    hash_table(unsigned int & HT_SIZE);
    ~hash_table();
    void print_HT();
    void store_curve_to_HT(curve * curve_ptr, unsigned int &HT_key);
    void get_all_curves_from_bucket(vector< pair<string, curve *> >  &matched_curves, unsigned int& HT_key);
};

struct multi_hash_table
{
    unsigned int L;
    vector <hash_table *> hash_table_id ; // basically a multi hash table is a vector that contains pointers to its hash tables
    multi_hash_table(unsigned int &L, unsigned &HT_SIZE);
    ~multi_hash_table();
    void print_MHT();
    void searching_MHT(vector< pair<string, curve *> >  &matched_curves, unsigned int& HT_key);
};

#endif // MULTI_HASH_TABLE_H
