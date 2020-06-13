#ifndef BST_H
#define BST_H

#include <iostream>
#include "curve.h"
#include <queue>

using namespace std;

class tree
{
    private:
        struct node
        {
           curve * curve_ptr ;
           node * left;
           node * right;
        };

        node * root ;

        void destroy(node * node_ptr);
        void post_order_print(node * leaf);
        curve* post_order_traversal(node * node_ptr, vector<curve *> &temp_centroids);


    public:
         tree(vector<curve *> curves);
        ~tree();
         void print();
         curve * get_overall_mean_DF_curve();


};

#endif
