#include "bst.h"
#include "lsh.h"

using namespace std;

void tree::post_order_print(node * leaf)
{
    if (leaf != NULL)
    {
        post_order_print(leaf->left);
        post_order_print(leaf->right);
        if ( leaf->curve_ptr != NULL)
        {
            cout<<leaf->curve_ptr->get_curve_name()<<endl;
        }
    }

}

void tree::print()
{
    post_order_print(root);
}


tree::tree(vector<curve *> curves) //create binary tree - all entries will be saved at the leafs
{
    queue<pair<unsigned int, node *>> nodes_queue; //use a queue as a helping structure

    do // first create a tree node for each curve and add them to the queue
    { //the queue stores pairs that include the number of sub-nodes and the tree node (as a pointer)
        node * leaf = new node ;
        leaf->left = NULL;
        leaf->right = NULL;
        leaf->curve_ptr = curves[curves.size() -1];
        curves.pop_back();
        nodes_queue.emplace(1, leaf);

        if(curves.size() == 0)
        {
            break;
        }
    }while(1);

    do //start poping nodes from the queue, creating their parent node and making them its children
    { //using pointers
        pair<unsigned int, node *> first_node = nodes_queue.front();
        nodes_queue.pop();
        if (nodes_queue.size()==0) //stop when queue is empty
        {
            root = first_node.second; //tree has been created
            break;
        }

        if (first_node.first < nodes_queue.front().first) //it must be a complete binary tree so keep in mind the number of each node's sub-nodes
        {
            nodes_queue.push(first_node);
        }
        else
        {
            pair<unsigned int, node *> second_node = nodes_queue.front();
            nodes_queue.pop();

            node * iner_node = new node;
            iner_node->left = first_node.second ;
            iner_node->right = second_node.second;
            iner_node->curve_ptr = NULL;

            nodes_queue.emplace(first_node.first + second_node.first, iner_node);
        }
    }while(1);
}

tree::~tree() //starting from the root recursively
{
    destroy(root);
}

void tree::destroy(node * node_ptr)
{
    if (node_ptr!=NULL)
    {
        destroy(node_ptr->left);
        destroy(node_ptr->right);
        delete node_ptr;
    }
}


curve* tree::post_order_traversal(node * node_ptr, vector<curve *> &temp_centroids)
{ //based on the pseudocode from eclass notes
    if (node_ptr->curve_ptr != NULL) //if tree node is a leaf
    {
        return node_ptr->curve_ptr;
    }
    else
    {
        curve * left_curve = post_order_traversal(node_ptr->left, temp_centroids);

        curve * right_curve;
        if (node_ptr->right != NULL)
        {
            right_curve = post_order_traversal(node_ptr->right, temp_centroids);
        }
        else
        {
            right_curve = NULL;
        }

        curve * temp_centroid = get_mean_DF_curve(left_curve, right_curve);
        node_ptr->curve_ptr = temp_centroid;
        if (right_curve != NULL)
        {
            temp_centroids.push_back(temp_centroid); //keep temporally created curves into a vector so that they can been deleted later avoiding memory leaks
        }
        return temp_centroid;
    }
}

curve * tree::get_overall_mean_DF_curve()
{
    vector<curve *> temp_centroids;
    curve * overall_centroid = post_order_traversal(root, temp_centroids);
    for(unsigned int i = 0; i < temp_centroids.size() - 1; i++)
    { //delete all temporally created curves except for the last one because it is the overall mean curve (the new possible centroid of the cluster)
        delete temp_centroids[i];
        temp_centroids[i] = NULL;
    }
    temp_centroids.clear();

    return overall_centroid;
}

