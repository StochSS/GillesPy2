#ifndef linked_list_h
#define linked_list_h

// Create data structure for a node of the list

typedef struct node_t {
	int data;
	struct node_t * next;
}node;

typedef struct node2_t {
	int val;
	int count;
	struct node2_t * next;
}node2;

// Data structure for the linked list type

typedef struct linked_list_t{
	node*head;
	node*tail;
	node2*head2;
	node2*tail2;
	int count;
	int count2;
}linked_list;

// Functions to manipulate the linked list

// constructor
linked_list* create_linked_list();

// destructor
void destroy_linked_list( linked_list* ll );

// add a new node to the end of the list
void linked_list_add( linked_list* ll, int data_in);

//Generate Histogram
void linked_list_update_histogram( linked_list* ll);

void add_to_histogram( linked_list* ll, int to_add_value);
// delete a node from the linked list
void linked_list_delete( linked_list* ll, node* to_delete);
void linked_list_delete2( linked_list* ll, node2* to_delete);
void linked_list_clean( linked_list* ll);

// Insertion sort ll
void insertion_sort(node2 **head, node2 **tail);

#endif /* linked_list_h*/
