#ifndef linked_list_h
#define linked_list_h

// Create data structure for a node of the list

typedef struct node_t {
	int data;
	struct node_t * next;
}node;

// Data structure for the linked list type

typedef struct linked_list_t{
	node*head;
	node*tail;
	int count;
}linked_list;

// Functions to manipulate the linked list

// constructor
linked_list* create_linked_list();

// destructor
void destroy_linked_list( linked_list* ll );

// add a new node to the end of the list
void linked_list_add( linked_list* ll, int data_in);

// delete a node from the linked list
void linked_list_delete( linked_list* ll, node* to_delete);

// search for a node by it's data field
node* linked_list_search( linked_list* ll, int search_value );

// get node by index
node* linked_list_get( linked_list* ll, int index);


#endif /* linked_list_h*/
