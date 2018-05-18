#include <string.h>  // for strcomp and strcpy
#include <stdlib.h> // for malloc and free
#include <stdio.h> // for printf
#include <unistd.h> // for testing with sleep
#include "linked_list.h"


void sorted_insert(node2** head, node2* new_node);

// constructor
linked_list* create_linked_list(){
	linked_list* ll = (linked_list*) malloc( sizeof(linked_list));
	ll->count = 0;
	ll->count2 = 0;
	ll->head = NULL;
	ll->tail = NULL;
	ll->head2 = NULL;
	ll->tail2 = NULL;
	return ll;
}

// destructor
void destroy_linked_list( linked_list* ll ){
	// empty the linked list
	while( ll->count > 0){
		linked_list_delete( ll, ll->head );
	}
	while( ll->count2 > 0){
		linked_list_delete2( ll, ll->head2 );
	}
// un-allocate the memory
	free(ll);
}

void linked_list_clean( linked_list* ll){
	//clear raw data list
	while(ll->count > 0){
		linked_list_delete( ll, ll->head );
	}
	ll->head = NULL;
	ll->tail = NULL;
}


// add a new node to the end of the list
void linked_list_add( linked_list* ll, int data_in){
	node* n = (node *) malloc( sizeof(node) );
	n->data = data_in;
	n->next = NULL;
	// Travesere the list to find the end note
	node * t = ll->head;
	if(t == NULL){
		ll->head = n;
		ll->tail = n;
	}else if(t->next == NULL){
		t->next = n;
		ll->tail = n;
	}else{
		ll->tail->next = n;
		ll->tail = n;
	}
	// increase the size
	ll->count++;
	add_to_histogram(ll, data_in);
}

void add_to_histogram( linked_list* ll, int to_add_value ){

	node2* n;
	node2* t = ll->head2;
	n = (node2*) malloc(sizeof(node2));
	n->val = to_add_value;
	n->count = 1;
	n->next = NULL;

	if( t == NULL){
		ll->head2 = n;
		ll->tail2 = n;
		ll->count2++;
		return;
	}else if(t->val != n->val && t->next == NULL){
		t->next = n;
		ll->tail2 = n;
		ll->count2++;
		return;
	}else{
		while(t != NULL){
			if(t->val == n->val){
				t->count++;
				return;
			}
			t = t->next;
		}
		ll->tail2->next = n;
		ll->tail2 = n;
		ll->count2++;
	}
}


// delete a node from the linked list
void linked_list_delete( linked_list* ll, node* to_delete){
	node* prev_node = ll->head;
	if (ll->head == NULL){
		printf("Error, linked_list_delete() empty list\n");
		return;
	}else if (to_delete == ll->head ){
		ll->head = ll->head->next;
	}else{
		for( prev_node = ll->head; prev_node->next!=NULL; prev_node=prev_node->next ){
			if(prev_node->next == to_delete){
				break;
			}
		}
		if( prev_node->next == NULL){
			printf("Error, linked_list_delete(), could not find item in list \n");
			return;
		}
		prev_node->next = to_delete->next; //connect the list
	}

	//free and reduce size
	ll->count--;
	//free(to_delete->data);
	free(to_delete);


}
void linked_list_delete2( linked_list* ll, node2* to_delete){
	node2* prev_node = ll->head2;
	if (ll->head2 == NULL){
		printf("Error, linked_list_delete() empty list\n");
		return;
	}else if (to_delete == ll->head2 ){
		ll->head2 = ll->head2->next;
	}else{
		for( prev_node = ll->head2; prev_node->next!=NULL; prev_node=prev_node->next ){
			if(prev_node->next == to_delete){
				break;
			}
		}
		if( prev_node->next == NULL){
			printf("Error, linked_list_delete(), could not find item in list \n");
			return;
		}
		prev_node->next = to_delete->next; //connect the list
	}

	//free and reduce size
	ll->count2--;
	//free(to_delete->data);
	free(to_delete);


}

//sort sourced and modified from https://www.geeksforgeeks.org/insertion-sort-for-singly-linked-list/
void insertion_sort(node2 **head)
{
	// Initialize sorted linked list
	node2 *sorted = NULL;

	// Traverse the given linked list and insert every
	// node to sorted
	node2 *current = *head;
	while (current != NULL)
	{
		// Store next for next iteration
		node2 *next = current->next;

		// insert current in sorted linked list
		sorted_insert(&sorted, current);

		// Update current
		current = next;
	}

	// Update head_ref to point to sorted linked list
	*head = sorted;
}

//Helper function for insertion_sort
void sorted_insert(node2** head, node2* new_node)
{
	node2* current;
	/* Special case for the head end */
	if (*head == NULL || (*head)->val >= new_node->val)
	{
		new_node->next = *head;
		*head = new_node;
	}
	else
	{
		/* Locate the node before the point of insertion */
		current = *head;
		while (current->next!=NULL &&
				current->next->val < new_node->val)
		{
			current = current->next;
		}
		new_node->next = current->next;
		current->next = new_node;
	}
}


