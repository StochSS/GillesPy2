#include <string.h>  // for strcomp and strcpy
#include <stdlib.h> // for malloc and free
#include <stdio.h> // for printf
#include <unistd.h> // for testing with sleep
#include "linked_list.h"


void sorted_insert(linked_list *ll, node2* new_node);
void binary_insert(linked_list *ll, node2* new_node);
node2* middle(node2* start, node2* last);

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
	while( ll->head != 0){
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
	while(ll->head != NULL){
		linked_list_delete( ll, ll->head );
	}
	ll->head = NULL;
	ll->tail = NULL;
}


// add a new node to the end of the list
void linked_list_add( linked_list* ll, int data_in){
	/*
	node* n = (node *) malloc( sizeof(node) );
	n->data = data_in;
	n->next = NULL;
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
	*/
	add_to_histogram(ll, data_in);
	ll->count++;
}

void add_to_histogram( linked_list* ll, int to_add_value ){

	node2* n;
	node2* t = ll->head2;
	n = malloc(sizeof(node2));
	n->val = to_add_value;
	n->count = 1;
	n->next = NULL;

	//Try to add to LL at O(1) if possible.
	if( t == NULL){
		ll->head2 = n;
		ll->tail2 = n;
		ll->count2++;
		return;
	}else if(t->val < n->val && t->next == NULL){
		t->next = n;
		ll->tail2 = n;
		ll->count2++;
		return;
	}else if(t->val == n->val){
		t->count++;
		return;
	}else if(ll->tail2->val == n->val){
		ll->tail2->count++;
		return;
	}else if(n->val < ll->head2->val){
		n->next = ll->head2;
		ll->head2 = n;
		t = n;
		ll->count2++;
		return;
	}else if (n->val > ll->tail2->val){
		ll->tail2->next = n;
		ll->tail2 = n;
		ll->count2++;
		return;
	}

		//if (ll->count2 < 10){
			sorted_insert(ll, n);
		//}else{
		//	binary_insert(ll, n);
		//}
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
		//ll->count--;
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
	//Helper function for insertion_sort
	void sorted_insert(linked_list *ll, node2* new_node)
	{
		node2* current;
		/* Locate the node before the point of insertion */
		current = ll->head2;
		while (current->next!=NULL && current->next->val < new_node->val){
			current = current->next;
		}
		if(current->next != NULL){
			//printf("current->val: %i, current->next->val: %i, new_node->val: %i\n", current->val, current->next->val, new_node->val);
			if (current->next->val == new_node->val){
				//printf("incrementing a value\n");
				current->count++;
				ll->count2++;
				return;
			}
		}
		new_node->next = current->next;
		current->next = new_node;
		//ll->count2++;
		if(new_node->next == NULL)
			ll->tail2 = new_node;

	}

node2* middle(node2* start, node2* last)
{
	if (start == NULL)
		return NULL;

	node2* slow = start;
	node2* fast = start -> next;

	while (fast != last)
	{
		fast = fast -> next;
		if (fast != last)
		{
			slow = slow -> next;
			fast = fast -> next;
		}
	}

	return slow;
}

// Function for implementing the Binary
// Search on linked list
void binary_insert(linked_list *ll, node2 *new_node)
{
	node2* start = ll->head2;
	node2* last = ll->tail2;
	node2* mid;

	do
	{
		if(start->val == new_node->val){
			start->count++;
			//printf("START: %i, INSERT: %i\n", start->val, new_node->val);
			//printf("INCREMENT START\n");
			return;
		}else if(last->val == new_node->val){
			last->count++;
			printf("LAST: %i, INSERT: %i\n", last->val, new_node->val);
			printf("INCREMENT LAST\n");
			return;
		}
		// Find middle
		mid = middle(start, last);
		//printf("START: %i MID: %i LAST: %i INSERT: %i MID->next: %i\n", start->val, mid->val, last->val, new_node->val, mid->next->val);

		// If middle is empty
		if (mid == NULL){
			printf("ERROR: NULL MID\n");
			return;
		}

		// If mid contains to-add
		if (mid -> val == new_node->val){
			mid->count++;
			//printf("INCREMENTED %i\n", new_node->val);
			return;
		}

		if(start == mid){
			if (start->next->val == new_node->val){
				start->next->count++;
				//printf("Incrementing start->next\n");
			}else{
			//printf("start= mid adding between %i and %i\n", start->val, start->next->val);
			new_node->next = start->next;
			start->next = new_node;
			ll->count2++;
			}
			return;
		}

		// If value is more than mid
		else if (mid->val < new_node->val){
			if(new_node->val < mid->next->val){
				//printf("Placing between %i and %i\n", mid->val, mid->next->val);
				new_node->next = mid->next;
				mid->next = new_node;
				ll->count2++;
				return;
			}
			start = mid->next;
		}

		// If value is less than mid.
		else
			last = mid;


	} while (last == NULL || last -> next != start);
	printf("Exit loop\n");
	return;
}

