#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <pthread.h>
#include <math.h>
#include <time.h>
#include "linked_list.h"

#define ALPHA 0.05
#define STR_TO_INT_BUF 8
#define READ_BUFFER_SIZE 512
#define BASE_TRAJECTORIES 100
#define TEST_OFFSET 312
pthread_mutex_t evens_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t odds_mutex = PTHREAD_MUTEX_INITIALIZER;

//Thread args
struct arg {
	int p_num;
	int num_runs;
	linked_list* arr;
	char** argv;
	pthread_mutex_t *mutex;
};
int* species, *num_lls;

void* c_solver_runner(void *targ_in);
linked_list* access_array(linked_list* ll_ptr, int x, int y, int num_species);
//void* sort_histogram( void* ll);
double calculate_ks_distance(linked_list *evens, linked_list *odds, int cycle);
void plot_bounds(linked_list* evens, linked_list* odds, size_t array_size);

int main (int argc, char** argv){

	int processes, num_runs, num_species, i, num_timesteps;
	double max_dist, elapsed;
	processes = atoi(argv[2]);
	num_species = atoi(argv[3])+1;
	species = &num_species;
	num_timesteps = atoi(argv[4]);
	pthread_t *thread_handle = malloc(sizeof(pthread_t)*processes);
	struct arg *targ = malloc(sizeof(struct arg)*processes);
	num_runs = BASE_TRAJECTORIES / processes;
	int run_count = BASE_TRAJECTORIES;
	int cycle = 1;
	time_t beg, mid, end;
	beg = clock();	

	size_t array_size = num_species * num_timesteps;
	num_lls = (int*) &array_size;
	printf("NUM LLS: %i\n", *num_lls);

	linked_list *evens = calloc(array_size, sizeof(linked_list));
	linked_list *odds = calloc(array_size, sizeof(linked_list));
	

	if (argc < 3){
		printf("Insufficient arguments\n");
		exit(1);
	}

	if (processes == 1)
		processes = 2;

	/*
	   Create threads and assign args
	 */
	for(i = 0; i < processes; i++){
		targ[i].p_num = i+1;
		if ( targ[i].p_num % 2 == 0){
			targ[i].arr = evens;
			targ[i].mutex = &evens_mutex;
		}
		else{
			targ[i].arr = odds;
			targ[i].mutex = &odds_mutex;
		}
		targ[i].argv = argv;
		targ[i].num_runs = num_runs;
		if (BASE_TRAJECTORIES % processes >= targ[i].p_num){
			targ[i].num_runs++;
		}
		pthread_create(&thread_handle[i], NULL, c_solver_runner, &targ[i]);
		//sleep(1);
	}

	for(int i = 0; i < processes; i++){
		pthread_join(thread_handle[i], NULL);
	}
	/**
	pthread_t evens_handler, odds_handler;
	pthread_create(&evens_handler, NULL, sort_histogram, &targ[1]);
	pthread_create(&odds_handler, NULL, sort_histogram, &targ[0]);

	pthread_join(evens_handler, NULL);
	pthread_join(odds_handler, NULL);
	**/
	printf("\n***TEST SUITE***\n");
	linked_list* tptr = odds;
	tptr += TEST_OFFSET;
	//for (int i = 0; i < 50; i++){
	/**
	printf("ODDS first ll size: %i\n", tptr->count2);
	node *nptr = tptr->head;
	printf("Elements in linked list:\n");
	while(nptr != NULL){
		printf("%i, ", nptr->data);
		nptr=nptr->next;
	}
	printf("\n");**/
	node2 *n2ptr = tptr->head2;
	printf("ODDS @ +%i location, # elements: %i\n", TEST_OFFSET, tptr->count);
	printf("Elements in linked list Histogram:\n");
	while(n2ptr != NULL){
		printf("[value: %i, count: %i] ,", n2ptr->val, n2ptr->count);
		n2ptr=n2ptr->next;
	}
	printf("\n");

	tptr = evens;
	tptr += TEST_OFFSET;
	printf("\nEVENS @ +%i location, # elements: %i\n", TEST_OFFSET, tptr->count);
	//for (int i = 0; i < 50; i++){
	/**
	printf("ODDS first ll size: %i\n", tptr->count2);
	node *nptr = tptr->head;
	printf("Elements in linked list:\n");
	while(nptr != NULL){
		printf("%i, ", nptr->data);
		nptr=nptr->next;
	}
	printf("\n");
	**/
	n2ptr = tptr->head2;
	printf("Elements in linked list Histogram:\n");
	while(n2ptr != NULL){
		printf("[value: %i, count: %i] ,", n2ptr->val, n2ptr->count);
		n2ptr=n2ptr->next;
	}
	printf("\n");
	//tptr++;
	//}
	
	/*
	*FIRST CHECK KS DIST
	*/
	max_dist = calculate_ks_distance(evens, odds, cycle);
	mid = clock();
	elapsed =(double) (mid-beg) / (double) CLOCKS_PER_SEC;
	
	printf("Run Count: %i, Max distance = %f\n", run_count, max_dist);
	printf("Time Elapsed This Cycle: %f seconds.\n", elapsed);

	
	/*
	*Repeat until fit reached
	*/
	while (max_dist > ALPHA){
		max_dist = 0;
		cycle++;
		for (i = 0; i < processes; i++){
			pthread_create(&thread_handle[i], NULL, c_solver_runner, &targ[i]);
			//sleep(1);
		}
	for(int i = 0; i < processes; i++){
		pthread_join(thread_handle[i], NULL);
	}
	printf("\n***TEST SUITE***\n");
	linked_list* tptr = odds;
	tptr += TEST_OFFSET;
	printf("ODDS @ +%i location, head: %i: %i tail %i: %i, # elements: %i\n", TEST_OFFSET, tptr->head2->val, tptr->head2->count, tptr->tail2->val,tptr->tail2->count, tptr->count);
	//for (int i = 0; i < 50; i++){
	/**
	printf("ODDS first ll size: %i\n", tptr->count2);
	node *nptr = tptr->head;
	printf("Elements in linked list:\n");
	while(nptr != NULL){
		printf("%i, ", nptr->data);
		nptr=nptr->next;
	}
	printf("\n");
	**/
	node2 *n2ptr = tptr->head2;
	printf("Elements in linked list Histogram:\n");
	while(n2ptr != NULL){
		printf("[value: %i, count: %i] ,", n2ptr->val, n2ptr->count);
		n2ptr=n2ptr->next;
	}
	printf("\n");

	tptr = evens;
	tptr += TEST_OFFSET;
	printf("\nEVENS @ +%i location, # elements: %i\n", TEST_OFFSET, tptr->count);
	//for (int i = 0; i < 50; i++){
	/**
	printf("ODDS first ll size: %i\n", tptr->count2);
	node *nptr = tptr->head;
	printf("Elements in linked list:\n");
	while(nptr != NULL){
		printf("%i, ", nptr->data);
		nptr=nptr->next;
	}
	printf("\n");
	**/
	n2ptr = tptr->head2;
	printf("Elements in linked list Histogram:\n");
	while(n2ptr != NULL){
		printf("[value: %i, count: %i] ,", n2ptr->val, n2ptr->count);
		n2ptr=n2ptr->next;
	}
	printf("\n");
	//tptr++;
	//}
	run_count+=BASE_TRAJECTORIES;
	max_dist = calculate_ks_distance(evens, odds, cycle);
	end = clock();
	elapsed =(double) (end - mid) / (double) CLOCKS_PER_SEC;
	mid = end;
	printf("Run Count: %i, Max distance = %f\n", run_count, max_dist);
	printf("Time Elapsed This Cycle: %f seconds.\n", elapsed);
	}


	max_dist = calculate_ks_distance(evens, odds, cycle);

	plot_bounds(evens, odds, array_size);
	printf("Run Count: %i, Max distance = %f\n", run_count, max_dist);
	end = clock();
	elapsed =(double) (end-beg) / (double) CLOCKS_PER_SEC;
	printf("TOTAL Time Elapsed: %f seconds.\n", elapsed);
	free(evens);
	free(odds);
	free(thread_handle);
	free(targ);
	return 0;
}
/*
 * Thread Function
 */
void *c_solver_runner(void *targ_in){

	int i;
	int filedes[2];
	pid_t pid;

	struct arg* targ =  (struct arg*) targ_in;
	//printf("p_num: %i, num_runs: %i\n", targ->p_num, targ->num_runs);
	for (i = 0; i < targ->num_runs; i++){
		/*
		   Create child processes from threads to exec, redirect stdout from executable to pipe and close it on child end
		 */
		if(pipe(filedes) == -1){
			perror("pipe");
			exit(1);
		}
		if ((pid = fork()) == 0){       //child
			char seed[10];
			sprintf(seed, "%i", getpid());
			char* args[] = {targ->argv[1], "-seed", seed};
			//printf("Seed: %s\n", seed);
			while ((dup2(filedes[1], STDOUT_FILENO) == -1) && (errno == EINTR)){}
			close(filedes[1]);
			close(filedes[0]);
			//printf("test\n");
			execv(targ->argv[1], args);
			perror("execv");

			_exit(1);
		}else{
			/*
			   Close Pipe output in thread and read pipe data to buffer
			   --*/
			close(filedes[1]);
			char buffer[READ_BUFFER_SIZE];
			int current_index = 0;
			int current_row = 0;
			int current_col = 0;
			linked_list *writer;
			writer = targ->arr;
			char builder[STR_TO_INT_BUF];
			memset(builder, 0, STR_TO_INT_BUF*sizeof(char));
			while (1){
				ssize_t count = read(filedes[0], buffer, sizeof(buffer));
				if (count == -1){
					if (errno == EINTR){
						continue;
					} else {
						perror("read");
						exit(1);
					}
				} else if (count == 0){
					break;
				} else {
					/*
					   Parse the buffer
					 */
					for (int i=0; i < READ_BUFFER_SIZE; i++){
						if (buffer[i] == '\n'){
							if(buffer[i+1] == '\n')
								break;
							current_index = 0;
							current_row++;
							current_col = 0;
							writer = access_array(targ->arr, 0, current_row, *species);
						}else if(buffer[i] == ' ') {
							pthread_mutex_lock(targ->mutex);
							linked_list_add(writer, atoi(builder));
							//printf("PTHREAD: %i LL ADDR: %p LL LOC: (%i, %i) LL Count: %i HEAD is: %i TAIL is: %i\n", targ->p_num, writer, current_col, current_row, writer->count, writer->head->data, writer->tail->data);
//fflush(NULL);
							pthread_mutex_unlock(targ->mutex);
							writer++;
							current_col++;
							memset(builder, 0, STR_TO_INT_BUF*sizeof(char));
							current_index = 0;
						}else if (buffer[i] > 0){
							builder[current_index] = buffer[i];
							current_index++;
						}
					}

				}
			}
		}
		close(filedes[0]);
	}
	pthread_exit(0);
}


linked_list* access_array(linked_list *ll_ptr, int x, int y, int num_species){
	return ll_ptr+x+(y*num_species);
}
/**
void* sort_histogram(void *targ_in){
	struct arg *targ = (struct arg*) targ_in;
	linked_list* ll_ptr = targ->arr;
	node* t = ll_ptr->head;
	while (ll_ptr->head != NULL){
		insertion_sort(&ll_ptr->head2, &ll_ptr->tail2);
		//linked_list_clean(ll_ptr);
		ll_ptr++;
	}
		
}
**/
double calculate_ks_distance(linked_list *evens, linked_list *odds, int cycle){

	linked_list* evens_ptr = evens;
	linked_list* odds_ptr = odds;
	int e_num_traj = evens->count * cycle;
	int o_num_traj = odds->count * cycle;
	double e_norm, o_norm;
	double max_distance = 0;
	double diff;
	int comparisons = 0;
	node2 *e, *o, *pe, *po;
	int evens_count = 0;
	int odds_count = 0;
	int number_of_linked_lists = 0;
	for(int i = 0; i < *num_lls; i++){
		//printf("evens_ptr->count2: %i\n", evens_ptr->count2);
		//printf("evens head2 val: %i count: %i\n", evens_ptr->head2->val, evens_ptr->head2->count);
		//printf("evens_ptr->count\n");
		//printf("HEAD2 VAL: %i, COUNT: %i\n", evens_ptr->head2->val, evens_ptr->head2->count);
		//printf("HEAD: %p, HEAD2: %p\n", evens_ptr->head, evens_ptr->head2);
		/** PRINT STUFF
		node2 *et, *ot;
		et = evens_ptr->head2;
		ot = odds_ptr->head2;
		printf("EVENS for this LL: ");
		while(et != NULL){
			printf("%i, ", et->val);
			et = et->next;
		}
		printf("\nODDS for this LL: ");
		while(ot != NULL){
			printf("%i, ", ot->val);
			ot = ot->next;
		}
		printf("\n");
		**/
		e = evens_ptr->head2;
		o = odds_ptr->head2;
		int iteration = 0;
		//RUN for head special case)
		if (e->val == o->val){
			e_norm = (double)e->count / e_num_traj;
			o_norm = (double)o->count / o_num_traj;
			diff = fabs(e_norm - o_norm);
			comparisons++;
			if (diff > max_distance){
				max_distance = diff;
			}
		}else if(e->val < o->val){
			e_norm = (double)e->count / e_num_traj;
			diff = e_norm; 
			comparisons++;
			if (diff > max_distance){
				max_distance = diff;
			}
		}else{
			o_norm = (double)o->count / o_num_traj;
			diff = o_norm;
			comparisons++;
			if (diff > max_distance){
				max_distance = diff;
			}
		}
		iteration++;
		pe = e;
		e = e->next;
		po = o;
		o = o->next;
		while(e != NULL && o != NULL){
			if (e->val == o->val){
				e_norm = (double)(e->count+pe->count) / e_num_traj;
				o_norm = (double)(o->count+po->count) / o_num_traj;
				diff = fabs(e_norm - o_norm);
				//printf("DIFF: %f\n", diff);
				comparisons++;
				if (diff > max_distance) max_distance = diff;
				e = e->next;
				o = o->next;
			}else if(e->val < o->val){
				e_norm = (double)e->count / e_num_traj;
				diff = e_norm; 
				//printf("DIFF: %f\n", diff);
				comparisons++;
				if (diff > max_distance) max_distance = diff;
				e = e->next;
			}else{
				o_norm = (double)o->count / o_num_traj;
				diff = o_norm;
				//printf("DIFF: %f\n", diff);
				comparisons++;
				if (diff > max_distance) max_distance = diff;
				o = o->next;
			}
			iteration++;
		}

		//printf("ITERATIONS IN THIS LL: %i\n", iteration);
		evens_ptr++;
		odds_ptr++;
		number_of_linked_lists++;
		//printf("Number of Linked Lists: %i, Number of Comparisons: %i\n", number_of_linked_lists, comparisons);
	}
	//printf("Comparisons at end of array: %i\n", comparisons);
	return max_distance;
}

void plot_bounds(linked_list* evens, linked_list* odds, size_t array_size){
	linked_list *ptr;
	ptr = evens;
	printf("\n\nUPPER AND LOWER BOUNDS:\nEVEN\n");
	for (int i = 0; i < array_size; i++){
		if(i%*species == 0)
			printf("%i: ", ptr->head2->val);
		else if(i%*species == *species-1)
			printf(" [L: %i, U: %i]\n", ptr->head2->val, ptr->tail2->val);
		else
			printf("[L: %i, U: %i], ", ptr->head2->val, ptr->tail2->val);
		ptr++;	
	}
	ptr = odds;
	printf("ODDS\n");
	for (int i = 0; i < array_size; i++){
		if(i%*species == 0)
			printf("%i: ", ptr->head2->val);
		else if(i%*species == *species-1)
			printf(" [L: %i, U: %i]\n", ptr->head2->val, ptr->tail2->val);
		else
			printf("[L: %i, U: %i], ", ptr->head2->val, ptr->tail2->val);
		ptr++;	
	}

}
