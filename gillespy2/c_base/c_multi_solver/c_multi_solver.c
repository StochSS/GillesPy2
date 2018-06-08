#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <pthread.h>
#include <math.h>
#include <time.h>
#include "linked_list.h"

#define STR_TO_INT_BUF 8
#define READ_BUFFER_SIZE 512
#define BASE_TRAJECTORIES 100
#define TEST_OFFSET 312
#define EVENS 1
#define ODDS 0

pthread_mutex_t evens_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t odds_mutex = PTHREAD_MUTEX_INITIALIZER;

//Thread args
struct arg {
	int p_num;
	int num_runs;
	linked_list* arr;
	char* exec;
	pthread_mutex_t *mutex;
	pthread_mutex_t *mutex_array;
};

//helper methods
void* c_solver_runner(void *targ_in);
void *update_histogram(void *targ_in);
linked_list* access_array(linked_list* ll_ptr, int x, int y);
double calculate_ks_distance(linked_list *evens, linked_list *odds);
void plot_bounds(linked_list* evens, linked_list* odds, size_t array_size);
void parse(int pipe, struct arg* targ);
void parse_binary(int pipe, struct arg* targ);
void print_status(linked_list *ll, int which);


//global variables from command line args
int *species, *timesteps, *num_lls, *num_processes;
double *end_time, *alpha;

int main (int argc, char** argv){

	int processes, num_runs, num_species, i, num_timesteps, run_count;
	double max_dist, end_t, alpha_value;
	char* executable;
	size_t array_size;
	//time_t beg, mid, end;
	struct arg *targ;
	linked_list *evens, *odds;

	//handle command line args
	executable = argv[1];
	processes = atoi(argv[2]);
	num_species = atoi(argv[3])+1;
	num_timesteps = atoi(argv[4]);
	end_t = atof(argv[5]);
	alpha_value = atof(argv[6]);
	array_size = num_species * num_timesteps;

	//make args global
	species = &num_species;
	timesteps = &num_timesteps;
	end_time = &end_t;
	alpha = &alpha_value;
	num_lls = (int*) &array_size;
	num_processes = &processes;

	//allocate thread memory
	pthread_t *thread_handle = malloc(sizeof(pthread_t)*processes);
	targ = malloc(sizeof(struct arg)*processes);

	//set initial states
	num_runs = BASE_TRAJECTORIES / processes;
	run_count = BASE_TRAJECTORIES;
	//beg = clock();	

	//allocate data memory
	evens = calloc(array_size, sizeof(linked_list));
	odds = calloc(array_size, sizeof(linked_list));

	// create a mutex for each linked list
	pthread_mutex_t evens_mutex_array[array_size];
	pthread_mutex_t odds_mutex_array[array_size];

	for(int i = 0; i < array_size; i++){

		pthread_mutex_init(&evens_mutex_array[i], NULL);
		pthread_mutex_init(&odds_mutex_array[i], NULL);
	}
	

	if (argc < 5){
		printf("Insufficient arguments\n");
		exit(1);
	}

	//need at least 2 processes for evens and odds
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
			targ[i].mutex_array = evens_mutex_array;
		}
		else{
			targ[i].arr = odds;
			targ[i].mutex = &odds_mutex;
			targ[i].mutex_array = odds_mutex_array;
		}
		targ[i].exec = executable;
		targ[i].num_runs = num_runs;
		if (BASE_TRAJECTORIES % processes >= targ[i].p_num){
			targ[i].num_runs++;
		}
		pthread_create(&thread_handle[i], NULL, c_solver_runner, &targ[i]);
	}

	for(i = 0; i < processes; i++){
		pthread_join(thread_handle[i], NULL);
	}

	for(i = 0; i < 2; i++){
		pthread_create(&thread_handle[i], NULL, update_histogram, &targ[i]);
	}
	for(i = 0; i < 2; i++){
		pthread_join(thread_handle[i], NULL);
	}

	
	//printf("Exited sort\n");

	//TEST PRINTING
	//print_status(odds, ODDS);
	//print_status(evens, EVENS);
	//fflush(NULL);

	/*
	*FIRST CHECK KS DIST
	*/
	max_dist = calculate_ks_distance(evens, odds);
	//printf("Calculated distance\n");
	//fflush(NULL);
	//calculate elapsed time for dev purposes
	//mid = clock();
	//elapsed =(double) (mid-beg) / (double) CLOCKS_PER_SEC;
	
	//print stats
	//printf("Run Count: %i, Max distance = %f\n", run_count, max_dist);
	//printf("Time Elapsed This Cycle: %f seconds.\n", elapsed);


	/*
	 *Repeat until fit reached
	 */
	while (max_dist > *alpha){
		max_dist = 0;
		for (i = 0; i < processes; i++){
			pthread_create(&thread_handle[i], NULL, c_solver_runner, &targ[i]);
		}
		for(int i = 0; i < processes; i++){
			pthread_join(thread_handle[i], NULL);
		}
		for(i = 0; i < 2; i++){
			pthread_create(&thread_handle[i], NULL, update_histogram, &targ[i]);
		}
		for(i = 0; i < 2; i++){
			pthread_join(thread_handle[i], NULL);
		}

		//TEST PRINTING
		//print_status(odds, ODDS);
		//print_status(evens, EVENS);

		//track total runs and update dist
		run_count+=BASE_TRAJECTORIES;
		max_dist = calculate_ks_distance(evens, odds);

		//calculate elapsed time for dev purposes
		
		   //end = clock();
		   //elapsed =(double) (end - mid) / (double) CLOCKS_PER_SEC;
		   //mid = end;
		 
		//print stats
		//printf("Run Count: %i, Max distance = %f\n", run_count, max_dist);
		//printf("Time Elapsed This Cycle: %f seconds.\n", elapsed);
	}



	//print final data
	//plot_bounds(evens, odds, array_size);
	//printf("Run Count: %i, Max distance = %f\n", run_count, max_dist);

	//print total time
	
	   //end = clock();
	   //elapsed =(double) (end-beg) / (double) CLOCKS_PER_SEC;
	   //printf("TOTAL Time Elapsed: %f seconds.\n", elapsed);
	 
	printf("\nTOTAL Runs: %i\n", run_count);

	//be responsible
	free(evens);
	free(odds);
	free(thread_handle);
	free(targ);
	return 0;
}
/*
 * Thread Functions
 */

void *update_histogram(void *targ_in){
	
	struct arg* targ = (struct arg*) targ_in;
	//time_t beg, end;
	//double elapsed;

	//beg = clock();
	/*I
	int sorter_number;
	int num_even_processes = *num_processes / 2;
	int num_odd_processes = (*num_processes / 2) + (*num_processes % 2);

	if (targ->p_num % 2)
		sorter_number = (targ->p_num / 2) + 1;
	else
		sorter_number = targ->p_num / 2;
	int sort_size;
	int offset;
	if (targ->p_num % 2 == 1){
		sort_size = *timesteps / num_odd_processes;
		offset = 0;
	}else{
		sort_size = *timesteps / num_even_processes;
		offset = 1;
	}

	if (*timesteps % *num_processes >= sorter_number)
		sort_size++;


	int start_location = (targ->p_num/2 - offset)*sort_size;
	int num_lists_to_sort = sort_size * *species;
	*/
	//access beginning of this threads sort locations
	//linked_list *sorter = access_array(targ->arr, 0, start_location);
	//printf("beginning sort at %i of array at %p\n", start_location, targ->arr);
	linked_list *sorter = targ->arr;
	for (int i=0; i < *num_lls; i++){
		build_histogram(sorter);
		sorter++;
	}
/*
	for (int i=0; i < num_lists_to_sort; i++){
		build_histogram(sorter);
		sorter++;
	}
*/
	//end = clock();
	//elapsed =(double) (end-beg) / (double) CLOCKS_PER_SEC;
	//printf("P_NUM %i: HISTOGRAM UPDATE TIME: %f\n",targ->p_num, elapsed);

}
void *c_solver_runner(void *targ_in){
	
	int i;
	int filedes[2];
	pid_t pid;
	struct arg* targ =  (struct arg*) targ_in;
		/*
		   Create child processes from threads to exec, redirect stdout from executable to pipe and close it on child end
		 */
		if(pipe(filedes) == -1){
			perror("pipe");
			exit(1);
		}
		if ((pid = vfork()) == 0){       //child
			char trajectories[10], steps[10], end[10], seed[10];
			sprintf(seed, "%i", getpid());
			sprintf(steps, "%i", *timesteps);
			sprintf(end, "%f", *end_time);
			sprintf(trajectories, "%i", targ->num_runs);
			char* args[] = {targ->exec, "-trajectories", trajectories, "-timesteps", steps, "-end", end, "-seed", seed, NULL};
			//printf("Seed: %s\n", seed);/
			while ((dup2(filedes[1], STDOUT_FILENO) == -1) && (errno == EINTR)){}
			close(filedes[1]);
			close(filedes[0]);
			//printf("test\n");
			execv(targ->exec, args);
			perror("execv");

			_exit(1);
		}else{
			/*
			   Close Pipe output in thread and read pipe data to buffer
			   --*/
			close(filedes[1]);
			parse_binary(filedes[0], targ);
		}
	close(filedes[0]);
	pthread_exit(0);
}


linked_list* access_array(linked_list *ll_ptr, int x, int y){
	return ll_ptr+x+(y* *species);
}


double calculate_ks_distance(linked_list *evens, linked_list *odds){

	//initialize everythin for beginning of array (0th linked list)
	linked_list* evens_ptr = evens;
	linked_list* odds_ptr = odds;
	int e_num_traj = (evens+1)->count;
	int o_num_traj = (odds+1)->count;
	//printf("e_num_traj = %i\n", e_num_traj);
	double e_norm, o_norm;
	double max_distance = 0;
	double diff;
	node2 *e, *o, *pe, *po;
	int evens_count = 0;
	int odds_count = 0;

	//iterate through linked lists
	for(int i = 0; i < *num_lls; i++){
		//iterate through histogram in linked list
		e = evens_ptr->head2;
		o = odds_ptr->head2;

		//RUN for head special case
		if (e->val == o->val){
			e_norm = (double)e->count / e_num_traj;
			o_norm = (double)o->count / o_num_traj;
			diff = fabs(e_norm - o_norm);
			if (diff > max_distance){
				max_distance = diff;
			}
		}else if(e->val < o->val){
			e_norm = (double)e->count / e_num_traj;
			diff = e_norm; 
			if (diff > max_distance){
				max_distance = diff;
			}
		}else{
			o_norm = (double)o->count / o_num_traj;
			diff = o_norm;
			if (diff > max_distance){
				max_distance = diff;
			}
		}

		//increment histogram ptrs and track previous
		pe = e;
		e = e->next;
		po = o;
		o = o->next;
		while(e != NULL && o != NULL){
			if (e->val == o->val){
				e_norm = (double)(e->count+pe->count) / e_num_traj;
				o_norm = (double)(o->count+po->count) / o_num_traj;
				diff = fabs(e_norm - o_norm);
				if (diff > max_distance) max_distance = diff;
				e = e->next;
				o = o->next;
			}else if(e->val < o->val){
				e_norm = (double)e->count / e_num_traj;
				diff = e_norm; 
				if (diff > max_distance) max_distance = diff;
				e = e->next;
			}else{
				o_norm = (double)o->count / o_num_traj;
				diff = o_norm;
				if (diff > max_distance) max_distance = diff;
				o = o->next;
			}
		}

		//move to next linked list
		evens_ptr++;
		odds_ptr++;
	}

	return max_distance;
}

void plot_bounds(linked_list* evens, linked_list* odds, size_t array_size){
	linked_list *ptr;
	ptr = evens;
	printf("\n\nUPPER AND LOWER BOUNDS:\nEVEN\n");
	for (int i = 0; i < array_size; i++){
		if(i%*species == 0)
			printf("%i:%f: ", ptr->head2->val, ptr->timestep);
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
			printf("%i:%f ", ptr->head2->val, ptr->timestep);
		else if(i%*species == *species-1)
			printf(" [L: %i, U: %i]\n", ptr->head2->val, ptr->tail2->val);
		else
			printf("[L: %i, U: %i], ", ptr->head2->val, ptr->tail2->val);
		ptr++;	
	}

}

void parse(int pipe, struct arg* targ){
	char buffer[READ_BUFFER_SIZE];
	int current_index = 0;
	int current_row = 0;
	int current_col = 0;
	linked_list *writer;
	writer = targ->arr;
	char builder[STR_TO_INT_BUF];
	memset(builder, 0, STR_TO_INT_BUF*sizeof(char));
	while (1){
		ssize_t count = read(pipe, buffer, sizeof(buffer));
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
					writer = access_array(targ->arr, 0, current_row);
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

void parse_binary(int pipe, struct arg*targ){
	double n;
	int v;
	int current_row = 0;
	int current_col = 0;
	char result[10];
	char results[targ->num_runs+1][*num_lls * 10];
	int current_trajectory = 0;
	int read_col = 0;
	int total_col = 1+(*species-1)*targ->num_runs;
	double timestamp = 0;

	strcpy(result, "");
	for (int i=0; i < targ->num_runs; i++){
		strcpy(results[i], "");
	}
	linked_list *writer;
	writer = targ->arr;
	while (1){
		ssize_t count = read(pipe, &n, sizeof(n));
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
			if (read_col == total_col){
				current_trajectory = 0;
				current_col = 0;
				read_col = 0;
				current_row++;
				//if(current_col != *species-1)
				//	strcat(results[current_trajectory], "\n");
				//printf("\n");
			}
			writer = access_array(targ->arr, current_col, current_row);
			v = (int) n;
			if (current_col == 0){
				pthread_mutex_lock(&targ->mutex_array[current_col+(current_row* *species)]);
				linked_list_add(writer, v);
				pthread_mutex_unlock(&targ->mutex_array[current_col+(current_row* *species)]);
				timestamp = n;
				sprintf(result, "\n%f ", timestamp);
				strcat(results[current_trajectory], result);
				//printf("%s: ", result);
			}else{
				pthread_mutex_lock(&targ->mutex_array[current_col+(current_row* *species)]);
				linked_list_add(writer, v);
				pthread_mutex_unlock(&targ->mutex_array[current_col+(current_row* *species)]);
				sprintf(result, "%i ", v);
				strcat(results[current_trajectory], result);
				//printf("%s ", result);
			}
			//if(targ->p_num == 2)
				//printf("Current col: %i Read col: %i Total col: %i\n", current_col, read_col, total_col);
			if(current_col == *species-1 && current_col != total_col){
				if(current_trajectory < targ->num_runs){
				current_trajectory++;
				sprintf(result, "\n%f ", timestamp);
				strcat(results[current_trajectory], result); 
				current_col = 1;
				}
				//printf("TRAJ: %i, Timestamp: %f\n",current_trajectory, timestamp);
				//fflush(NULL);
			}
			else
				current_col++;
			read_col++;
		}
	}
	//printf("%s\n", &results);
	for (int i = 0; i < targ->num_runs; i++){
		printf("%s", results[i]);
	}
}

void print_status(linked_list *ll, int which){

	/*
	 *TEST SUITE FOR TESTING PURPOSES ONLY
	 */
	char* which_set;
	if (which == 0)
		which_set = "ODDS";
	else
		which_set = "EVENS";

	printf("\n***TEST SUITE***\n");
	linked_list* tptr = ll;
	tptr += TEST_OFFSET;
	node2 *n2ptr = tptr->head2;
	printf("%s @ +%i location, # elements: %i\n", which_set, TEST_OFFSET, tptr->count);
	printf("Elements in linked list Histogram:\n");
	while(n2ptr != NULL){
		printf("[value: %i, count: %i] ,", n2ptr->val, n2ptr->count);
		n2ptr=n2ptr->next;
	}
	printf("\n");
	node *nptr = tptr->head;
	printf("%s @ +%i location, # elements: %i\n", which_set, TEST_OFFSET, tptr->count);
	printf("RAW ELEMENTS:\n");
	while(nptr != NULL){
		printf("[data: %i] ,", nptr->data);
		nptr=nptr->next;
	}
	printf("\n");
}
