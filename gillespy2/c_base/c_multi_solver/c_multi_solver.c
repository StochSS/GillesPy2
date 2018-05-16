#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <pthread.h>
#include "linked_list.h"

#define STR_TO_INT_BUF 8
#define READ_BUFFER_SIZE 512
#define BASE_TRAJECTORIES 100

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
int* species;

void* c_solver_runner(void *targ_in);
//void display_results(linked_list *evens, linked_list *odds);
linked_list* access_array(linked_list* ll_ptr, int x, int y, int num_species);

int main (int argc, char** argv){

	int processes, num_runs, num_species, i, num_timesteps;
	processes = atoi(argv[2]);
	num_species = atoi(argv[3])+1;
	species = &num_species;
	num_timesteps = atoi(argv[4]);
	pthread_t *thread_handle = malloc(sizeof(pthread_t)*processes);
	struct arg *targ = malloc(sizeof(struct arg)*processes);
	num_runs = BASE_TRAJECTORIES / processes;
	
	size_t array_size = num_species * num_timesteps;

	linked_list *evens = calloc(array_size, sizeof(linked_list));
	linked_list *odds = calloc(array_size, sizeof(linked_list));
	

	if (argc < 3){
		printf("Insufficient arguments\n");
		exit(1);
	}

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
	}

	for(int i = 0; i < processes; i++){
		pthread_join(thread_handle[i], NULL);
	}

	//display_results(evens, odds);	

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

	for (i = 0; i < targ->num_runs; i++){
		/*
		   Create child processes from threads to exec, redirect stdout from executable to pipe and close it on child end
		 */
		if(pipe(filedes) == -1){
			perror("pipe");
			exit(1);
		}
		if ((pid = fork()) == 0){       //child
			while ((dup2(filedes[1], STDOUT_FILENO) == -1) && (errno == EINTR)){}
			close(filedes[1]);
			close(filedes[0]);
			execv(targ->argv[1], targ->argv);
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
							printf("PTHREAD: %i LL ADDR: %p LL LOC: (%i, %i) LL Count: %i HEAD is: %i TAIL is: %i\n", targ->p_num, writer, current_col, current_row, writer->count, writer->head->data, writer->tail->data);
//fflush(NULL);
							pthread_mutex_unlock(targ->mutex);
							writer++;
							current_col++;
							memset(builder, 0, STR_TO_INT_BUF*sizeof(char));
							current_index = 0;
						}else if (buffer[i] > 0){
	//printf("SEGFAULT HUNTER\n");
	//fflush(NULL);
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

/*
 *Sample Printing of Arrays
 */
void display_results(linked_list *evens, linked_list *odds){

	printf("Evens array:\n");
	for (int i = 0; i < 101; i++){
		printf("row %i[", i);
		for (int j = 0; j < 40; j++){
			//TODO
		}
		printf("]\n");
	}
	printf("\n");
	printf("Odds array:\n");
	for (int i = 0; i < 101; i++){
		printf("row %i[", i);
		for (int j = 0; j < 40; j++){
			//TODO
		}
		printf("]\n");
	}
}

linked_list* access_array(linked_list *ll_ptr, int x, int y, int num_species){
	return ll_ptr+x+(y*num_species);
}
