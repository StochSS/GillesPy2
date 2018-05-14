#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <pthread.h>

#define STR_TO_INT_BUF 8
#define READ_BUFFER_SIZE 512
#define ROWS 1024
#define COLUMNS 64
#define MAX_ARRAY_SIZE sizeof(int)*ROWS*COLUMNS
#define ACCESS_ARRAY(target, x, y) (target+x+(y*COLUMNS))

pthread_mutex_t evens_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t odds_mutex = PTHREAD_MUTEX_INITIALIZER;

struct arg {
	int p_num;
	int num_runs;
	int* arr;
	char** argv;
	pthread_mutex_t mutex;
};

void* c_solver_runner(void *targ_in);

int main (int argc, char** argv){

	int processes, total_num_runs, i;
	int *evens = malloc(MAX_ARRAY_SIZE);
	int *odds = malloc(MAX_ARRAY_SIZE);
	processes = atoi(argv[2]);
	pthread_t thread_handle[processes];
	struct arg targ[processes];
	total_num_runs = 100;
	int num_runs = total_num_runs / processes;

	if (argc < 3){
		printf("Insufficient arguments\n");
		exit(1);
	}

	/*
	   Create threads and assign values
	 */
	for(i = 0; i < processes; i++){
		targ[i].p_num = i+1;
		if ( targ[i].p_num % 2 == 0){
			targ[i].arr = evens;
			targ[i].mutex = evens_mutex;
		}
		else{
			targ[i].arr = odds;
			targ[i].mutex = odds_mutex;
		}
		targ[i].argv = argv;
		targ[i].num_runs = num_runs;
		if (total_num_runs % processes >= targ[i].p_num){
			targ[i].num_runs++;
		}
		pthread_create(&thread_handle[i], NULL, c_solver_runner, &targ[i]);
	}

	for(int i = 0; i < processes; i++){
		pthread_join(thread_handle[i], NULL);
	}

	printf("Evens array:\n");
	for (int i = 0; i < 10; i++){
		printf("row %i[", i);
		for (int j = 0; j < 40; j++){
			printf("%i, ", *ACCESS_ARRAY(evens, j, i));
		}
		printf("]\n");
	}
	printf("\n");
	printf("Odds array:\n");
	for (int i = 0; i < 10; i++){
		printf("row %i[", i);
		for (int j = 0; j < 40; j++){
			printf("%i, ", *ACCESS_ARRAY(odds, j, i));
		}
		printf("]\n");
	}

	return 0;
}
/*
 * Thread Function
 */
void *c_solver_runner(void *targ_in){

	int i;
	int *ptr;
	int filedes[2];
	pid_t pid;

	struct arg* targ =  (struct arg*) targ_in;

	char* filename = targ->argv[1];

	printf("Thread %i running with %i runs\n", targ->p_num, targ->num_runs);

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
			execv(filename, targ->argv);
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
			char *reader;
			int *writer;
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
					reader = &buffer[0];
					while(reader <= &buffer[512]){
						if (*reader == '\n'){
							current_index = 0;
							current_row++;
							current_col = 0;
							writer = ACCESS_ARRAY(targ->arr, current_col, current_row);
						}else if(*reader == ' ') {
							pthread_mutex_lock(&targ->mutex);
							*writer += atoi(builder);
							pthread_mutex_unlock(&targ->mutex);
							writer++;
							current_col++;
							memset(builder, 0, STR_TO_INT_BUF*sizeof(char));
							current_index = 0;
						}else if (*reader > 0){
							builder[current_index] = *reader;
							current_index++;
						}
						reader++;
					}

				}
			}
		}
		close(filedes[0]);
	}
	pthread_exit(0);
}
