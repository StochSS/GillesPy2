#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include "semaphore_lib.h"

#define STR_TO_INT_BUF 8
#define READ_BUFFER_SIZE 512
#define ROWS 1024
#define COLUMNS 64
#define MAX_ARRAY_SIZE sizeof(int)*ROWS*COLUMNS
#define ACCESS_ARRAY(target, x, y) (target+x+(y*COLUMNS))

int main (int argc, char** argv){

	int processes, p_num, status, total_num_runs, num_runs, i, evens_shmid, odds_shmid, evens_sem, odds_sem, sem;
	int *evens, *odds, *ptr;
	int filedes[2];
	pid_t pid, parent;
	processes = atoi(argv[2]);
	pid_t mid_pids[processes];
	p_num = 0;
	total_num_runs = 5;
	evens = NULL;
	odds = NULL;
	num_runs = total_num_runs / processes;
	parent = getpid();
	key_t evens_key = ftok(".", 'x');
	key_t odds_key = ftok(".", 'y');

	if (argc < 3){
		printf("Insufficient arguments\n");
		exit(1);
	}

	/*
	   Initialize and attach shared memory for evens and odds
	   */
	if ((evens_shmid = shmget(evens_key, MAX_ARRAY_SIZE, IPC_CREAT | 0666)) < 0){
		perror("***FAILED TO ALLOCATE SHARED MEMORY FOR EVENS***\n");
		exit(1);
	}

	if ((odds_shmid = shmget(odds_key, MAX_ARRAY_SIZE, IPC_CREAT | 0666)) < 0){
		perror("***FAILED TO ALLOCATE SHARED MEMORY FOR ODDS***\n");
		exit(1);
	}
	if (( evens = (int*)shmat(evens_shmid, NULL, 0)) == (void*) -1){
		perror("***FAILED to attach to shared memory***\n");
		exit(1);
	}
	if (( odds = (int*)shmat(odds_shmid, NULL, 0)) == (void*) -1){
		perror("***FAILED to attach to shared memory***\n");
		exit(1);
	}

	evens_sem = create_semaphore();
	odds_sem = create_semaphore();

	/*
	   Create middle processes and assign p_num
	   */
	for(i = 0; i < processes; i++){
		if (parent == getpid()){
			p_num = i+1;
			mid_pids[i] = fork();
		}
	}

	if (parent != getpid()){ //middle processes
		/*
		 * Assign generic pointers and semaphores for middle processes
		 */

		if (p_num % 2 == 0){
			ptr = evens;
			sem = evens_sem;
		}else{
			ptr = odds;
			sem = odds_sem;
		}



		/*
		   If total number of runs is not divisible by the number of processes, add the remainder runs to the first processes
		   */

		if ((total_num_runs % processes) == p_num)
			num_runs++;

		for (i = 0; i < num_runs; i++){
			/*
			   Create child processes from middle processes to exec, redirect stdout from executable to pipe and close it on child end
			   */
			/*
			   Create a pipe in each middle process
			   */
			if(pipe(filedes) == -1){
				perror("pipe");
				exit(1);
			}
			if ((pid = fork()) == 0){	//child
				while ((dup2(filedes[1], STDOUT_FILENO) == -1) && (errno == EINTR)){}
				close(filedes[1]);
				close(filedes[0]);
				execv(argv[1], argv);
				perror("execv");

				_exit(1);
			}else{
				/*
				   Close Pipe output in middle processes and read pipe data to buffer
				   --*/
				close(filedes[1]);
				char buffer[READ_BUFFER_SIZE];
				int current_index = 0;
				int current_row = 0;
				int current_col = 0;
				char *reader;
				int *writer;
				writer = ptr;
				char builder[STR_TO_INT_BUF];
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
								writer = ACCESS_ARRAY(ptr, current_col, current_row);
							}else if(*reader == ' ') {
								semaphore_lock(sem);
								*writer += atoi(builder);
								semaphore_unlock(sem);
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
		exit(0);
	}//end of Middle Processes

	//wait for middle processes to die
	for(int i = 0; i < processes; i++){
		waitpid(mid_pids[i], &status, 0);
	}
	printf("Evens array: ");
	for (int i = 0; i < 10; i++){
		printf("row %i[", i);
		for (int j = 0; j < 40; j++){
			printf("%i, ", *ACCESS_ARRAY(evens, j, i));
		}
		printf("]\n");
	}
	printf("\n");
	printf("Odds array: ");
	for (int i = 0; i < 10; i++){
		printf("row %i[", i);
		for (int j = 0; j < 40; j++){
			printf("%i, ", *ACCESS_ARRAY(odds, j, i));
		}
		printf("]\n");
	}
	if ((shmdt(evens)) < 0){
		perror("SHMDT ERROR EVENS");
	}
	if ((shmdt(odds)) < 0){
		perror("SHMDT ERROR ODDS");
	}
	shmctl(evens_shmid, IPC_RMID, NULL);
	shmctl(odds_shmid, IPC_RMID, NULL);
	return 0;
}
