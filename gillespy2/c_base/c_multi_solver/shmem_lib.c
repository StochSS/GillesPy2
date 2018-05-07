#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/ipc.h>
#include <sys/types.h>
#include <sys/shm.h>
#include "semaphore_lib.h"

#define ROWS 100000
#define COLUMNS 50
#define ACCESS_ARRAY(target, x, y) *(&target+x*sizeof(int)+(y*sizeof(int)*ROWS))

int shmid, sem;
key_t shmkey;
int* ptr;
int row = 0;
int col = 0;

void shmem_init(int *target){
	
	sem  = create_semaphore();
	shmkey = ftok(".", 'x');
	ptr = target;

	if ((shmid = shmget(shmkey, sizeof(int[ROWS][COLUMNS]), IPC_CREAT | 0666)) < 0){
		perror("*** failed to create shared memory!***");
		exit(1);
	}else
		printf("SUCCESSFUL INIT shmid = %i\n", shmid);

}

void attach_to_shmem(){
	printf("[%i] attempting to attach to shmid %i at %p...\n", getpid(), shmid, ptr);
	if (( ptr = (int*)shmat(shmid, NULL, 0)) == (void*) -1){
		perror("***Failed to attach to shmem***");
		exit(1);
	}else{
		printf("[%i] SUCCESSFUL ATTACH to shmid %i\n",getpid(), shmid);
	}
}

void detach_from_shmem(){
	printf("[%i] attempting to detatch shmid %i from %p...\n", getpid(), shmid, ptr);
	if (shmdt(ptr) == -1){
		perror("***Failed to detatch from shmem***");
	}else{
		//printf("[%i] SUCCESSFUL DETACH from shmid %i\n",getpid(), shmid);
	}
}

void add_value(int value){
	semaphore_lock(sem);
	//printf("value %i to be added\n", value);
	*ptr = value;
	ptr += sizeof(int);
	semaphore_unlock(sem);

	//printf("Adding value %i to location %p\n", value, ptr);
	
}

int destroy_shmem(){
	shmctl(shmid, IPC_RMID, NULL);
	return *ptr;
}
