#include "semaphore_lib.h"
#include <unistd.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>

int all_locks[1];

void semaphore_init(){

	all_locks[0] = 0;

}

void semaphore_destroy(){
	
	int i;
	for(i=0;i<1;i++){
		
	}
}

int semaphore_lock(int lock_num){

	struct sembuf op;
	op.sem_num = 0;
	op.sem_op = -1;
	op.sem_flg = IPC_NOWAIT;

	return semop(lock_num, &op, 1);

}

int semaphore_unlock(int lock_num){

	struct sembuf op;
	op.sem_num = 0;
	op.sem_op = 1;
	op.sem_flg = 0;

	return semop(lock_num, &op, 1);

}

int create_semaphore(){
	key_t semkey = ftok(".", (1234));

	int lock = semget(semkey, 1, IPC_PRIVATE | IPC_CREAT | 0666);

	union semun {
		int val;
		struct semid_ds *buf;
		ushort *array;
	}argument;
	argument.val = 1;

	semctl( lock, 0, SETVAL, argument );

	all_locks[0] = lock;

	return lock;
}

