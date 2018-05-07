#ifndef semaphore_lib_h
#define semaphore_lib_h

extern int all_locks[1];

int semaphore_lock(int lock_num);

int semaphore_unlock(int lock_num);

int create_semaphore();

void semaphore_init();

void semaphore_destroy();



#endif // semaphore_lib_h
