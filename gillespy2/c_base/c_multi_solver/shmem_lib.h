#ifndef shmem_lib_h
#define shmem_lib_h

void shmem_init(int *target);
void attach_to_shmem();
void add_value(int value);
void detach_from_shmem();
int destroy_shmem();


#endif
