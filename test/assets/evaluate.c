#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#ifndef EXP1
#define EXP1(x) (0)
#endif
#ifndef EXP2
#define EXP2(x,y) (0)
#endif
#ifndef EXP3
#define EXP3(x,y,z) (0)
#endif
#ifndef EXP4
#define EXP4(x,y,z,t) (0)
#endif

int main(int argc, char **argv)
{
	if (argc < 2) return -1;

	double result = 0.0f;
	double *args = (double*) malloc(sizeof(double) * (argc - 1));
	for (int i = 1; i < argc; ++i)
	{
		args[i-1] = atof(argv[i]);
	}

	switch (argc-1)
	{
	case 1:
		result = EXP1(args[0]);
		break;
	case 2:
		result = EXP2(args[0], args[1]); break;
	case 3:
		result = EXP3(args[0], args[1], args[2]); break;
	case 4:
		result = EXP4(args[0], args[1], args[2], args[3]); break;
	default:
		result = 0.0f;
	}

	free(args);
	printf("%f", result);
	return 0;
}
