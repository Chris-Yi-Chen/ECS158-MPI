#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>


#define MANAGER_NODE 0

/* Global variable */
int* map;

/* Struct for input args */
struct mandelbrot_args {
	int order;
    int cutoff;
	double xcenter;
    double ycenter;
    double zoom;
};

void worker_loop(void) {
    struct mandelbrot_args margs;
    /* Grab CLI from broadcast */
	MPI_Bcast(&margs, sizeof(struct mandelbrot_args), MPI_BYTE, MANAGER_NODE, MPI_COMM_WORLD);

    while (1) {

        /* Recieve chunk to calculate */

        /* Stop order: break loop and quit DO WE NEED THIS? */



    }
}

void write_pgm(char *filename, int* map, size_t N, int cutoff) {
    FILE* fp;
    size_t i;
    char* pixels;

    pixels = malloc(N * N);

    for (i = 0; i < N * N; i++) {
        pixels[i] = map[i];
    }

    /* Open file */
	fp = fopen(filename, "wb");
	if (!fp) {
		fprintf(stderr, "Error: cannot open file %s", filename);
		exit(1);
	}
    fprintf(fp, "P5\n%ld %ld\n%d\n", N, N, cutoff);
    fwrite(pixels, sizeof(char), N * N, fp);

    free(pixels);
    fclose(fp);

}

void manager_loop(int nworkers, struct mandelbrot_args *margs) {

	MPI_Bcast(margs, sizeof(struct mandelbrot_args), MPI_BYTE, MANAGER_NODE, MPI_COMM_WORLD);


}
/*
 * CLI argument parsing
 */
void parse_int(char *str, char *val, int min, int max, int *num)
{
	int n = atoi(str);
	if (n < min || n > max) {
		fprintf(stderr, "Error: wrong %s (%d <= N <= %d)", val, min, max);
		exit(1);
	}
	*num = n;
}

void parse_double(char *str, char *val, double min, double max, double *num)


{
	double n = atof(str);
	if (n < min || n > max) {
		fprintf(stderr, "Error: wrong %s (%lf <= N <= %lf)", val, min, max);
		exit(1);
	}
	*num = n;
}

void manager_main(int argc, char *argv[])
{
    struct mandelbrot_args margs;
    char *filename;
    int ntasks;


	/* Command line arguments */
	if (argc < 6) {
		fprintf(stderr, "Usage: %s order xcenter ycenter zoom cutoff\n",
				argv[0]);
		exit(1);
	}

    /* Get # of tasks total */
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    if (ntasks == 1) {
		fprintf(stderr, "Error: not enough tasks\n");
        exit(1);
    }

    /* Parse the arguments */
	parse_int(argv[1], "order", 128, 8192, &margs.order);
	parse_double(argv[2], "x-coordinate", -10, 10, &margs.xcenter);
	parse_double(argv[3], "y-coordinate",  -10, 10, &margs.ycenter);
	parse_double(argv[4], "zoom", 0, 100, &margs.zoom);
	parse_int(argv[5], "cutoff", 10, 255, &margs.cutoff);

    /* Allocate Global Variable */
    map = aligned_alloc(64, order * order * sizeof(int));

    /* Call manager loop */
    manager_loop(ntasks - 1, &margs);


    /* Save output image */
	filename = malloc(PATH_MAX);
	sprintf(filename, "mandel_%d_%.3lf_%.3lf_%.3lf_%d",
			margs.order, margs.xcenter, margs.ycenter, margs.zoom, margs.cutoff);
	if (argc > 6)
		strcat(filename, argv[6]);
	strcat(filename, ".pgm");
    write_pgm(filename, map, margs.order, margs.cutoff);

    /* Free resources */
	free(filename);

}


int main(int argc, char *argv[])
{
	int rank;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == MANAGER_NODE)
		manager_main(argc, argv);
	else
		worker_loop();

	MPI_Finalize();
	return 0;
}