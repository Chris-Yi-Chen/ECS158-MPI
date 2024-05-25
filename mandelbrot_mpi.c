#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>

// mpicc -Wall -Wextra -Werror -O2 -o mandelbrot_mpi mandelbrot_mpi.c -lm
// mpirun -nolocal --hostfile csif_hostfile --host pc27,pc25 -n 16 ./mandelbrot_mpi 1024 -0.722 0.246 15 255

#define MANAGER_NODE 0

/* Struct for input args */
struct mandelbrot_args {
	int order;
    int cutoff;
	double xcenter;
    double ycenter;
    double zoom;
};

MPI_Datatype mpi_mandelbrot_args;

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

void mandelbrot(int* map, size_t N, double xcenter, double ycenter, double zoom, int cutoff, int start, int end) {
    size_t i,j;

    // z: distance between points
    double dist = pow(2.0, -zoom);

    // Calculate the half width and half height of the plot area
    double half_width = (N / 2) * dist;
    double half_height = (N / 2) * dist;

    // Calculate the bounds of the plot area
    double x_min = xcenter - half_width;
    double y_max = ycenter + half_height;
    
    // iterate through x points 
    for (i = 0; i < N; i++) {
        // iterate through y points
        for (j = start; j < (size_t)end; j++) {
            // x and y coordinates
            double x_coords = x_min + (i * dist);
            double y_coords = y_max - (j * dist); 

            double zx = 0.0, zy = 0.0;
            
            int iterations = 0;
            // cx * cx + cy * cy <= 4.0 AKA |c| <= 2
            while (zx * zx + zy * zy <= 4.0 && iterations < cutoff) {
                // temp because zy uses the old zx value
                // z^2 = (a^2 - b^2) + (2ab)i
                double temp_zx = (zx * zx - zy * zy) + x_coords;
                zy = 2.0 * zx * zy + y_coords;
                zx = temp_zx;
                iterations++;
            }
            map[(j-start) * N + i] = (int)(iterations);
        }
    }
}

void worker_loop(void) {
    struct mandelbrot_args margs;
    /* Grab CLI from broadcast */
	MPI_Bcast(&margs, sizeof(struct mandelbrot_args), MPI_BYTE, MANAGER_NODE, MPI_COMM_WORLD);

    int rank, num_proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    fprintf(stderr, "Rank: %d\n", rank);
    fprintf(stderr, "Number of processes: %d\n\n", num_proc);

    int rows_per_worker = margs.order / (num_proc - 1);
    int start = (rank - 1) * rows_per_worker;

    // if last worker, finish the remaining (to N)

    int end = (rank == num_proc - 1) ? margs.order : start + rows_per_worker;

    int num_rows = end - start;

    fprintf(stderr, "Rows per worker: %d\n", rows_per_worker);
    fprintf(stderr, "Start row: %d\n", start);
    fprintf(stderr, "End row: %d\n", end);
    fprintf(stderr, "Number of rows: %d\n", num_rows);

    int* work_map = malloc(margs.order * num_rows * sizeof(int));
    mandelbrot(work_map, margs.order, margs.xcenter, margs.ycenter, margs.zoom, margs.cutoff, start, end);
    fprintf(stderr, "before gather\n");

    MPI_Gatherv(work_map, margs.order * num_rows, MPI_INT, NULL, NULL, NULL, MPI_INT, MANAGER_NODE, MPI_COMM_WORLD);
    fprintf(stderr, "after gather\n");

    free(work_map);

}



void manager_loop(int num_workers, struct mandelbrot_args *margs, int* map) {
	MPI_Bcast(margs, sizeof(struct mandelbrot_args), MPI_BYTE, MANAGER_NODE, MPI_COMM_WORLD);

    int rows_per_worker = margs->order / num_workers;
    int remaining_rows = margs->order % num_workers;

    fprintf(stderr, "rows_per_worker: %d\n", rows_per_worker);
    fprintf(stderr, "remaining_rows: %d\n", remaining_rows);
    fprintf(stderr, "num_workers: %d\n", num_workers);

    int* recvcounts = malloc(num_workers * sizeof(int));
    int* displs = malloc(num_workers * sizeof(int));

    

    for (int i = 0; i < num_workers; i++) {
        // # of rows for each process, last proc may have a different value
        recvcounts[i] = (i == num_workers - 1) ? (rows_per_worker + remaining_rows) * margs->order : rows_per_worker * margs->order;
        displs[i] = i * rows_per_worker * margs->order;
        fprintf(stderr, "i: %d\nrecvcounts: %d\n", i, recvcounts[i]);
        fprintf(stderr, "displs: %d\n\n", displs[i]);
    }


    MPI_Gatherv(NULL, 0, MPI_INT, map, recvcounts, displs, MPI_INT, MANAGER_NODE, MPI_COMM_WORLD);
    fprintf(stderr, "After manager gathers\n");

    free(recvcounts);
    free(displs);
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
    int num_proc;
    int* map;

	/* Command line arguments */
	if (argc < 6) {
		fprintf(stderr, "Usage: %s order xcenter ycenter zoom cutoff\n",
				argv[0]);
		exit(1);
	}

    /* Get # of tasks total */
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    fprintf(stderr, "from manager: num_proc: %d\n", num_proc);
    if (num_proc == 1) {
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
    map = aligned_alloc(64, margs.order * margs.order * sizeof(int));


    /* Initialize custom MPI args */
    int blocklen[] = {1, 1, 1, 1, 1};
    MPI_Aint disp[] = {offsetof(struct mandelbrot_args, order), offsetof(struct mandelbrot_args, cutoff),
                       offsetof(struct mandelbrot_args, xcenter), offsetof(struct mandelbrot_args, ycenter),
                       offsetof(struct mandelbrot_args, zoom)};
    MPI_Datatype type[] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Type_create_struct(5, blocklen, disp, type, &mpi_mandelbrot_args);
    MPI_Type_commit(&mpi_mandelbrot_args);


    /* Call manager loop */
    manager_loop(num_proc - 1, &margs, map);


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
    free(mpi_mandelbrot_args);

}


int main(int argc, char *argv[])
{
	int rank;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == MANAGER_NODE) {
		manager_main(argc, argv);
    }
	else {
		worker_loop();
    }

	MPI_Finalize();
	return 0;
}