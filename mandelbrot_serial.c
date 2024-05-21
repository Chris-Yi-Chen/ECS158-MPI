#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <complex.h>

void mandelbrot(int* map, size_t N, double xcenter, double ycenter, double zoom, int cutoff) {
    double complex z_0 = CMPLX(0,0);
    double complex z_n_1, z_n;
    size_t i,j;

    double z = pow(2.0, -zoom);

    // Calculate the half width and half height of the plot area
    double half_width = (N / 2) * z;
    double half_height = (N / 2) * z;

    // Calculate the bounds of the plot area
    double x_min = xcenter - half_width;
    double x_max = xcenter + half_width;
    double y_min = ycenter - half_height;
    double y_max = ycenter + half_height;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            double x = x_min + i * z;
            double y = y_max - j * z; // Note the inversion in the y-axis

            double zx = 0.0, zy = 0.0;
            int iterations = 0;
            while (zx * zx + zy * zy <= 4.0 && iterations < cutoff) {
                double temp = zx * zx - zy * zy + x;
                zy = 2.0 * zx * zy + y;
                zx = temp;
                iterations++;
            }

            // Store the result in the image array
            map[i * N + j] = get_color(iterations, cutoff);
        }
    }
    // while (order < cutoff || cabs(z_n) < 2.0) {
    // }

}
void write_pgm(char *filename, int* map, size_t order, int cutoff) {
    FILE* fp;
    size_t i;
    char* pixels;

    pixels = malloc(order * order);

    for (i = 0; i < order * order; i++) {
        pixels[i] = map[i];
    }

    /* Open file */
	fp = fopen(filename, "wb");
	if (!fp) {
		fprintf(stderr, "Error: cannot open file %s", filename);
		exit(1);
	}
    fprintf(fp, "P5\n%ld %ld\n%ld\n", order, order, cutoff);
    // fwrite(pixels);

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


int main(int argc, char *argv[])
{
	int order, cutoff;
	double xcenter, ycenter, zoom;
    char *filename;
    int* map;

	/* Command line arguments */
	if (argc < 6) {
		fprintf(stderr, "Usage: %s order xcenter ycenter zoom cutoff\n",
				argv[0]);
		exit(1);
	}

	parse_int(argv[1], "order", 128, 8192, &order);
	parse_double(argv[2], "x-coordinate", -10, 10, &xcenter);
	parse_double(argv[3], "y-coordinate",  -10, 10, &ycenter);
	parse_double(argv[4], "zoom", 0, 100, &zoom);
	parse_int(argv[5], "cutoff", 10, 255, &cutoff);

    map = aligned_alloc(64, order * order * sizeof(double complex));

    /* Call implementation */
    mandelbrot(map, order, xcenter, ycenter, zoom, cutoff); 

    /* Save output image */
	filename = malloc(PATH_MAX);
	sprintf(filename, "hmap_%d_%.3lf_%.3lf_%.3lf_%d",
			order, xcenter, ycenter, zoom, cutoff);
	if (argc > 6)
		strcat(filename, argv[6]);
	strcat(filename, ".pgm");
    write_pgm(filename, map, order, cutoff);

    /* Free resources */
	free(filename);
	free(map);

    return 0;
}