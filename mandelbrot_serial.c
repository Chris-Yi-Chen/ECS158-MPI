#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

//gcc -g -Wall -Wextra -Werror -O2 .\mandelbrot_serial.c -o mandelbrot_serial -lm
void mandelbrot(int* map, size_t N, double xcenter, double ycenter, double zoom, int cutoff) {
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
        for (j = 0; j < N; j++) {
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
            int greyscaleVal = (int)iterations;
            // Store the result in the image array
            map[j * N + i] = greyscaleVal;
        }
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

    map = aligned_alloc(64, order * order * sizeof(int));

    /* Call implementation */
    mandelbrot(map, order, xcenter, ycenter, zoom, cutoff); 

    /* Save output image */
	filename = malloc(PATH_MAX);
	sprintf(filename, "mandel_%d_%.3lf_%.3lf_%.3lf_%d",
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