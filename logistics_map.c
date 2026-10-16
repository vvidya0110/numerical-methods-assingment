//input to command line (using argc, argv) four values of rmin rmax rstep and Niter
//output 3 columns of numbers: r n x from the logistic map x[n+1] = r x[n]*(1-x[n]
//Starting with a random value for x[0], for each value of r generate Niter values of x by iteration, throwing away the first few (~1000) values of x
//plot the map of x versus r, showing the period doubling bifurcation and the chaotic regime

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {

    if (argc != 5) {
        fprintf(stderr, "Usage: %s r_min r_max step_value nitr\n", argv[0]);
        return 1;
    }


    double r_min = atof(argv[1]);
    double r_max = atof(argv[2]);
    double step_value = atof(argv[3]);
    int nitr = atoi(argv[4]);

    if (step_value <= 0 || r_max <= r_min || nitr <= 1000) {
        fprintf(stderr, "Need step_value > 0, r_max > r_min and nitr > 1000\n");
        return 1;
    }
    FILE *fptr = fopen("data.dat", "w");
    if (fptr == NULL) {
        fprintf(stderr, "ERROR OPENING FILE!\n");
        return 1;
    }


    int nr = (int)((r_max - r_min) / step_value + 0.5);
    for (int k = 0; k <= nr; k++) {
        double r = r_min + k * step_value;
        double p = 0.5;


        for (int t = 1; t <= nitr; t++) {
            p = r * p * (1 - p);


            if (t > 1000) {
                fprintf(fptr, "%f %f\n", r, p);
            }
        }
    }
    fclose(fptr);

    FILE *plot = popen("gnuplot", "w");
    if (plot == NULL) {
        fprintf(stderr, "ERROR: could not start gnuplot!\n");
        return 1;
    }
    fprintf(plot, "set terminal jpeg size 1200,800\n");
    fprintf(plot, "set output 'EE25B152_logistics.jpg'\n");
    fprintf(plot, "set xlabel 'r'\n");
    fprintf(plot, "set ylabel 'x[n]'\n");
    fprintf(plot, "unset key\n");
    fprintf(plot, "plot 'data.dat' using 1:2 with dots\n");
    fprintf(plot, "set output\n");
    pclose(plot);

    return 0;
}
