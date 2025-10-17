#include "Plots.h"
#include <math.h>
#include <stdio.h>

void plot_polymere(double** R, int N) {
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set title 'Polymere'\n");
        fprintf(gnuplotPipe, "set xlabel 'X'\n");
        fprintf(gnuplotPipe, "set ylabel 'Y'\n");
        fprintf(gnuplotPipe, "plot '-' with linespoints title 'Polymere'\n");

        for (int i = 0; i < N; i++) {
            fprintf(gnuplotPipe, "%f %f\n", R[i][0], R[i][1]);
        }

        fprintf(gnuplotPipe, "e\n");
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    } else {
        printf("Error: Could not open gnuplot pipe.\n");
    }
}

#include "Plots.h"
#include <math.h>
#include <stdio.h>

void Splot_polymere(double** R, int N) {
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set title 'Polymere'\n");
        fprintf(gnuplotPipe, "set xlabel 'X'\n");
        fprintf(gnuplotPipe, "set ylabel 'Y'\n");
        fprintf(gnuplotPipe, "set zlabel 'Z'\n");
        fprintf(gnuplotPipe, "splot '-' with linespoints title 'Polymere'\n");

        for (int i = 0; i < N; i++) {
            fprintf(gnuplotPipe, "%f %f %f\n", R[i][0], R[i][1], R[i][2]);
        }

        fprintf(gnuplotPipe, "e\n");
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    } else {
        printf("Error: Could not open gnuplot pipe.\n");
    }
}

void plotendtoend(double ** stock, int Tp){
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set title 'End to end distance'\n");
        fprintf(gnuplotPipe, "set xlabel 'Time'\n");
        fprintf(gnuplotPipe, "set ylabel 'End to end distance'\n");
        fprintf(gnuplotPipe, "plot '-' with linespoints title 'End to end distance'\n");

        for (int i = 0; i < Tp; i++) {
            fprintf(gnuplotPipe, "%f %f\n", stock[i][1], stock[i][0]);
        }

        fprintf(gnuplotPipe, "e\n");
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    } else {
        printf("Error: Could not open gnuplot pipe.\n");
    }
}



