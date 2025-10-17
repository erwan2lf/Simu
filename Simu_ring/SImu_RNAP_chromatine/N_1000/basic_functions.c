#include "basic_functions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>
#include "potentiels.h"
// #include </Users/erwan/Documents/These/MTwister/mt19937ar.h>
#include </home/elefloch/Simulation/MT/mt19937ar.c>
#include <string.h>
#include <sys/stat.h> // Pour mkdir()
#include <errno.h> 
#include <omp.h>



#define FILENAME_SIZE 100

#define PI 3.14159265358979323846


double randn() {
    //unsigned long seed=97714454678; init_genrand(seed);
    static int hasSpare = 0;
    static double spare;
    if (hasSpare) {
        hasSpare = 0;
        return spare;
    }
    hasSpare = 1;
    double u, v, s;
    do {
        u = genrand_real2() * 2 - 1;
        v = genrand_real2() * 2 - 1;
        //u = (rand() / ((double) RAND_MAX)) * 2.0 - 1;
        //v = (rand() / ((double) RAND_MAX)) * 2.0 - 1;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);
    s = sqrt(-2.0 * log(s) / s);
    spare = v * s;
    return u * s;
}

// Fonction pour calculer la distance entre deux points en 3D
double distance(double* p1, double* p2) {
    return sqrt(pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2) + pow(p1[2] - p2[2], 2));
}

// Fonction pour calculer la distance au carré entre deux points en 3D
double distance2(double* p1, double* p2) {
    return pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2) + pow(p1[2] - p2[2], 2);
}

double norm(double* vec) {
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}
double norm2(double* vec) {
    return vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
}

// Fonction pour afficher la structure
void afficher_structure(double** R_matrix, int N) {
    //printf("Nombre de particules : %d\n", N);
    for (int i = 0; i < N; i++) {
        printf("Particule %d : x = %lf, y = %lf, z = %lf\n", i+1, R_matrix[i][0], R_matrix[i][1], R_matrix[i][2]);
    }
}

double** allocate_matrix(size_t rows, size_t cols) {
    double** matrix = (double**)malloc(rows * sizeof(double*));
    if (matrix == NULL) {
        fprintf(stderr, "Memory allocation failed for matrix\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double*)malloc(cols * sizeof(double));
        if (matrix[i] == NULL) {
            fprintf(stderr, "Memory allocation failed for matrix row\n");
            exit(EXIT_FAILURE);
        }
    }
    return matrix;
}

bool is_allocated(void *ptr) {
    return ptr != NULL;
}

void free_if_allocated(void **ptr) {
    if (is_allocated(*ptr)) {
        free(*ptr);
        *ptr = NULL;
    }
}

void free_matrix_if_allocated(double ***matrix, int rows) {
    if (is_allocated(*matrix)) {
        for (int i = 0; i < rows; i++) {
            if (is_allocated((*matrix)[i])) {
                free((*matrix)[i]);
                (*matrix)[i] = NULL; // Prevent double free
            }
        }
        free(*matrix);
        *matrix = NULL; // Prevent double free
    }
}

void free_matrix_cube_if_allocated(double ***matrix, int rows, int cols) {
    if (is_allocated(matrix)) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (is_allocated(matrix[i][j])) {
                    free(matrix[i][j]);
                    matrix[i][j] = NULL; // Prevent double free
                }
            }
            free(matrix[i]);
            matrix[i] = NULL; // Prevent double free
        }
        free(matrix);
        matrix = NULL; // Prevent double free
    }
}


void free_matrix_4D(double ****matrix, int Nm, int T_msd) {
    if (matrix != NULL) {
        for (int i = 0; i < Nm; i++) {
            if (matrix[i] != NULL) {
                for (int j = 0; j < T_msd; j++) {
                    if (matrix[i][j] != NULL) {
                        for (int l = 0; l < 3; l++) {
                            if (matrix[i][j][l] != NULL) {
                                free(matrix[i][j][l]);
                                matrix[i][j][l] = NULL; // Prevent double free
                            }
                        }
                        free(matrix[i][j]);
                        matrix[i][j] = NULL; // Prevent double free
                    }
                }
                free(matrix[i]);
                matrix[i] = NULL; // Prevent double free
            }
        }
        free(matrix);
        matrix = NULL; // Prevent double free
    }
}



int* particules_outside(double** R, double* origine, double rayon, double epaisseur, int N) {
    int* outside = (int*)malloc(N * sizeof(int));
    for (int part = 0; part < N; ++part) {
        double distance_origine = sqrt(pow(R[part][0] - origine[0], 2) + pow(R[part][1] - origine[1], 2));
        double d_z = fabs(R[part][2] - origine[2]);

        if (distance_origine > rayon && d_z < epaisseur) {
            outside[part] = 1;
        } else {
            outside[part] = 0;
        }
    }
    return outside;
}

//typedef struct {double moyenne;double std;} Mesures;
Mesures calcul_mesures (double ** R, int N){
    double liste_distance[N-1]; double moyenne = 0.0 ; double std = 0.0;
    for (int part=0; part<N-1; part ++){
    double d = distance(R[part], R[part+1]);
    liste_distance[part]= d;
    moyenne += d; }
    moyenne = moyenne / (N-1) ;
    for (int part=0; part<N-1; part ++){std += (moyenne - liste_distance[part]) * (moyenne - liste_distance[part]);}
    std = sqrt(std/(N-1));

    //double mesures[2] ; mesures [0] =moyenne ; mesures[1] = std ;
    //return moyenne, std ;return mesures ;
    Mesures result;
    result.moyenne = moyenne;
    result.std = std;
    //printf("%f %f \n", moyenne, result.moyenne);
    return result;
    }

/*Mesures_temp calcul_mesures_temp (double ** R, int N,  int t){
    
}*/

void enregistrement(FILE* fichier, double** R_matrix, int N, int t) {
    double TT = 1e+3;
    double cdm[3] = {0, 0, 0};
    fprintf(fichier, "ITEM: TIMESTEP\n");
    fprintf(fichier, "%d\n", t);
    fprintf(fichier, "ITEM: NUMBER OF ATOMS\n");
    fprintf(fichier, "%d\n", N);
    fprintf(fichier, "ITEM: BOX BOUNDS ss ss ss\n");
    fprintf(fichier, "%lf %lf\n", -TT, TT);
    fprintf(fichier, "%lf %lf\n", -TT, TT);
    fprintf(fichier, "%lf %lf\n", -TT, TT);
    fprintf(fichier, "ITEM: ATOMS id type xs ys zs\n");

    for(int i = 0; i < N; i++){
        cdm[0] += R_matrix[i][0];
        cdm[1] += R_matrix[i][1];
        cdm[2] += R_matrix[i][2];
    }
    cdm[0] /= N;
    cdm[1] /= N;
    cdm[2] /= N;

    for (int particle = 0; particle < N; particle++) {
        fprintf(fichier, "%d 1 %lf %lf %lf -2\n", particle, (R_matrix[particle][0]) / 1000, (R_matrix[particle][1]) / 1000, (R_matrix[particle][2]) / 1000);
    }
}

double dot_product(double *vec1, double *vec2) {
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

double calculate_mean_for_correlation(double **R, int i, int j) {
    double sum = 0.0;
    for (int k = i; i < i + j; k++) {
        sum += norm(R[k]);
    }
    return sum / (i+j);

}

void plot_columns_from_file(const char *filename, const char *title, const char *xlabel, const char *ylabel) {
    // Open the file for reading
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error: Could not open file %s for reading.\n", filename);
        return;
    }

    // Create a temporary file to store the data for Gnuplot
    FILE *tempFile = fopen("temp_data.txt", "w");
    if (tempFile == NULL) {
        fprintf(stderr, "Error: Could not create temporary file.\n");
        fclose(file);
        return;
    }

    // Read data from the input file and write to the temporary file
    char line[1024];
    while (fgets(line, sizeof(line), file)) {
        fprintf(tempFile, "%s", line);
    }

    fclose(file);
    fclose(tempFile);

    // Open a pipe to Gnuplot
    FILE *gnuplotPipe = popen("gnuplot -persist", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set title '%s'\n", title);
        fprintf(gnuplotPipe, "set xlabel '%s'\n", xlabel);
        fprintf(gnuplotPipe, "set ylabel '%s'\n", ylabel);
        fprintf(gnuplotPipe, "plot for [col=2:*] 'temp_data.txt' using 1:col with linespoints title columnheader\n");
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    } else {
        fprintf(stderr, "Error: Could not open gnuplot pipe.\n");
    }
    // Supprime le fichier temporaire@
    if (remove("temp_data.txt") != 0) {
        fprintf(stderr, "Error: Could not delete temporary file.\n");
    }
}

void calculate_autocorrelation(double** Rbb, double *stock_correl, int T_correlation, int N) {

    double mean_rbb[3] = {0.0, 0.0, 0.0};

    
    for(int i = 0; i < T_correlation; i++)
    {
        mean_rbb[0] += Rbb[i][0];
        mean_rbb[1] += Rbb[i][1];
        mean_rbb[2] += Rbb[i][2];
    }
    mean_rbb[0] = mean_rbb[0] / T_correlation;
    mean_rbb[1] = mean_rbb[1] / T_correlation;
    mean_rbb[2] = mean_rbb[2] / T_correlation;

    for (int i = 0; i < T_correlation; i++) 
    {
        double correlation = 0;

        for (int j = 0; j < T_correlation - i; j++)
        {
            double adjusted_vec1[3] = {0, 0, 0};
            double adjusted_vec2[3] = {0, 0, 0};
            for (int k = 0; k < 3; k++)
            {
                //adjusted_vec1[k] = Rbb[i + j][k] - mean_rbb[k];
                //adjusted_vec2[k] = Rbb[j][k] - mean_rbb[k];
                adjusted_vec1[k] = Rbb[i + j][k];
                adjusted_vec2[k] = Rbb[j][k];
            }
            correlation += dot_product(adjusted_vec1, adjusted_vec2);

        }
        stock_correl[i] += correlation / (T_correlation - i);
    }
}


void plot_autocorrelation(double *time, double *autocorrelation, int length, const char *filename, const char *title, const char *xlabel, const char *ylabel) {
    // Create a file to store the data for Gnuplot
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error: Could not open file %s for writing\n", filename);
        return;
    }

    // Write the time and autocorrelation data to the file
    for (int i = 0; i < length; i++) {
        fprintf(file, "%f %f\n", time[i], autocorrelation[i]);
    }
    fclose(file);

    // Plot data using Gnuplot
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set title '%s'\n", title);
        fprintf(gnuplotPipe, "set xlabel '%s'\n", xlabel);
        fprintf(gnuplotPipe, "set ylabel '%s'\n", ylabel);
        fprintf(gnuplotPipe, "plot '%s' using 1:2 with linespoints title 'Autocorrelation'\n", filename);
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    } else {
        fprintf(stderr, "Error: Could not open gnuplot pipe\n");
    }
}

void plot_msd(double **stock_msd, int Tp, int Nm, int *list_monomere, const char *base_filename, const char *base_log_filename){
    /////////////////////////////// Plot ////////////////////////////
    for (int j = 0; j < Nm; j++) {
        char msd_filename[FILENAME_SIZE];
        snprintf(msd_filename, FILENAME_SIZE, "%s_%d.txt", base_filename, list_monomere[j]);

        FILE *msd_file = fopen(msd_filename, "w");
        if (msd_file == NULL) {
            fprintf(stderr, "Error: Could not open file %s for writing\n", msd_filename);
            return;
        }

        for (int i = 0; i < Tp; i++) {
            fprintf(msd_file, " %d %f \n", i, stock_msd[i][j]);
        }
        fclose(msd_file);

        char title[50];
        snprintf(title, sizeof(title), "MSD du monomere %d", list_monomere[j]);
        const char *xlabel = " Temps";
        const char *ylabel = " MSD ";

        //plot_columns_from_file(msd_filename, title, xlabel, ylabel);
    }
    ////////////////////// log log plot ///////////////////////////////
     for (int j = 0; j < Nm; j++) {
        char msd_log_filename[FILENAME_SIZE];
        snprintf(msd_log_filename, FILENAME_SIZE, "%s_%d.txt", base_log_filename, list_monomere[j]);

        FILE *msd_log_file = fopen(msd_log_filename, "w");
        if (msd_log_file == NULL) {
            fprintf(stderr, "Error: Could not open file %s for writing\n", msd_log_filename);
            return;
        }

        for (int i = 0; i < Tp; i++) {
            fprintf(msd_log_file, " %f %f \n", log(i), log(stock_msd[i][j]));
        }
        fclose(msd_log_file);

        char title[50];
        snprintf(title, sizeof(title), "MSD du monomere %d", list_monomere[j]);
        const char *xlabel = " Temps";
        const char *ylabel = " MSD ";

        //plot_columns_from_file(msd_log_filename, title, xlabel, ylabel);
    }
    
}

void save_centre_de_masse(double** R_centre_de_masse, int N, int T_centre_de_masse, int nbr_simu, int periode_centre_de_masse){
    const char* folder_name = "CENTRE_DE_MASSE";
    const char* folder_name_2 = "VITESSE_CENTRE_DE_MASSE";
    if(nbr_simu == 0){ 
        if (mkdir(folder_name, 0777) == -1){
            if (errno == EEXIST){
                printf("Le dossier '%s' existe déjà \n", folder_name);
            }
            else {
                perror("Erreur lors de la création du dossier \n");
            }
        }
        else {
            printf("Le dossier '%s' a été crée avec succès \n", folder_name);
        }

        if (mkdir(folder_name_2, 0777) == -1){
            if (errno == EEXIST){
                printf("Le dossier '%s' existe déjà \n", folder_name_2);
            }
            else {
                perror("Erreur lors de la création du dossier \n");
            }
        }
        else {
            printf("Le dossier '%s' a été crée avec succès \n", folder_name_2);
        }
    }

    char file_path[512]; 
    snprintf(file_path, sizeof(file_path), "%s/centre_de_masse_%d.txt", folder_name, nbr_simu);

    FILE *file = fopen(file_path, "w");
    if(file == NULL){
        perror("Erreur lors de l'ouverture du fichier \n");
    }

    for(int i = 0; i < T_centre_de_masse; i++){
        double t = i * periode_centre_de_masse; 
        fprintf(file, "%f %f %f %f \n", t, R_centre_de_masse[i][0], R_centre_de_masse[i][1], R_centre_de_masse[i][2]);
    }

    fclose(file);


    char file_path_vitesse[512]; 
    snprintf(file_path_vitesse, sizeof(file_path_vitesse), "%s/vitesse_centre_de_masse_%d.txt", folder_name_2 ,nbr_simu);

    FILE *file_vitesse = fopen(file_path_vitesse, "w");
    if(file_vitesse == NULL){
        perror("Erreur lors de l'ouverture du fichier vitesse\n");
        return; 
    }

    for (int i = 0; i < T_centre_de_masse-1; i++){
        double t = i * periode_centre_de_masse; 
        double v_x = R_centre_de_masse[i+1][0] - R_centre_de_masse[i][0]; 
        double v_y = R_centre_de_masse[i+1][1] - R_centre_de_masse[i][1]; 
        double v_z = R_centre_de_masse[i+1][2] - R_centre_de_masse[i][2]; 

        fprintf(file_vitesse, "%f %f %f %f \n", t, v_x, v_y, v_z);
    }

    fclose(file_vitesse);


}

void calculate_msd(double ****R_monomere_arrays,
                   double ***stock_msd,
                   int Tp, int Nm,
                   int *list_monomere,
                   int nbr_simu, int nbr_total_simu)
{
    const char *folder_name = "MSD";

    // Création du dossier si nécessaire
    if (mkdir(folder_name, 0777) == -1) {
        if (errno == EEXIST)
            printf("Le dossier '%s' existe déjà.\n", folder_name);
        else
            perror("Erreur lors de la création du dossier");
    } else {
        printf("Le dossier '%s' a été créé avec succès.\n", folder_name);
    }

    // Prépare le fichier de sortie
    char file_path[256];
    snprintf(file_path, sizeof(file_path), "%s/msd_file_%d.txt", folder_name, nbr_simu);
    FILE *msd_file = fopen(file_path, "w");
    if (!msd_file) {
        fprintf(stderr, "Erreur : impossible d’ouvrir %s\n", file_path);
        return;
    }

    printf("[Simulation %d] Calcul de la MSD avec OpenMP (%d threads)\n",
           nbr_simu, omp_get_max_threads());

    // --- Boucle principale ---
    for (int i = 0; i < Tp; i++) {
        fprintf(msd_file, "%d ", i);

        double msd_total = 0.0;

        // Parallélisation sur les monomères
        #pragma omp parallel for reduction(+:msd_total) schedule(static)
        for (int k = 0; k < Nm; k++) {

            double msd = 0.0;

            for (int j = 0; j < Tp - i; j++) {
                msd += distance2(
                    R_monomere_arrays[nbr_simu][i + j][k],
                    R_monomere_arrays[nbr_simu][j][k]);
            }

            msd /= (Tp - i);   // moyenne temporelle

            // Stockage thread-safe : chaque thread écrit dans sa propre case [k]
            stock_msd[nbr_simu][i][k] = msd;

            #pragma omp critical  // écriture séquentielle dans le fichier
            fprintf(msd_file, "%f ", msd);

            msd_total += msd;
        }

        // Moyenne sur tous les monomères
        fprintf(msd_file, "%f\n", msd_total / Nm);
    }

    fclose(msd_file);
    printf("✅ Fichier sauvegardé : %s\n", file_path);
}

// Function to check if a number is already in the list
bool is_in_list(int *list, int size, int num) {
    for (int i = 0; i < size; i++) {
        if (list[i] == num) {
            return true;
        }
    }
    return false;
}



void calculate_msd_serial(double ****R, double ***stock_msd,
                          int Tp, int Nm, int nbr_simu)
{
    for (int i = 0; i < Tp; i++) {
        double msd_total = 0.0;
        for (int k = 0; k < Nm; k++) {
            double msd = 0.0;
            for (int j = 0; j < Tp - i; j++) {
                msd += distance2(R[nbr_simu][i + j][k], R[nbr_simu][j][k]);
            }
            msd /= (Tp - i);
            stock_msd[nbr_simu][i][k] = msd;
            msd_total += msd;
        }
        stock_msd[nbr_simu][i][0] = msd_total / Nm;
    }
}


void calculate_msd_parallel(double ****R, double ***stock_msd,
                            int Tp, int Nm, int nbr_simu)
{
    for (int i = 0; i < Tp; i++) {
        double msd_total = 0.0;

        #pragma omp parallel for reduction(+:msd_total) schedule(static)
        for (int k = 0; k < Nm; k++) {
            double msd = 0.0;
            for (int j = 0; j < Tp - i; j++) {
                msd += distance2(R[nbr_simu][i + j][k], R[nbr_simu][j][k]);
            }
            msd /= (Tp - i);
            stock_msd[nbr_simu][i][k] = msd;
            msd_total += msd;
        }

        stock_msd[nbr_simu][i][0] = msd_total / Nm;
    }
}




// Function to generate a list of unique random numbers
void generate_unique_random_numbers(int *list, int N, int count) {
    int index = 0;
    while (index < count) {
        int num = rand() % N;
        if (!is_in_list(list, index, num)) {
            list[index] = num;
            index++;
        }
    }
}

// Fonction pour calculer la somme
double somme(double* data, int taille) {
    double sum = 0.0;
    for (int i = 1; i < taille; i++) {
        sum += data[i];
    }
    
    return sum;
}

// Fonction de régression linéaire
void regression_lineaire(double* x, double* y, int taille, double* alpha, double* logA) {
    double sumX = somme(x, taille);
    double sumY = somme(y, taille);
    double sumXY = 0.0;
    double sumXX = 0.0;

    for (int i = 1; i < taille; i++) {
        sumXY += x[i] * y[i];
        sumXX += x[i] * x[i];
    }
    double denominateur = taille * sumXX - sumX * sumX;
    

    *alpha = (taille * sumXY - sumX * sumY) / denominateur;
    *logA = (sumY * sumXX - sumX * sumXY) / denominateur;
}

void print_time_remaining(int current_iteration, int total_iterations, time_t start_time){
    time_t current_time = time(NULL);
    double elapsed_time = difftime(current_time, start_time);

    double avg_time_per_iteration = elapsed_time / (current_iteration + 1);

    double remaining_time = avg_time_per_iteration * (total_iterations - current_iteration - 1);

    int remaining_minutes = (int) remaining_time / 60;
    int remaining_seconds = (int) remaining_time % 60;

    printf("Temps restant estime : %d min %d sec \n", remaining_minutes, remaining_seconds);
}

void plot_cols(double **array, int rows, int cols, int column_index, const char *filename, const char *title, const char *xlabel, const char *ylabel) {
    if (column_index >= cols) {
        fprintf(stderr, "Error: column index out of bounds\n");
        return;
    }

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error (plot_cols): Could not open file %s for writing\n", filename);
        return;
    }

    for (int i = 0; i < rows; i++) {
        fprintf(file, "%f\n", array[i][column_index]);
    }
    fclose(file);

    // Plot data using gnuplot
    FILE *gnuplotPipe = popen("gnuplot -persistent","w"); 
    if (gnuplotPipe){
        fprintf(gnuplotPipe, "set title '%s'\n", title); 
        fprintf(gnuplotPipe, "set xlabel '%s'\n", xlabel); 
        fprintf(gnuplotPipe, "set ylabel '%s'\n", ylabel); 
        fprintf(gnuplotPipe, "plot '%s' with linespoints title 'Data'\n", filename);
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    } else {
        fprintf(stderr, "Error: Could not open gnuplot pipe\n");
    }
}

double calculate_angle(double *vec1, double *vec2){
    double dot_product_value = dot_product(vec1, vec2);
    double norm_v1 = norm(vec1);
    double norm_v2 = norm(vec2);
    
    return acos(dot_product_value / (norm_v1 * norm_v2));
}

void update_link_vectors(double **r_new,double **t_link, int N){
     for (int i = 1; i < N; i++){
        for(int j = 0; j < 3; j++){
            t_link[i-1][j] = r_new[i][j] - r_new[i-1][j];
        }
    }
}

void compteur_grands_deplacements(int N, int T, double **R, double **r_new, int compteur){
    compteur = 0;
    for(int i = 0; i < 0; i++){
        if(distance(R[i], r_new[i]) > 0.1){
            compteur++;
        }
        
    }
    
}

void create_histogram(double *data, int data_size, int *bins, int bin_count, double min_value, double max_value){
    double bin_width = (max_value -min_value)/bin_count;

    // Initialisation des intervalles
    for(int i = 0; i < bin_count; i++){
        bins[i] = 0;
    }

    for(int i = 0; i < data_size; i++){
        if(data[i]>=min_value && data[i] < max_value){
            int bin_index = (int)((data[i]-min_value)/bin_width);
            bins[bin_index]++;
        }
    }

    FILE *histogram = fopen("histogram.txt","w");
    if (histogram == NULL) {
        fprintf(stderr, "Error: Could not open file for writing\n");
        return;
    }

    for(int i = 0; i < bin_count; i++){
       double bin_center = min_value + (i +0.5) * bin_width;
        double probability_density = (double)bins[i]/(data_size*bin_width);
        fprintf(histogram,"%.2f\t%f\n", bin_center, probability_density);
    }
    fclose(histogram);
}

bool is_present(int *liste, int taille, int valeur){
    for(int i = 0; i < taille; i++){
        if(liste[i] == valeur){
            return true;
        }
    }
    return false;
}

void echanger_deux_particules(double **R, int i, int j){
    double *temp = R[i];
    R[i] = R[j];
    R[j] = temp;
}

void reogarnisation(double **R, int N){
    FILE *fichier = fopen("test.txt", "w");

    for(int i = 0; i < N - 1; i++){
        double min_dist = 1000000000.0;
        int index_proche = -1; 

        for(int j = i + 1; j < N; j++){
            if(j!=N-1){
                double dist = distance(R[i],R[j]);
                if(dist < min_dist){
                    min_dist = dist; 
                    index_proche = j;
                }
            }
            
        }
        if(index_proche != -1){
            echanger_deux_particules(R, i+1, index_proche);
        }
    }

}

//////////// // Remplace data par la transformer de fourier discrete si Isign = 1, et inverse si isign = -1 ////////////////
void four1(double* data, int n, int Isign){ 
    int nn, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;
    if(n<2 || n&(n-1)) {
        fprintf(stderr, "Nombre de points n doit etre une puissance de 2 dans four1\n");
        exit(EXIT_FAILURE);
    }
    nn = n << 1;
    j = 1;
    for(i = 1; i < nn; i += 2) {
        if(j > i) {
            tempr = data[j]; data[j] = data[i]; data[i] = tempr;
            tempr = data[j+1]; data[j+1] = data[i+1]; data[i+1] = tempr;
        }
        m = n;
        while(m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    mmax = 2;
    while(nn > mmax) {
        istep = mmax << 1;
        theta = Isign * (6.28318530717959 / mmax);
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for(m = 1; m < mmax; m += 2) {
            for(i = m; i <= nn; i += istep) {
                j = i + mmax;
                tempr = wr * data[j] - wi * data[j+1];
                tempi = wr * data[j+1] + wi * data[j];
                data[j] = data[i] - tempr;
                data[j+1] = data[i+1] - tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
}
/////////////// Calcul la transformée de Fourier de n valeurs réelles, remplace ces valeurs par 
//////////////  la fréquence/2 de la transformée de Fourier, n doit être une puissance de 2
/////////////// Calcul la transformée de Fourier inverse si les datas sont deja transformée, dans ce cas 
/////////////// le resultat doit être multiplier par 2/n
void realft(double* vec_double, double* data, int Isign, int data_size){
    int i, i1, i2, i3, i4, n = data_size; 
    double c1 = 0.5, c2, h1r, h1i, h2r, h2i, wr, wi, wpr, wpi, wtemp;
    double theta = 3.141592653589793238 / (double) (n >> 1);
    if (Isign == 1) {
        c2 = -0.5;
        four1(data, n, 1);
    } else {
        c2 = 0.5;
        theta = -theta;
    }
    wtemp = sin(0.5 * theta);
    wpr = -2.0 * wtemp * wtemp;
    wpi = sin(theta);
    wr = 1.0 + wpr;
    wi = wpi;

    for(i = 1; i < (n >> 2); i++) {
        i2 = 1 + (i1 = i + i);
        i4 = 1 + (i3 = n - i1);
        h1r = c1 * (data[i1] + data[i3]);
        h1i = c1 * (data[i2] - data[i4]);
        h2r = -c2 * (data[i2] + data[i4]);
        h2i = c2 * (data[i1] - data[i3]);
        data[i1] = h1r + wr * h2r - wi * h2i;
        data[i2] = h1i + wr * h2i + wi * h2r;
        data[i3] = h1r - wr * h2r + wi * h2i;
        data[i4] = -h1i + wr * h2i + wi * h2r;
        wr = (wtemp = wr) * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
    }

    if (Isign == 1) {
        data[0] = (h1r = data[0]) + data[1];
        data[1] = h1r - data[1];
    } else {
        data[0] = c1 * ((h1r = data[0]) + data[1]);
        data[1] = c1 * (h1r - data[1]);
        four1(data, n, -1);
    }
}

void calculate_gyration_radius(double** R_centre_de_masse, double* gyration_radius, double** R, int N, int T_centre_de_masse, double** stock_cdm){
    for( int i = 0; i < T_centre_de_masse; i++){
        for(int j = 0; j < N; j++){
            gyration_radius[i]+= distance2(stock_cdm[i],R_centre_de_masse[i]); 
        }
        gyration_radius[i]/= N; 
    }
}


void save_gyration(double* gyration_radius, int nbr_simu, int T_centre_de_masse, int periode_centre_de_masse){
    const char* folder_name = "GYRATION_RADIUS"; 
    if(nbr_simu == 0){
        if(mkdir(folder_name, 0777) == -1){
            if( errno = EEXIST){
                printf("Dossier '%s' existe déja, \n", folder_name);
            }
            else{
                perror("Erreur lors de la création du dossier \n");
            }
        }
        else{
            printf("Dossier '%s' a été crée avec succès \n", folder_name);
        }

        char file_path[512];
        snprintf(file_path, sizeof(file_path), "%s/gyration_radius_%d.txt", folder_name, nbr_simu);

        FILE* file = fopen(file_path, "w");
        if(file == NULL){
            perror("Erreur lors de l'ouverture du fichier \n");
        }

        for(int i = 0; i < T_centre_de_masse; i++){
            double t = i * periode_centre_de_masse; 
            fprintf(file, "%f %f \n", t, gyration_radius[i]);
            //printf("rayon de giration = %f \n", gyration_radius[i]);
        }
        fclose(file);

    }
}

void creation_polymere_aleatoire(int N, double a, double **R){
    R[0][0]= genrand_real2() * 100;
    R[0][1]= genrand_real2() * 100;
    R[0][2]= genrand_real2() * 100;
    for(int i = 1; i < N; i++){
        double phi = genrand_real2() * 2 * PI;
        double theta = genrand_real2() * PI;
        R[i][0] = R[i-1][0] + a * sin(theta) * cos(phi);
        R[i][1] = R[i-1][1] + a * sin(theta) * sin(phi);
        R[i][2] = R[i-1][2] + a * cos(theta);
    }
}
