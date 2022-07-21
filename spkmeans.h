#ifndef SPKMEANS_H
#define SPKMEANS_H

double **run_algo(char *file_name, char *goal, int *N, int *k);
double **run_spkmeans(double **data, int N, int d, int *k, char *goal);
void run_kmeans(double **T, double** centroids, int N, int d, int k);
void free_matrix(double **matrix, int N);

#endif