#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define INT_MIN (pow(-2, 31))
#define EPSILON (pow(10,-15))
int const MAX_ITER_JACOBI = 100;
int const MAX_ITER_KMEANS = 300;

/***************************************** DECLERATIONS *****************************************/

typedef struct eigenValue{
   double value;
   int index;
}eigenValue;


double **process_data(int *d, int *N, char *FILE_NAME);

double **compute_wam(double **data, int N, int d);
double **compute_ddg(double **wam, int N);
void compute_sqrt_diagonal(double **ddg, int N);
double **compute_lnorm(double **wam, double** sqrt_ddg, int N);

double **Jacobi(double **A, int N);
eigenValue *generate_sorted_eigenValues (double **A, int N);
double **compute_U(double **V, eigenValue *eigenValues, int N, int k);
double **compute_T(double **U, int N, int k);

int compare(const void *a, const void *b);

double calculate_off_A(double **A, int N);
void update_V(double **V, int p, int q, double s, double c, int N);
double **Identity_matrix(int N);
void pivot(double **A, int *max_indices, int N);
double calc_theta(double **A, int i, int j);
double calc_t(double theta);
double calculte_norm2(double *vec1, double *vec2, int d);

double sum_of_pows(double *vec, int d);
double sum_row(double *row, int N);
void print_jacobi(double **eigenVectors, double **eigenValues, int N);

int eigen_gap_heuristic(eigenValue *eigenValues, int N);
double **initiate_centroids(double **T, int N, int k);
double **run_algo(char *file_name, char *goal, int *N, int *k);
double **run_spkmeans(double **data, int N, int d, int *k, char *goal);

void run_kmeans(double **T, double** centroids, int N, int d, int k);
int assign(double vec[], double** centroids, int K, int d);
double calculate_distance(double vec1[], double vec2[], int d);
void update_centroids(double** cSum, int* cNum, double **centroids, int K, int d); 
void sum_vec(double vec1[], double vec2[], int d);
void sub_vec(double vec1[], double vec2[], int d);
void print_matrix(double **matrix, int N, int k, int transpose);
void free_kmeans_data(int N, int K, double** T,double** centroids,double** cSum,int* cNum,int* T_index);

void free_all(int N, double **wam, double **ddg, double **lnorm, double **V, eigenValue* eigenValues);
void free_matrix(double **matrix, int N);

/***************************************** PROCESS DATA *****************************************/
double **process_data(int *d, int *N, char *FILE_NAME){
   char c;
   int scanf_retvalue= 1;
   double n1, *vec, **data;
   FILE *file;
   int len_max=128, cur_vec_size, cur_data_size, i=0, j=0, first_line=1;
   cur_vec_size = len_max;
   cur_data_size = len_max;
   file = fopen(FILE_NAME, "r");
   data =(double**) malloc(cur_data_size*sizeof(double*));
   if(data==NULL){
       printf("An Error Has Occured");
       assert(data != NULL);
    }  
    vec = (double *) malloc(cur_vec_size*sizeof(double));
    if(vec==NULL){
        printf("An Error Has Occured");
        assert(vec != NULL);
    }  
    while((scanf_retvalue = fscanf(file, "%lf%c", &n1, &c)) > 0){
      vec[i++] = n1;
      if(first_line){
         if(i == cur_vec_size){
            cur_vec_size = i + len_max;
            vec = (double*) realloc(vec, cur_vec_size*sizeof(double));
            if(vec==NULL){
               printf("An Error Has Occured");
               assert(vec != NULL);
            }         
        }
        if(scanf_retvalue < 2 || c == '\n'){
            *d = i;
            vec = (double*) realloc(vec, *d*sizeof(double));
            if(vec==NULL){
               printf("An Error Has Occured");
               assert(vec != NULL);
            }   
            data[j++] = vec;
            i = 0;
            first_line = 0;
            vec = (double*) malloc(*d*sizeof(double));
            if(vec==NULL){
               printf("An Error Has Occured");
               assert(vec != NULL);
            }   
         }
      }else{
         if(j == cur_data_size){
            cur_data_size += len_max;
            data = (double**) realloc(data, cur_data_size*sizeof(double*));
            if(data==NULL){
               printf("An Error Has Occured");
               assert(data != NULL);
            }   
         }
         if(scanf_retvalue < 2 || c == '\n'){
            data[j++] = vec;
            i = 0;
            vec = (double*) malloc(*d*sizeof(double));
            if(vec==NULL){
               printf("An Error Has Occured");
               assert(vec != NULL);
            }   
         }
      }
   }
   fclose(file);
   free(vec);  
   *N = j;
   return data;
}


/***************************************** MATRIX CALCULATION *****************************************/
double **compute_wam(double **data, int N, int d){
   int i, j;
   double norm;
   double **wam;
   wam = (double **) malloc(N*sizeof(double*));
   if(wam == NULL){
         printf("An Error Has Occured");
         assert(wam != NULL);
      }
   for(i=0; i<N; i++){
      wam[i] =(double*)calloc(N, sizeof(double));
      if(wam[i] == NULL){
         printf("An Error Has Occured");
         assert(wam[i] != NULL);
      }
   }
   for(i=0; i<N; i++){
      for(j=i+1; j<N; j++){
         norm = calculte_norm2(data[i], data[j], d);
         wam[i][j] = exp(-norm/2.0);
         wam[j][i] = wam[i][j];
      }
   }
   return wam;
}

double **compute_ddg(double **wam, int N){
   int i;
   double **ddg;
   ddg = (double **) malloc(N*sizeof(double*));
   if(ddg == NULL){
      printf("An Error Has Occured");
      assert(ddg != NULL);
   }
   for(i=0; i<N; i++){
      ddg[i] =(double*)calloc(N, sizeof(double));
      if(ddg[i] == NULL){
         printf("An Error Has Occured");
         assert(ddg[i] != NULL);
      }
   }
   for(i=0; i<N; i++){
      ddg[i][i] = sum_row(wam[i], N);
   }
   return ddg;
}

void compute_sqrt_diagonal(double **ddg, int N){
   int i;
   for(i=0; i<N; i++){
      ddg[i][i] = 1/sqrt(ddg[i][i]);
   }
   return;
}

double **compute_lnorm(double **wam, double** sqrt_ddg, int N){
   int i, j;
   double** lnorm;
   lnorm = (double **) malloc(N*sizeof(double*));
   if(lnorm == NULL){
      printf("An Error Has Occured");
      assert(lnorm != NULL);
   }
   for(i=0; i<N; i++){
      lnorm[i] =(double*) malloc(N*sizeof(double));
      if(lnorm[i] == NULL){
         printf("An Error Has Occured");
         assert(lnorm[i] != NULL);
      }
   }
   for(i=0; i<N; i++){
      lnorm[i][i] = 1;
      for(j=i+1; j<N; j++){
         lnorm[i][j] = -1*wam[i][j]*sqrt_ddg[i][i]*sqrt_ddg[j][j];
         lnorm[j][i] = lnorm[i][j];
      }
   }
   return lnorm;
}

/**************** JACOBI ALGORITHM ****************/
double **Jacobi(double **A, int N){
   double theta, t, c, s, **V, off_A, off_A_new=0,tmp;
   int i, p, q, max_indices[2], iteration=0;
   V = Identity_matrix(N);
   off_A = calculate_off_A(A, N);
   while(((off_A - off_A_new)>EPSILON) && (iteration<MAX_ITER_JACOBI)){
      if(iteration==0){
         off_A_new = off_A;
      }
      off_A = off_A_new;
      pivot(A, max_indices, N);
      p = max_indices[0], q = max_indices[1];
      if((fabs(A[p][q]) == 0) && (iteration ==0)){
         break;
      } 
      theta = calc_theta(A, p, q);
      t = calc_t(theta);
      c = 1/(sqrt(pow(t, 2)+1));
      s = t*c;
      if(iteration == 0){
         V[p][p] = V[q][q] = c;
         V[p][q] = s;
         V[q][p] = -s;
      }
      for(i=0; i<N; i++){
         if(i!=p && i!=q){
            off_A_new = off_A_new - 2*pow(A[i][p],2) + 2*pow((c*A[i][p] - s*A[i][q]),2);
            tmp = A[i][p];
            A[i][p] = c*A[i][p] - s*A[i][q];
            A[p][i] = A[i][p];

            off_A_new = off_A_new - 2*pow(A[i][q],2) + 2*pow((c*A[i][q] + s*tmp),2);
            A[i][q] = c*A[i][q] + s*tmp;
            A[q][i] = A[i][q];
         }
      }
      tmp = A[p][p];
      A[p][p] = pow(c, 2)*A[p][p] + pow(s, 2)*A[q][q] - 2*s*c*A[p][q];
      A[q][q] = pow(s, 2)*tmp + pow(c,2)*A[q][q] + 2*s*c*A[p][q];
      off_A_new = off_A_new - 2*pow(A[p][q],2);
      A[p][q] = 0;
      A[q][p] = 0;
      if(iteration != 0){
         update_V(V, p, q, s, c, N);
      }
      iteration++;
   }
   return V;
}

eigenValue *generate_sorted_eigenValues(double **A, int N){
   eigenValue *eigenValues;
   int i;
   eigenValues = (eigenValue*)malloc(N*sizeof(eigenValue));
   if(eigenValues==NULL){
      printf("An Error Has Occured");
      assert(eigenValues != NULL);
   }
   for(i=0; i<N; i++){
      eigenValues[i].value = A[i][i]; 
      eigenValues[i].index = i;
   }
   qsort(eigenValues, N, sizeof(eigenValue), compare);
   return eigenValues;
}

double **compute_U(double **V, eigenValue *eigenValues, int N, int k){
   double **U;
   int i, j;
   U = (double **) malloc(N*sizeof(double *));
   if(U == NULL){
      printf("An Error Has Occured");
      assert(U != NULL);
   }
   for(i=0; i<N; i++){
      U[i] = (double *)malloc(k*sizeof(double));
      if(U[i] == NULL){
         printf("An Error Has Occured");
         assert(U[i] != NULL);
      }
   }
   for(j=0; j<k; j++){
      for(i=0; i<N; i++){
         U[i][j] = V[i][eigenValues[j].index];
      }
   }
   return U;
}

double **compute_T(double **U, int N, int k){
   int i, j;
   double row_norm;
   for(i=0; i<N; i++){
      row_norm = sqrt(sum_of_pows(U[i], k));
      for(j=0; j<k; j++){
         U[i][j] = U[i][j]/row_norm;
      }
   }
   return U;
}

int compare(const void *a, const void *b){
   int i_a, i_b;
   double value_a, value_b; 
   value_a = ((eigenValue *)a)->value;
   value_b = ((eigenValue *)b)->value;
   i_a = ((eigenValue *)a)->index;
   i_b = ((eigenValue *)b)->index;
   if(value_a < value_b){
      return -1;
   }
   if(value_b < value_a){
      return 1;
   }
   if(i_a < i_b){
      return -1;
   }
   return 1;
}


/*********************************************** HELPER FUNCTIONS ***********************************************/
double calculate_off_A(double **A, int N){
   int i, j;
   double off_A=0;
   for(i=0; i<N; i++){
      for(j=i+1; j<N; j++){
         off_A += 2*pow(A[i][j],2);
      }
   }
   return off_A;
}

void update_V(double **V, int p, int q, double s, double c, int N){
   int i;
   double tmp;
   for(i=0; i<N; i++){
      tmp = V[i][p];
      V[i][p] = c*V[i][p]-s*V[i][q];
      V[i][q] = c*V[i][q] + s*tmp;
   }
   return;
}

/*Generates an identity matrix, and returns pointer to the matrix*/
double **Identity_matrix(int N){
   double **I;
   int i;
   I = (double **) malloc(N*sizeof(double*));
   if(I == NULL){
      printf("An Error Has Occured");
      assert(I != NULL);
   }
   for(i=0; i<N; i++){
      I[i] =(double*)calloc(N, sizeof(double));
      if(I[i] == NULL){
         printf("An Error Has Occured");
         assert(I[i] != NULL);
      }
   }
   for(i=0; i<N; i++){
      I[i][i] = 1;
   }
   return I;
}

void pivot(double **A, int *max_indices, int N){
   double max;
   int i = 0, j = 1;
   max_indices[0] = i;
   max_indices[1] = j;
   max = fabs(A[0][1]);
   for(i=0; i<N; i++){
      for(j=i+1; j<N; j++){
         if(fabs(A[i][j]) > max){
            max = fabs(A[i][j]);
            max_indices[0] = i, max_indices[1] = j;
         }
      }
   }
   return;
}

double calc_theta(double **A, int i, int j){
   double theta;
   theta = (A[j][j]-A[i][i])/(2*A[i][j]);
   return theta;
}

double calc_t(double theta){
   double t, sign; 
   if(theta==0){
      sign = 1;
   }
   else{
      sign = theta > 0 ? 1 : -1;
   }
   t = sign/(fabs(theta) +sqrt(pow(theta,2)+1));
   return t;
}

double calculte_norm2(double *vec1, double *vec2, int d){
   int i;
   double norm2, *sub;
   sub = (double *) malloc(d*sizeof(double));
   if(sub == NULL){
      printf("An Error Has Occured");
      assert(sub != NULL);
   }
   for(i=0; i<d; i++){
      sub[i] = vec1[i] - vec2[i];
   }
   norm2 = sqrt(sum_of_pows(sub, d));
   free(sub);
   return norm2;
}

/*********************************************** HELPER FUNCTIONS - MATRIX CALCULATION ***********************************************/

double sum_of_pows(double *vec, int d){
   int i;
   double result=0;
   for(i=0; i<d; i++){
      result += pow(vec[i], 2);
   }
   return result;
}

double sum_row(double *row, int N){
   int i;
   double sum=0;
   for(i=0; i<N; i++){
      sum += row[i];
   }
   return sum;
}

void print_matrix(double **matrix, int N, int k, int transpose){
   int i, j;
   double x;
   for(i=0; i<N; i++){
      for(j=0; j<k; j++){
         x = transpose ? matrix[j][i]:matrix[i][j];
         printf("%.4f", (x>-0.00005 && x<0) ? 0 : x);
         if(j<k-1){
            printf(",");
         }
      }
      if(i != (N-1)){
         printf("\n");
      }
   }
}

void print_jacobi(double **eigenVectors, double **eigenValues, int N){
   int i;
   double x;
   for(i=0; i<N; i++){
      x = eigenValues[i][i];
      printf("%.4f", (x>-0.00005 && x<0) ? 0 : x);
         if(i<N-1){
            printf(",");
         }
   }
   printf("\n");
   print_matrix(eigenVectors, N, N, 1);
}


/**************** EIGEN VALUE HURESTIC ****************/
int eigen_gap_heuristic(eigenValue *eigenValues, int N){
   /*check case k = 0*/
   int i, gap_index = 0;
   double max_gap, new_gap;
   max_gap = INT_MIN;
   for(i=0; i<(N/2)-1; i++){
      new_gap = fabs(eigenValues[i].value - eigenValues[i+1].value);
      if(new_gap > max_gap){
         max_gap = new_gap;
         gap_index = i;
      }
   }
   return gap_index+1;
}

/*********************************************** SPKMEANS ***********************************************/

double **initiate_centroids(double **T, int N, int K){
   int i,j;
   double** centroids;
   centroids = (double**) malloc(N*sizeof(double*));
   if(centroids == NULL){
         printf("An Error Has Occured");
         assert(centroids != NULL);
      }
   for(i=0; i<K; i++){
      centroids[i]= (double*)malloc(K*sizeof(double));
      if(centroids[i] == NULL){
         printf("An Error Has Occured");
         assert(centroids[i] != NULL);      }
   }
   for(i=0; i < K ; i++){
      for(j=0; j<K; j++){
         centroids[i][j] = T[i][j];
      }
   }
   return centroids;
}

/*this function runs the algorithm until T*/
double **run_algo(char *file_name, char *goal, int *N, int *k){
   int d;
   double **data, **T;
   data = process_data(&d, N, file_name);
   if(*N<=*k && (strcmp(goal,"spk") == 0)){
      printf("Invalid Input!");
      return NULL;
   }
   T = run_spkmeans(data, *N, d, k, goal);
   free_matrix(data, *N);
   return T;
}

double **run_spkmeans(double **data, int N, int d, int *k, char *goal){
   double **wam, **ddg, **V, **lnorm, **U, **T;
   eigenValue *eigenValues;
   if(!strcmp(goal, "jacobi")){
      /*NOTE: if goal==jacobi, data is A_0*/
      V = Jacobi(data, N);
      print_jacobi(V, data, N);
      free_matrix(V, N);
      return NULL;
   }
   wam = compute_wam(data, N, d);
   if(!strcmp(goal, "wam")){
      print_matrix(wam, N, N,0);
      free_matrix(wam,N);
      return NULL;
   }
   ddg = compute_ddg(wam, N);
   if(!strcmp(goal, "ddg")){
      print_matrix(ddg, N, N, 0);
      free_matrix(wam, N);
      free_matrix(ddg, N);
      return NULL;
   }
   compute_sqrt_diagonal(ddg, N);
   lnorm = compute_lnorm(wam, ddg, N);
   if(!strcmp(goal, "lnorm")){
      print_matrix(lnorm, N, N,0);
      free_matrix(wam, N);
      free_matrix(ddg, N);  
      free_matrix(lnorm,N);    
      return NULL;
   }

   /*NOTE: goal == spk*/
   V = Jacobi(lnorm, N);
   eigenValues = generate_sorted_eigenValues(lnorm, N);
   if(*k==0){
      *k=eigen_gap_heuristic(eigenValues, N);
   }
   U = compute_U(V, eigenValues, N, *k);
   T = compute_T(U, N, *k); 
   free_all(N, wam, ddg, lnorm, V, eigenValues);
   return T;
}

/*********************************************** KMEANS ***********************************************/

void run_kmeans(double **T, double**centroids, int N, int d, int K){
    int counter=0, changed=1, i, *cNum, *T_index;
    double **cSum;
    cNum = (int*)calloc(K,sizeof(int));
    if(cNum==NULL){
        printf("An Error Has Occured");
        assert(cNum != NULL);
    } 
    cSum =(double**)malloc(K*sizeof(double*));
    if(cSum==NULL){
        printf("An Error Has Occured");
        assert(cSum != NULL);
    } 
    for(i=0; i<K; i++){
        cSum[i] =(double*)calloc(d,sizeof(double));
        if(cSum[i]==NULL){
            printf("An Error Has Occured");
            assert(cSum[i] != NULL);
        }
    }
    T_index = (int*)malloc(N*sizeof(int));
    if(T_index==NULL){
        printf("An Error Has Occured");
        assert(T_index != NULL);
    } 
    for(i=0; i<N; i++){
        T_index[i] = -1;
    }
    while((counter < MAX_ITER_KMEANS) & changed){
        changed = 0;
        counter += 1;
        for(i=0; i<N; i++){
            int index = assign(T[i],centroids, K, d);
            if(T_index[i] != index){
                changed = 1;
                cNum[index] += 1;
                sum_vec(cSum[index], T[i], d);
                if(T_index[i] != -1){
                    cNum[T_index[i]] -= 1;
                    sub_vec(cSum[T_index[i]], T[i], d);
                }
            }T_index[i] = index;
        }update_centroids(cSum, cNum, centroids, K, d);  
    }
    print_matrix(centroids, K, K,0);
    free_kmeans_data(N, K, T, centroids, cSum, cNum, T_index);
    return;
}


int assign(double vec[], double** centroids, int K, int d){
    double min = 0;
    int i, centroid_index = -1;
    for(i=0; i<K; i++){
        double new_dis = calculate_distance(vec, centroids[i], d);
        if((new_dis < min) || (centroid_index == -1)){
            min = new_dis;
            centroid_index = i;
        }
    }return centroid_index;
}

void update_centroids(double** cSum, int* cNum, double** centroids, int K, int d){
    int i, j;
    for(i=0; i<K; i++){
        double * centroid = centroids[i];
        for(j=0; j<d; j++){
            if(cNum[i] != 0){
                centroid[j] = cSum[i][j]/cNum[i];
            }
        }
    }
}

/*********************************************** MATRIX CALCULATION ***********************************************/
double calculate_distance(double vec1[], double vec2[], int d){
    double sum = 0;
    int i;
    for(i=0; i<d; i++){
        sum += (vec1[i]-vec2[i])*(vec1[i]-vec2[i]);
    }return sum;
}

void sum_vec(double vec1[], double vec2[], int d){
    int i;
    for(i=0; i<d; i++){
        vec1[i] += vec2[i];
    }
}

void sub_vec(double vec1[], double vec2[], int d){
    int i;
    for(i=0; i<d; i++){
        vec1[i] -= vec2[i];
    }
}

/*********************************************** MEMORY ***********************************************/
void free_matrix(double **matrix, int N){
    int i;
    for(i=0 ; i<N ; i++){
        free(matrix[i]);
    }free(matrix);
    return;
}

void free_kmeans_data(int N, int K, double** T,double** centroids,double** cSum,int* cNum,int* T_index){
    int i;
    for(i=0 ; i<N ; i++){
        free(T[i]);
    }free(T);
    for(i=0; i<K; i++){
        free(centroids[i]);
        free(cSum[i]);
    }
    free(centroids);
    free(cSum);
    free(cNum);
    free(T_index);
}

void free_all(int N, double **wam, double **ddg, double **lnorm, double **V, eigenValue* eigenValues){
   int i;
   for(i=0 ; i<N ; i++){
      free(wam[i]);
      free(ddg[i]);
      free(lnorm[i]);
      free(V[i]);
   }
   free(wam);
   free(ddg);
   free(lnorm);
   free(V);
   free(eigenValues);
   return;
}


/*********************************************** MAIN ***********************************************/
int main(int argc, char* argv[]){
   int k, N, i;
   double **T, **centroids;
   char *goal;
   char *goals[] = {"wam", "ddg", "lnorm", "jacobi", "spk"};
   if(argc != 4){
      printf("Invalid Input!");
      return -1;
   }
   k = atoi(argv[1]);
   goal = argv[2];
   for(i = 0; i<5; i++){
      if(strcmp(goal, goals[i]) == 0){
         break;
      }
   }
   if(i==5){
      printf("Invalid Input!");
      return -1;
   }
   T = run_algo(argv[3], goal, &N, &k);
  if(T != NULL){
      centroids = initiate_centroids(T, N, k);
      run_kmeans(T, centroids, N, k, k);
   }
   return 1;
}

