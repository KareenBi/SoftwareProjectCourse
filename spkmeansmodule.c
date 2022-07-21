#include <assert.h>
#include <Python.h>
#include "spkmeans.h"
#define PY_SSIZE_T_CLEANS

/* functions declarations */
static PyObject* spkmeans_capi(PyObject *self, PyObject *args);
static PyObject* kmeans_capi(PyObject *self, PyObject *args);
static double **list_to_matrix(PyObject* _list, int num_row, int num_col);
static PyObject *link_kmeans(PyObject* Py_data, PyObject* Py_centroids, int N, int d, int K); 
static PyObject *extract_T(const char* FILE_NAME, const char* goal,int k);
static PyObject* matrix_to_list(double ** matrix, int K, int d);


/*********************************************** C API ***********************************************/

/* the wrapping function for the run_algo in spkmeans.c - parses variables*/
static PyObject* spkmeans_capi(PyObject *self, PyObject *args){
    PyObject* result = NULL;
    const char *FILE_NAME, *goal;
    int k;
    if(!PyArg_ParseTuple(args, "ssi", &FILE_NAME, &goal, &k)){
       return NULL;
   }
   result = extract_T(FILE_NAME, goal, k);
   if(result == NULL){
       Py_RETURN_NONE;
   }
   return Py_BuildValue("O", result);
}


/* the wrapping function for the link_kmeans - parses PyObjects */
static PyObject* kmeans_capi(PyObject *self, PyObject *args){
    PyObject *data, *centroids;
    int N, d, K;
   if(!PyArg_ParseTuple(args, "OOiii", &data, &centroids, &N, &d, &K)){
       return NULL;
   }
   return Py_BuildValue("O", link_kmeans(data, centroids, N, d, K));
}

/* functino that parses the data and puts them in arrays */
static double **list_to_matrix(PyObject* _list, int num_row, int num_col) {
    int i, j;
    Py_ssize_t Py_i, Py_j;
    double **parsed_data;
    parsed_data = malloc(num_row * sizeof(double*));
    if(parsed_data == NULL){
        printf("An Error Has Occured");
        assert(parsed_data != NULL);
    }    
    PyObject* item; PyObject* num;
    for (i = 0; i < num_row; i++) {
        Py_i = (Py_ssize_t)i;
        parsed_data[Py_i] = malloc(num_col * sizeof(double));
        if(parsed_data[Py_i] == NULL){
            printf("An Error Has Occured");
            assert(parsed_data[Py_i] != NULL);
        }
        item = PyList_GetItem(_list, Py_i);
        if (!PyList_Check(item)){ /* Skips non-lists */
            continue;
        }
        for (j = 0; j < num_col; j++) {
            Py_j = (Py_ssize_t)j;
            num = PyList_GetItem(item, Py_j);
            if (!PyFloat_Check(num)) continue; /* Skips non-floats */
            parsed_data[Py_i][Py_j] = PyFloat_AsDouble(num);
        }
    }return parsed_data;
}

/* this array tells python what methods this module has */
static PyMethodDef capiMethods[] = {
                                    {"fit_spkmeans",
                                    (PyCFunction) spkmeans_capi,
                                    METH_VARARGS,
                                    PyDoc_STR("computes T if goal==spk, else returns NULL and prints appropriate output")},

                                    {"fit_kmeans",
                                    (PyCFunction) kmeans_capi,
                                    METH_VARARGS,
                                    PyDoc_STR("calculates the centroids using kmeans algorithm")},
                                    
                                    {NULL, NULL, 0, NULL}
        };  



/* This struct initiates the module using the above definition. */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,  
    "spkmeans",           
    NULL,                   
    -1,                     
    capiMethods
};

/* Module Creation */
PyMODINIT_FUNC
PyInit_spkmeans(void) {
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}
static PyObject *extract_T(const char* FILE_NAME, const char* goal,int k){
    int N, valid_k = k;
    double **T;
    PyObject* T_list;
    T = run_algo((char*)FILE_NAME, (char*)goal, &N, &valid_k);
    if(T == NULL){
        return NULL;
    }
    T_list = matrix_to_list(T, N, valid_k);
    free_matrix(T, N);
    return T_list;
}

/* initializing data and running the algorithm */
static PyObject *link_kmeans(PyObject* Py_data, PyObject* Py_centroids, int N, int d, int K){
    double **data, **centroids;
    /*data and centroids are freed in kmeans.c after run_kmeans().*/
    data = list_to_matrix(Py_data, N, d);
    centroids = list_to_matrix(Py_centroids, K, d);
    run_kmeans(data, centroids, N, d, K);
    Py_RETURN_NONE;
}

static PyObject* matrix_to_list(double ** matrix, int K, int d){
    int  i, j;
    PyObject *lst_obj, *vec, *num;
    lst_obj = PyList_New(K);
    if (!lst_obj){
        return NULL;
    }
    for(i=0; i<K; i++){
        vec = PyList_New(d);
        if (!vec){
            return NULL;
        }
        for (j = 0; j < d; j++) {
            num = PyFloat_FromDouble(matrix[i][j]);
            if (!num){
                Py_DECREF(vec);
                return NULL;
            }PyList_SET_ITEM(vec, j, num); 
        }PyList_SET_ITEM(lst_obj, i, vec); 
    }
    return lst_obj;
}



