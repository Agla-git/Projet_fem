
/*
 *  fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#ifndef _FEM_H_
#define _FEM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gmshc.h"
#include <complex.h>


#define ErrorScan(a)   femErrorScan(a,__LINE__,__FILE__)
#define ErrorGmsh(a)   femErrorGmsh(a,__LINE__,__FILE__)
#define Error(a)       femError(a,__LINE__,__FILE__)
#define Warning(a)     femWarning(a,  __LINE__, __FILE__)
#define FALSE 0 
#define TRUE  1
#define MAXNAME 256

typedef enum {FEM_TRIANGLE,FEM_QUAD,FEM_EDGE} femElementType;
typedef enum {DIRICHLET_X,DIRICHLET_Y,NEUMANN_X,NEUMANN_Y} femBoundaryType;
typedef enum {PLANAR_STRESS,PLANAR_STRAIN,AXISYM} femElasticCase;


typedef struct {
    int nNodes;
    double *X;
    double *Y;
} femNodes;

typedef struct {
    int nLocalNode;
    int nElem;
    int *elem;
    femNodes *nodes;
} femMesh;

typedef struct {
    femMesh *mesh;
    int nElem;
    int *elem;
    char name[MAXNAME];
} femDomain;

typedef struct {
    //double LxPlate, LyPlate;
    double h;
    femElementType elementType;
    double (*geoSize)(double x, double y);
    femNodes *theNodes;
    femMesh  *theElements;
    femMesh  *theEdges;
    int nDomains;
    femDomain **theDomains;
    double h_Min;
    double h_Max;
    double d_Max;
    double rayonTrou;
    double xStart;
    double yPos;
    double espace_trous;
    double holePositions[3][2];
} femGeo;

typedef struct {
    int n;
    femElementType type;
    void (*x2)(double *xsi, double *eta);
    void (*phi2)(double xsi, double eta, double *phi);
    void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta);
    void (*x)(double *xsi);
    void (*phi)(double xsi, double *phi);
    void (*dphidx)(double xsi, double *dphidxsi);
} femDiscrete;
    
typedef struct {
    int n;
    const double *xsi;
    const double *eta;
    const double *weight;
} femIntegration;

typedef struct {
    double *B;
    double **A;
    int size;
} femFullSystem;


typedef struct {
    femDomain* domain;
    femBoundaryType type; 
    double value;
} femBoundaryCondition;


typedef struct {
    double E,nu,rho,g;
    double A,B,C;
    int planarStrainStress;
    int nBoundaryConditions;
    femBoundaryCondition **conditions;  
    int *constrainedNodes; 
    double *soluce;
    double *residuals;
    femGeo *geometry;
    femDiscrete *space;
    femIntegration *rule;
    femDiscrete *spaceEdge;
    femIntegration *ruleEdge;
    femFullSystem *system;
} femProblem;

typedef struct {
    double *B;
    double **A;        
    int size;
    int band;        
} femBandSystem;

typedef struct {
    double *R;
    double *D;
    double *S;
    double *X; 
    double error;      
    int size;
    int iter;        
} femIterativeSolver;

// --- Structures pour RCM aka optimisation de node placement 

// Pour la liste d'adjacence de chaque noeud
typedef struct femNodeAdj {
    int neighborIndex;        // Indice du noeud voisin
    struct femNodeAdj *next;  // Pointeur vers le voisin suivant
} femNodeAdj;

// Pour représenter le graphe du maillage
typedef struct {
    int nNodes;             // Nombre total de noeuds
    femNodeAdj **adjLists;  // Tableau de pointeurs vers les têtes des listes d'adjacence (taille nNodes)
    int *degree;            // Tableau des degrés de chaque noeud (taille nNodes)
} femGraph;


void femApplyNodePermutation(femGeo* geometry, int* permutationP, int* old_to_new);
int* femComputeOldToNewMap(int* permutationP, int nNodes);
int* femComputeRcmPermutation(femGraph* graph, int startNode);
int compareNeighbors(const void* a, const void* b);
int femFindMinDegreeNode(femGraph* graph);
void femFreeGraph(femGraph* graph);
femGraph* femBuildAdjacencyGraph(femMesh* elements, femMesh* edges, int nNodes);
void addEdge(femGraph* graph, int u, int v);


// --- Fin Structures pour RCM ---



void                geoInitialize();
femGeo*             geoGetGeometry();
double              geoSize(double x, double y);
double              geoSizeDefault(double x, double y);
double              computeMeshSize(double x, double y, double h_Min, double h_Max, double d_Max, double rayonTrou, double xStart, double yPos, double holePositions[3][2]);
void                geoSetSizeCallback(double (*geoSize)(double x, double y));
int                 count_lines(const char *filename);
//double            (*transform_joukovski(const char *filename, int num_lines))[2];
double              (*transform_NACA(const char *filename, int num_lines))[2];
void                geoMeshGenerate(const char *filename, int num_lignes, double rayonTrou, double xStart, double yPos);
void                geoMeshImport();
void                geoMeshPrint();
void                geoMeshWrite(const char *filename);
void                geoMeshRead(const char *filename);
void                geoSetDomainName(int iDomain, char *name);
int                 geoGetDomain(char *name);
void                geoFinalize();

femProblem*         femElasticityCreate(femGeo* theGeometry, 
                                      double E, double nu, double rho, double g, femElasticCase iCase);
void                femElasticityFree(femProblem *theProblem);
void                femElasticityPrint(femProblem *theProblem);
void                femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value);
void                femElasticityAssembleElements(femProblem *theProblem);
void                femElasticityAssembleNeumann(femProblem *theProblem);
double*             femElasticitySolve(femProblem *theProblem);
double*             femElasticityForces(femProblem *theProblem);
double              femElasticityIntegrate(femProblem *theProblem, double (*f)(double x, double y));


femIntegration*     femIntegrationCreate(int n, femElementType type);
void                femIntegrationFree(femIntegration *theRule);

femDiscrete*        femDiscreteCreate(int n, femElementType type);
void                femDiscreteFree(femDiscrete* mySpace);
void                femDiscretePrint(femDiscrete* mySpace);
void                femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta);
void                femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi);
void                femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);
void                femDiscreteXsi(femDiscrete* mySpace, double *xsi);
void                femDiscretePhi(femDiscrete* mySpace, double xsi, double *phi);
void                femDiscreteDphi(femDiscrete* mySpace, double xsi, double *dphidxsi);

femFullSystem*      femFullSystemCreate(int size);
void                femFullSystemFree(femFullSystem* mySystem);
void                femFullSystemPrint(femFullSystem* mySystem);
void                femFullSystemInit(femFullSystem* mySystem);
void                femFullSystemAlloc(femFullSystem* mySystem, int size);
double*             femFullSystemEliminate(femFullSystem* mySystem);
void                femFullSystemConstrain(femFullSystem* mySystem, int myNode, double value);

double              femMin(double *x, int n);
double              femMax(double *x, int n);
void                femError(char *text, int line, char *file);
void                femErrorScan(int test, int line, char *file);
void                femErrorGmsh(int test, int line, char *file);
void                femWarning(char *text, int line, char *file);

femBandSystem* femBandSystemCreate(int size, int band); 
void           femBandSystemFree(femBandSystem *myBandSystem); 
void           femBandSystemInit(femBandSystem *myBandSystem); 
void           femBandSystemPrint(femBandSystem *myBand); 
void           femBandSystemPrintInfos(femBandSystem *myBand); 
double         femBandSystemGet(femBandSystem* myBandSystem, int myRow, int myCol); 
double*        femBandSystemSolve(femBandSystem *sys);

void femIterativeSolverInit(femIterativeSolver *mySolver); 





#endif
