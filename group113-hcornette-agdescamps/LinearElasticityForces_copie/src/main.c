/*
 *  main.c
 *  Library for EPL1110 : Finite Elements for dummies
 *  Elasticite lineaire plane
 *  Calcul des densités de force aux noeuds contraints
 *
 *  Copyright (C) 2024 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "glfem.h"

double fun(double x, double y) 
{
    return 1;
}

int main(void)
{  
    printf("\n\n    V : Mesh and displacement norm \n");
    printf("    D : Domains \n");
    printf("    X : Horizontal residuals for unconstrained equations \n");
    printf("    Y : Horizontal residuals for unconstrained equations \n");
    printf("    N : Next domain highlighted\n\n\n");


    //double Lx = 1.0;
    //double Ly = 1.0;
      
    geoInitialize();
    
    femGeo* theGeometry = geoGetGeometry();
    
    //theGeometry->LxPlate     =  Lx;
    //theGeometry->LyPlate     =  Ly;     
    theGeometry->h           =  1 * 0.075;    
    theGeometry->elementType = FEM_TRIANGLE;

    //Déclaration des variables globales
    theGeometry->h_Min = 0.005;
    theGeometry->h_Max;
    theGeometry->d_Max = 0.1;
    theGeometry->rayonTrou = 0.025;
    theGeometry->xStart = 0.2;
    theGeometry->yPos = 0.01;
    theGeometry->espace_trous = 0.2;
    theGeometry->holePositions[3][2];


   
    theGeometry->h_Max = theGeometry->h; // h_Max dépend de theGeometry.h
    for (int i = 0; i < 3; i++) {
        theGeometry->holePositions[i][0] = theGeometry->xStart + i * theGeometry->espace_trous;
        theGeometry->holePositions[i][1] = theGeometry->yPos;
    }


    const char *filename = "../NACA_points.csv";  // Vérifie l'orthographe correcte du fichier
    
    int num_lines_3 = count_lines(filename);  // Passe l'adresse de numPoints

    if (num_lines_3 <= 1) {
        printf("Erreur : Pas assez de points dans le fichier.\n");
        return 1;  // Quitter le programme en cas d'erreur
    }
    printf("Nombre de lignes : %d\n", num_lines_3);
    // Récupérer le tableau de points transformés
    double (*points_NACA)[2] = transform_NACA(filename, num_lines_3);

    if (!points_NACA) {
        printf("Erreur lors de la transformation de Joukovski.\n");
        return 1;
    }

    
    geoMeshGenerate(filename, num_lines_3,  theGeometry->rayonTrou, theGeometry->xStart, theGeometry->yPos);
    geoMeshImport();
    free(points_NACA);
    geoSetDomainName(0, "cercle1"); 
    geoSetDomainName(1,"cercle2");
    geoSetDomainName(2,"cercle3");


    for (int i = 0; i < theGeometry->nDomains; i++) {
    printf("Domain %d: %s\n", i, theGeometry->theDomains[i]->name);
    }
    

    geoMeshWrite("../data/elasticity.txt");
    
        
//
//  -2- Creation probleme 
//
    
    double E   = 211.e9;
    double nu  = 0.3;
    double rho = 7.85e3; 
    double g   = 9.81;
    //double g   = 0;
    femProblem* theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRAIN);
    
    femElasticityAddBoundaryCondition(theProblem,"cercle1",DIRICHLET_X,0.0);
    femElasticityAddBoundaryCondition(theProblem,"cercle1",DIRICHLET_Y,0.0);
    femElasticityAddBoundaryCondition(theProblem,"cercle2",DIRICHLET_X,0.0);
    femElasticityAddBoundaryCondition(theProblem,"cercle2",DIRICHLET_Y,0.0);
    femElasticityAddBoundaryCondition(theProblem,"cercle3",DIRICHLET_X,0.0);
    femElasticityAddBoundaryCondition(theProblem,"cercle3",DIRICHLET_Y,0.0);
    
    for (int i = 0; i < theGeometry->nDomains; i++) {
        if(i%2 == 0){
            femElasticityAddBoundaryCondition(theProblem,theGeometry->theDomains[i]->name,NEUMANN_Y,500.0);//intrados
        }
        if(i%2 !=0){
            femElasticityAddBoundaryCondition(theProblem,theGeometry->theDomains[i]->name,NEUMANN_Y,-4000.0);//extrados
        }
    }

    
    femElasticityPrint(theProblem);

//
//  -3- Resolution du probleme et calcul des forces
//

    double *theSoluce = femElasticitySolve(theProblem);
    double *theForces = femElasticityForces(theProblem);
    double area = femElasticityIntegrate(theProblem, fun);   
   
//
//  -4- Deformation du maillage pour le plot final
//      Creation du champ de la norme du deplacement
//
    
    femNodes *theNodes = theGeometry->theNodes;
    double deformationFactor = 1e5;
    double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
    double *forcesX = malloc(theNodes->nNodes * sizeof(double));
    double *forcesY = malloc(theNodes->nNodes * sizeof(double));
    
    for (int i=0; i<theNodes->nNodes; i++){
        theNodes->X[i] += theSoluce[2*i+0]*deformationFactor;
        theNodes->Y[i] += theSoluce[2*i+1]*deformationFactor;
        normDisplacement[i] = sqrt(theSoluce[2*i+0]*theSoluce[2*i+0] + 
                                   theSoluce[2*i+1]*theSoluce[2*i+1]);
        forcesX[i] = theForces[2*i+0];
        forcesY[i] = theForces[2*i+1]; }
  
    double hMin = femMin(normDisplacement,theNodes->nNodes);  
    double hMax = femMax(normDisplacement,theNodes->nNodes);  
    printf(" ==== Minimum displacement          : %14.7e [m] \n",hMin);
    printf(" ==== Maximum displacement          : %14.7e [m] \n",hMax);

//
//  -5- Calcul de la force globaleresultante
//

    double theGlobalForce[2] = {0, 0};
    for (int i=0; i<theProblem->geometry->theNodes->nNodes; i++) {
        theGlobalForce[0] += theForces[2*i+0];
        theGlobalForce[1] += theForces[2*i+1]; }
    printf(" ==== Global horizontal force       : %14.7e [N] \n",theGlobalForce[0]);
    printf(" ==== Global vertical force         : %14.7e [N] \n",theGlobalForce[1]);
    printf(" ==== Weight                        : %14.7e [N] \n", area * rho * g);

//
//  -6- Visualisation du maillage
//  
    
    int mode = 1; 
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[MAXNAME];
   
 
    GLFWwindow* window = glfemInit("EPL1110 : Recovering forces on constrained nodes");
    glfwMakeContextCurrent(window);

    do {
        int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theGeometry->theNodes,w,h);

        t = glfwGetTime();  
        if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
        if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
        if (glfwGetKey(window,'X') == GLFW_PRESS) { mode = 2;}
        if (glfwGetKey(window,'Y') == GLFW_PRESS) { mode = 3;}
        if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}
        if (t-told > 0.5) {freezingButton = FALSE; }
        
        if (mode == 0) {
            domain = domain % theGeometry->nDomains;
            glfemPlotDomain( theGeometry->theDomains[domain]); 
            sprintf(theMessage, "%s : %d ",theGeometry->theDomains[domain]->name,domain);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 1) {
            glfemPlotField(theGeometry->theElements,normDisplacement);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 2) {
            glfemPlotField(theGeometry->theElements,forcesX);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 3) {
            glfemPlotField(theGeometry->theElements,forcesY);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
         glfwSwapBuffers(window);
         glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed

    free(normDisplacement);
    free(forcesX);
    free(forcesY);
    femElasticityFree(theProblem) ; 
    geoFinalize();
    glfwTerminate(); 
    
    exit(EXIT_SUCCESS);
    return 0;  
}


