
#include "glfem.h"
#include <time.h>


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

    double x_center_hole1 [8] = {0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2};
    double x_center_hole2 [8] = {0.3, 0.4, 0.5, 0.5, 0.6, 0.6, 0.6, 0.6, 0.6};
    double x_center_hole3 [8] = {0.5, 0.7, 0.8, 0.9, 1, 0.95, 0.9, 0.85, 0.825};
    // vu que le dernier trou dépasse on va diminuer l'ordonée ==> 2e colonne des tableau y_center_hole1,2,3
    // tout compte fait ça reste trop étroit mais je crois que les rayon sont déjà bien assez petits OK dcp je vais rester avant 0.85
    // et à cemoment-là c'est hyper fort déformé punaise. bon par grave je garde ça. Puis je peux monter qu'avec des plus grands rayons ça se déforme beaucoup moins. l'idéal ce serait de faire un plot python avec les déplacement max  et force globale max 

    double y_center_hole1 [3] = {0.01, 0.005, 0.01};
    double y_center_hole2 [3] = {0.01, 0.005, 0.01};
    double y_center_hole3 [3] = {0.01, 0.005, 0.01};

    int indice_center_x = 7;
    int indice_center_y = 1;


     femHole true_holes[3] = {
    {x_center_hole1[indice_center_x], y_center_hole1[indice_center_y], 0.01},  // Trou 1
    {x_center_hole2[indice_center_x], y_center_hole2[indice_center_y], 0.01},  // Trou 2
    {x_center_hole3[indice_center_x], y_center_hole3[indice_center_y], 0.01}   // Trou 3
};

   
    geoInitialize(true_holes);
    femGeo* theGeometry = geoGetGeometry();  
    
    //theGeometry->LxPlate     =  Lx;
    //theGeometry->LyPlate     =  Ly;     
    theGeometry->h           =  1 * 0.075;    
    theGeometry->elementType = FEM_TRIANGLE;

    //Déclaration des variables globales
    theGeometry->h_Min = 0.005;
    theGeometry->h_Max = theGeometry->h;
    theGeometry->d_Max = 0.1;


   
    
    const char *filename = "../NACA_points.csv";  
    
    int num_lines_3 = count_lines(filename);  

    if (num_lines_3 <= 1) {
        printf("Erreur : Pas assez de points dans le fichier.\n");
        return 1;  
    }
    printf("Nombre de lignes : %d\n", num_lines_3);
    // Récupérer le tableau de points transformés
    double (*points_NACA)[2] = transform_NACA(filename, num_lines_3);

    if (!points_NACA) {
        printf("Erreur lors de la transformation de Joukovski.\n");
        return 1;
    }

    
    geoMeshGenerate(filename, num_lines_3);
    geoMeshImport();

    clock_t begin = clock();  // Start timing here


    // --- Renumérotation RCM ---
    // printf("Starting RCM node renumbering...\n");
    femGraph* meshGraph = femBuildAdjacencyGraph(theGeometry->theElements, theGeometry->theEdges, theGeometry->theNodes->nNodes);
    if (meshGraph) {
        int startNode = femFindMinDegreeNode(meshGraph);
        int* rcmPermutationP = femComputeRcmPermutation(meshGraph, startNode);
        if (rcmPermutationP) {
            int* old_to_new = femComputeOldToNewMap(rcmPermutationP, theGeometry->theNodes->nNodes);
            if (old_to_new) {
                femApplyNodePermutation(theGeometry, rcmPermutationP, old_to_new);
                free(old_to_new);
                // printf("RCM renumbering applied successfully.\n");
            } else {
                printf("Error computing old_to_new map.\n");
            }
            free(rcmPermutationP);
        } else {
            printf("RCM permutation computation failed (check graph connectivity?).\n");
        }
        femFreeGraph(meshGraph); // Libérer le graphe après usage
    } else {
        printf("Failed to build adjacency graph.\n");
    }
    // --- Fin Renumérotation RCM ---




    free(points_NACA);
    geoSetDomainName(0, "cercle1"); 
    geoSetDomainName(1,"cercle2");
    geoSetDomainName(2,"cercle3");


    // for (int i = 0; i < theGeometry->nDomains; i++) {
    // printf("Domain %d: %s\n", i, theGeometry->theDomains[i]->name);
    // }
    

    geoMeshWrite("../data/elasticity.txt");
    
        
//
//  -2- Creation probleme 
//
    //Acier
    //double E   = 211.e9;
    //double nu  = 0.3;
    //double rho = 7.85e3; 

    //Aluminium
    double E = 60.e9; 
    double nu = 0.33; 
    double rho = 2.7e3;
    double g   = 9.81;
    //double g   = 0;
    femProblem* theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRESS);
    //femProblem* theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRAIN);
    
    femElasticityAddBoundaryCondition(theProblem,"cercle1",DIRICHLET_X,0.0);
    femElasticityAddBoundaryCondition(theProblem,"cercle1",DIRICHLET_Y,0.0);
    femElasticityAddBoundaryCondition(theProblem,"cercle2",DIRICHLET_X,0.0);
    femElasticityAddBoundaryCondition(theProblem,"cercle2",DIRICHLET_Y,0.0);
    femElasticityAddBoundaryCondition(theProblem,"cercle3",DIRICHLET_X,0.0);
    femElasticityAddBoundaryCondition(theProblem,"cercle3",DIRICHLET_Y,0.0);


    double pressures_intrados[] = {500.0, 750.0, 1000.0, 1250.0, 1500.0, 2000.0, 2250.0, };   // Pressions pour l'intrados
    double pressures_extrados[] = {-1000.0, -1500.0, -2000.0, -2500.0, -3000.0, -3250.0, -4000.0}; // Pressions pour l'extrados

   
    int numPressuresIntrados = sizeof(pressures_intrados) / sizeof(pressures_intrados[0]);
    int numPressuresExtrados = sizeof(pressures_extrados) / sizeof(pressures_extrados[0]);

   
    int currentPressureIntradosIndex = 2;
    int currentPressureExtradosIndex = 4;

    
    for (int i = 0; i < theGeometry->nDomains; i++) {
        if(i%2 == 0){
            femElasticityAddBoundaryCondition(theProblem,theGeometry->theDomains[i]->name,NEUMANN_Y,pressures_intrados[currentPressureIntradosIndex]);//intrados
        }
        if(i%2 !=0){
            femElasticityAddBoundaryCondition(theProblem,theGeometry->theDomains[i]->name,NEUMANN_Y,pressures_extrados[currentPressureExtradosIndex]);//extrados
        }
    }

    
    femElasticityPrint(theProblem);

//
//  -3- Resolution du probleme et calcul des forces
//

    double *theSoluce = femElasticitySolve(theProblem);
    double *theForces = femElasticityForces(theProblem);
    double area = femElasticityIntegrate(theProblem, fun);   

    clock_t end = clock();    // End timing here
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Time taken: %f seconds\n", time_spent);
   
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

    // Trouver l'index du noeud avec le déplacement minimal et maximal
    int minNodeIndex = 0;
    int maxNodeIndex = 0;
    for (int i = 1; i < theNodes->nNodes; i++) {
        if (normDisplacement[i] < normDisplacement[minNodeIndex]) {
            minNodeIndex = i;
        }
        if (normDisplacement[i] > normDisplacement[maxNodeIndex]) {
            maxNodeIndex = i;
        }
    }

    // Afficher les noeuds avec les déplacements min et max
    printf(" ==== Node with minimum displacement: Node %d, Displacement: %14.7e [m]\n", minNodeIndex, normDisplacement[minNodeIndex]);
    printf(" ==== Node with maximum displacement: Node %d, Displacement: %14.7e [m]\n", maxNodeIndex, normDisplacement[maxNodeIndex]);

    // Afficher les coordonnées des noeuds avec les déplacements min et max
    printf(" ==== Coordinates of Node with minimum displacement: (%14.7e, %14.7e) [m]\n", theNodes->X[minNodeIndex], theNodes->Y[minNodeIndex]);
    printf(" ==== Coordinates of Node with maximum displacement: (%14.7e, %14.7e) [m]\n", theNodes->X[maxNodeIndex], theNodes->Y[maxNodeIndex]);


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








/*
#include "glfem.h"
#include <time.h>
#include "pressure_manager.h"  
double fun(double x, double y) 
{
    return 1;
}

int enableAnimation = 0;  // 1 pour activer l'animation, 0 pour la désactiver

int main(void)
{  
    printf("\n\n    V : Mesh and displacement norm \n");
    printf("    D : Domains \n");
    printf("    X : Horizontal residuals for unconstrained equations \n");
    printf("    Y : Horizontal residuals for unconstrained equations \n");
    printf("    N : Next domain highlighted\n\n\n");

    femHole true_holes[3] = {
        {0.2, 0.01, 0.01},  // Trou 1
        {0.4, 0.01, 0.01},  // Trou 2
        {0.6, 0.01, 0.01}   // Trou 3
    };

    geoInitialize(true_holes);
    femGeo* theGeometry = geoGetGeometry();  
    
    theGeometry->h           = 1 * 0.075;    
    theGeometry->elementType = FEM_TRIANGLE;

    // Déclaration des variables globales
    theGeometry->h_Min = 0.005;
    theGeometry->h_Max = theGeometry->h;
    theGeometry->d_Max = 0.1;

    const char *filename = "../NACA_points.csv";  
    
    int num_lines_3 = count_lines(filename);  

    if (num_lines_3 <= 1) {
        printf("Erreur : Pas assez de points dans le fichier.\n");
        return 1;  
    }
    printf("Nombre de lignes : %d\n", num_lines_3);

    // Récupérer le tableau de points transformés
    double (*points_NACA)[2] = transform_NACA(filename, num_lines_3);

    if (!points_NACA) {
        printf("Erreur lors de la transformation de Joukovski.\n");
        return 1;
    }

    geoMeshGenerate(filename, num_lines_3);
    geoMeshImport();

    clock_t begin = clock();  // Start timing here

    // --- Renumérotation RCM ---
    printf("Starting RCM node renumbering...\n");
    femGraph* meshGraph = femBuildAdjacencyGraph(theGeometry->theElements, theGeometry->theEdges, theGeometry->theNodes->nNodes);
    if (meshGraph) {
        int startNode = femFindMinDegreeNode(meshGraph);
        int* rcmPermutationP = femComputeRcmPermutation(meshGraph, startNode);
        if (rcmPermutationP) {
            int* old_to_new = femComputeOldToNewMap(rcmPermutationP, theGeometry->theNodes->nNodes);
            if (old_to_new) {
                femApplyNodePermutation(theGeometry, rcmPermutationP, old_to_new);
                free(old_to_new);
                printf("RCM renumbering applied successfully.\n");
            } else {
                printf("Error computing old_to_new map.\n");
            }
            free(rcmPermutationP);
        } else {
            printf("RCM permutation computation failed (check graph connectivity?).\n");
        }
        femFreeGraph(meshGraph); // Libérer le graphe après usage
    } else {
        printf("Failed to build adjacency graph.\n");
    }
    // --- Fin Renumérotation RCM ---

    free(points_NACA);
    geoSetDomainName(0, "cercle1"); 
    geoSetDomainName(1,"cercle2");
    geoSetDomainName(2,"cercle3");

    for (int i = 0; i < theGeometry->nDomains; i++) {
        printf("Domain %d: %s\n", i, theGeometry->theDomains[i]->name);
    }

    geoMeshWrite("../data/elasticity.txt");
    
    // --- Création du problème ---
    double E = 60.e9; 
    double nu = 0.33; 
    double rho = 2.7e3;
    double g = 9.81;
    femProblem* theProblem = femElasticityCreate(theGeometry, E, nu, rho, g, PLANAR_STRESS);
    
    femElasticityAddBoundaryCondition(theProblem,"cercle1", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem,"cercle1", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem,"cercle2", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem,"cercle2", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem,"cercle3", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem,"cercle3", DIRICHLET_Y, 0.0);

    // Tableaux de pressions pour l'intrados et l'extrados
    double pressures_intrados[] = {500.0, 1000.0, 1500.0};   // Pressions pour l'intrados
    double pressures_extrados[] = {-2000.0, -3000.0, -4000.0}; // Pressions pour l'extrados

    // Nombre de pressions dans chaque tableau
    int numPressuresIntrados = sizeof(pressures_intrados) / sizeof(pressures_intrados[0]);
    int numPressuresExtrados = sizeof(pressures_extrados) / sizeof(pressures_extrados[0]);

    // Index pour naviguer dans les tableaux
    int currentPressureIntradosIndex = 1;
    int currentPressureExtradosIndex = 2;


    // Appel de la fonction pour changer et appliquer les pressions
    
    for (int i = 0; i < theGeometry->nDomains; i++) {
        if(i%2 == 0){
            femElasticityAddBoundaryCondition(theProblem,theGeometry->theDomains[i]->name,NEUMANN_Y,pressures_intrados[currentPressureIntradosIndex]);//intrados
        }
        if(i%2 !=0){
            femElasticityAddBoundaryCondition(theProblem,theGeometry->theDomains[i]->name,NEUMANN_Y,pressures_extrados[currentPressureExtradosIndex]);//extrados
        }
    }
    //applyPressureChanges(theProblem, theGeometry, pressures_intrados, pressures_extrados, numPressuresIntrados, numPressuresExtrados, &currentPressureIntradosIndex, &currentPressureExtradosIndex);
    femElasticityPrint(theProblem);

    // Résolution du problème et calcul des forces
    double *theSoluce = femElasticitySolve(theProblem);
    double *theForces = femElasticityForces(theProblem);
    double area = femElasticityIntegrate(theProblem, fun);   

    clock_t end = clock();    // End timing here
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Time taken: %f seconds\n", time_spent);
   
    // Deformation du maillage pour le plot final
    femNodes *theNodes = theGeometry->theNodes;
    double deformationFactor = 1e5;
    double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
    double *forcesX = malloc(theNodes->nNodes * sizeof(double));
    double *forcesY = malloc(theNodes->nNodes * sizeof(double));
    
    for (int i = 0; i < theNodes->nNodes; i++) {
        theNodes->X[i] += theSoluce[2*i+0]*deformationFactor;
        theNodes->Y[i] += theSoluce[2*i+1]*deformationFactor;
        normDisplacement[i] = sqrt(theSoluce[2*i+0]*theSoluce[2*i+0] + 
                                   theSoluce[2*i+1]*theSoluce[2*i+1]);
        forcesX[i] = theForces[2*i+0];
        forcesY[i] = theForces[2*i+1];
    }

    // Calcul de la force globale résultante
    double theGlobalForce[2] = {0, 0};
    for (int i = 0; i < theProblem->geometry->theNodes->nNodes; i++) {
        theGlobalForce[0] += theForces[2*i+0];
        theGlobalForce[1] += theForces[2*i+1];
    }
    printf(" ==== Global horizontal force       : %14.7e [N] \n", theGlobalForce[0]);
    printf(" ==== Global vertical force         : %14.7e [N] \n", theGlobalForce[1]);
    printf(" ==== Weight                        : %14.7e [N] \n", area * rho * g);

    // Visualisation du maillage
    int mode = 1; 
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[MAXNAME];

     // Code d'animation
    
    if (enableAnimation) {
   
    GLFWwindow* window = glfemInit("EPL1110 : Recovering forces on constrained nodes");
    glfwMakeContextCurrent(window);

    // Boucle d'affichage
    do {
        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glfemReshapeWindows(theGeometry->theNodes, w, h);
        t = glfwGetTime();  
        
        // Gestion des différentes actions dans la fenêtre
        if (glfwGetKey(window, 'D') == GLFW_PRESS) { mode = 0; }
        if (glfwGetKey(window, 'V') == GLFW_PRESS) { mode = 1; }
        if (glfwGetKey(window, 'X') == GLFW_PRESS) { mode = 2; }
        if (glfwGetKey(window, 'Y') == GLFW_PRESS) { mode = 3; }
        if (glfwGetKey(window, 'N') == GLFW_PRESS && freezingButton == FALSE) { 
            domain++; 
            freezingButton = TRUE; 
            told = t;
        }
        if (t - told > 0.5) { freezingButton = FALSE; }

        if (mode == 0) {
            domain = domain % theGeometry->nDomains;
            glfemPlotDomain(theGeometry->theDomains[domain]); 
            sprintf(theMessage, "%s : %d ", theGeometry->theDomains[domain]->name, domain);
            glColor3f(1.0, 0.0, 0.0); glfemMessage(theMessage); 
        }
        if (mode == 1) {
            glfemPlotField(theGeometry->theElements, normDisplacement);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ", theGeometry->theElements->nElem);
            glColor3f(1.0, 0.0, 0.0); glfemMessage(theMessage); 
        }
        if (mode == 2) {
            glfemPlotField(theGeometry->theElements, forcesX);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ", theGeometry->theElements->nElem);
            glColor3f(1.0, 0.0, 0.0); glfemMessage(theMessage); 
        }
        if (mode == 3) {
            glfemPlotField(theGeometry->theElements, forcesY);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ", theGeometry->theElements->nElem);
            glColor3f(1.0, 0.0, 0.0); glfemMessage(theMessage); 
        }

        glfwSwapBuffers(window);
        glfwPollEvents();
    } while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(window) != 1);
    glfwTerminate();
    }
    

    free(normDisplacement);
    free(forcesX);
    free(forcesY);
    femElasticityFree(theProblem);
    geoFinalize();
    glfwTerminate(); 
    
    exit(EXIT_SUCCESS);
    return 0;  
}
*/