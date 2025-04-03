#include "fem.h"

//initialisation de la matrice A et du vecteur B
double **A_copy = NULL;
double *B_copy  = NULL;
char* file = "homework.c";

void femElasticityAssembleElements(femProblem *theProblem){
    femFullSystem  *theSystem = theProblem->system;
                     //femBandSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    femMesh        *theEdges = theGeometry->theEdges;
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    int nLocal = theMesh->nLocalNode;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];} 
        
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }            
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + 
                                            dphidy[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + 
                                            dphidy[i] * c * dphidx[j]) * jac * weight;                                                                                           
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + 
                                            dphidx[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + 
                                            dphidx[i] * c * dphidx[j]) * jac * weight; }}
             for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * g * rho * jac * weight; }}} 
}


void femElasticityAssembleNeumann(femProblem *theProblem){
    femFullSystem  *theSystem = theProblem->system;

    // femBandSystem  *theSystem = theProblem->system;

    femIntegration *theRule = theProblem->ruleEdge;
    femDiscrete    *theSpace = theProblem->spaceEdge;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theEdges = theGeometry->theEdges;
    double x[2],y[2],phi[2];
    int iBnd,iElem,iInteg,iEdge,i,j,d,map[2],mapU[2];
    int nLocal = 2;
    double *B  = theSystem->B;


    //on commence par parcourir toutes les conditions définies aux limites
    for(iBnd=0; iBnd < theProblem->nBoundaryConditions; iBnd++){
        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
        femBoundaryType type = theCondition->type;
        femDomain *domain = theCondition->domain;
        double value = theCondition->value;

        // Ignorer les conditions de Dirichlet (contraintes de déplacement)
        if (type == DIRICHLET_X || type == DIRICHLET_Y) { continue; }
        
        int shift = (type == NEUMANN_X) ? 0 : 1;

        //on parcourt toutes les arêtes du maillage
        for(iEdge = 0; iEdge < domain->nElem; iEdge++){

            //on prend l'indice de l'élément associé à l'arête
            iElem = domain->elem[iEdge];

            //on récupère les coordonnées des noeuds de l'élément
            for(j=0; j<nLocal; j++){
                map[j]=theEdges->elem[iElem * nLocal + j];
                mapU[j] = 2 * map[j] + shift;
                x[j] = theNodes->X[map[j]];
                y[j] = theNodes->Y[map[j]];
            }


            double dx = x[1] - x[0];
            double dy = y[1] - y[0];
            double length = sqrt(dx*dx + dy*dy);    
            double jac = length/2.0; //jacobien d'une arrête en 2d

            //intégration sur l'arête (The Rule)
            for(iInteg = 0;iInteg < theRule->n; iInteg++){

                //coos du point d'intégration
                double xsi = theRule->xsi[iInteg];
                double weight = theRule->weight[iInteg];

                femDiscretePhi(theSpace,xsi,phi);

                for (int i = 0; i<theSpace->n; i++){
                    B[mapU[i]] += phi[i] * value * jac * weight;
                }
            }
        }
    }
}

void femBandSystemSet(femBandSystem* myBandSystem, int myRow, int myCol, double value) {
    // Vérifie si l'élément est dans la bande supérieure stockée
    // (band est le nombre total de diagonales = p + 1)
    if (myCol >= myRow && myCol < myRow + myBandSystem->band) {
        // Vérifie les limites globales (sécurité additionnelle)
        if (myRow >= 0 && myRow < myBandSystem->size && myCol >= 0 && myCol < myBandSystem->size) {
             // Accès direct A[row][col] (basé sur fem.c) - RISQUÉ
             myBandSystem->A[myRow][myCol] = value;
        } else {
             // Optionnel: Gérer l'erreur si les indices sont hors limites globales
             femError("Indices hors limites dans femBandSystemSet", __LINE__, file);
        }
    } else {
         // Optionnel: Gérer l'erreur si on essaie d'écrire hors de la bande stockée
         // Si la valeur qu'on essaie d'écrire n'est pas nulle, c'est peut-être une erreur.
         if (fabs(value) > 1e-16) {
              // Warning ou Error si on tente d'écrire une valeur significative hors bande
              // femWarning("Tentative d'écriture hors bande dans femBandSystemSet", __LINE__, file);
         }
    }
}


femBandSystem *theBandSystem_copy = NULL;
double *femElasticitySolve(femProblem *theProblem){

    femGeo *theGeometry = theProblem->geometry;
    femMesh *theMesh = theGeometry->theElements;
    femNodes *theNodes = theGeometry->theNodes;
    int nTotalNodes = theNodes->nNodes;
    int size = nTotalNodes * 2;
    int spaceDim = 2; // 2D Elasticity

    printf("Using BAND SOLVER workflow\n");

    // ===== ÉTAPE 1 : Calculer la largeur de bande MANUELLEMENT =====
    printf("Computing Bandwidth Manually...\n");
    int maxDifference = 0;

    for (int iElem = 0; iElem < theMesh->nElem; iElem++) {
        int nLoc = theMesh->nLocalNode;
        int minDOF = size;
        int maxDOF = -1;
        for (int j = 0; j < nLoc; j++) {
            int nodeIndex = theMesh->elem[iElem * nLoc + j];
            for (int k = 0; k < spaceDim; k++) {
                int globalDOF = nodeIndex * spaceDim + k;
                if (globalDOF < minDOF) minDOF = globalDOF;
                if (globalDOF > maxDOF) maxDOF = globalDOF;
            }
        }
        int diff = maxDOF - minDOF;
        if (diff > maxDifference) maxDifference = diff;
    }

    int semiBandWidth_p = maxDifference;
    // D'après femBandSystemGet, 'band' est p+1 (nb total de diagonales stockées)
    int bandWidth_for_alloc = semiBandWidth_p + 1;

    if (bandWidth_for_alloc <= 0) {
        femError("Error computing band width manually (result <= 0)", __LINE__, file);
        return NULL;
    }
    printf("  Manually Computed Semi-Bandwidth p = %d\n", semiBandWidth_p);
    printf("  Using bandWidth = %d for allocation (p+1)\n", bandWidth_for_alloc);
    // =======================================================================

    // ===== ÉTAPE 2 : Allouer le système BANDE =====
    printf("Allocating Band System (size %d, band %d)\n", size, bandWidth_for_alloc);
    femBandSystem *theBandSystem = femBandSystemCreate(size, bandWidth_for_alloc);
    if (theBandSystem == NULL) {
        femError("Could not allocate band system", __LINE__, file);
        return NULL;
    }
    // =====================================================

    // ===== ÉTAPE 3 : Initialiser le système BANDE =====
    femBandSystemInit(theBandSystem); // Met A et B à zéro
    // ================================================

    // ===== ÉTAPE 4 : Intégrer l'assemblage des éléments (BANDE) =====
    printf(" ----- Assembling elements (Band) -----\n");
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    // theMesh, theNodes sont déjà définis
    int nLoc = theSpace->n;
    int nElem = theMesh->nElem;
    double x[nLoc], y[nLoc], phi[nLoc], dphidxsi[nLoc], dphideta[nLoc], dphidx[nLoc], dphidy[nLoc];
    int map[nLoc], mapX[nLoc], mapY[nLoc]; // mapX/mapY indices globaux pour DDLs
    double a = theProblem->A;
    double b = theProblem->B;
    double c = theProblem->C;
    double rho = theProblem->rho;
    double g = theProblem->g;
    int nLocal = theSpace->n;

    for (int iElem = 0; iElem < nElem; iElem++) {
        // Récupérer infos élément (map, x, y)
        for (int j = 0; j < nLoc; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];     // Indice global DDL X du noeud local j
            mapY[j] = 2*map[j] + 1; // Indice global DDL Y du noeud local j
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];
        }

        // Boucle sur les points d'intégration
        for (int iInteg = 0; iInteg < theRule->n; iInteg++) {
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];
            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

            // Calcul Jacobien et dérivées globales
            double dxdxsi = 0.0, dxdeta = 0.0, dydxsi = 0.0, dydeta = 0.0;
            for (int i = 0; i < nLocal; i++) {
                dxdxsi += x[i]*dphidxsi[i]; dxdeta += x[i]*dphideta[i];
                dydxsi += y[i]*dphidxsi[i]; dydeta += y[i]*dphideta[i];
            }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            if (jac < 1e-14) femError("Jacobian null",__LINE__,file);

            for (int i = 0; i < nLocal; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }

            // Assemblage dans A (matrice bande) et B
            for (int i = 0; i < nLocal; i++) {
                for (int j = 0; j < nLocal; j++) {
                    double Aij_xx = (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * jac * weight;
                    double Aij_xy = (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * jac * weight;
                    double Aij_yx = (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * jac * weight; // Terme pour symétrie si nécessaire
                    double Aij_yy = (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * jac * weight;

                    // Ajout à la matrice bande (lecture puis écriture)
                    double oldVal_xx = femBandSystemGet(theBandSystem, mapX[i], mapX[j]);
                    double oldVal_xy = femBandSystemGet(theBandSystem, mapX[i], mapY[j]);
                    double oldVal_yx = femBandSystemGet(theBandSystem, mapY[i], mapX[j]);
                    double oldVal_yy = femBandSystemGet(theBandSystem, mapY[i], mapY[j]);

                    femBandSystemSet(theBandSystem, mapX[i], mapX[j], oldVal_xx + Aij_xx);
                    femBandSystemSet(theBandSystem, mapX[i], mapY[j], oldVal_xy + Aij_xy);
                    // Si A est symétrique, Aij_yx = Aji_xy. Si on stocke que la bande sup,
                    // on n'assemble pas Aij_yx explicitement SI on est sûr que le solveur
                    // gère la symétrie. Sinon, il faut l'assembler :
                    femBandSystemSet(theBandSystem, mapY[i], mapX[j], oldVal_yx + Aij_yx);
                    femBandSystemSet(theBandSystem, mapY[i], mapY[j], oldVal_yy + Aij_yy);
                }
                // Assemblage dans B (forces volumiques : gravité)
                theBandSystem->B[mapY[i]] -= phi[i] * g * rho * jac * weight;
            }
        } // Fin boucle intégration
    } // Fin boucle éléments
    // ====================================================================

    // ===== ÉTAPE 5 : Intégrer l'assemblage Neumann (BANDE) =====
    printf(" ----- Assembling Neumann contributions (Band) -----\n");
    femIntegration *theRuleEdge = theProblem->ruleEdge;
    femDiscrete    *theSpaceEdge = theProblem->spaceEdge;
    femMesh        *theEdges = theProblem->geometry->theEdges;
    int nLocalEdge = theSpaceEdge->n; // Devrait être 2

    for (int iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++) {
        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
        femBoundaryType type = theCondition->type;
        femDomain *domain = theCondition->domain;
        double value = theCondition->value;

        if (type == DIRICHLET_X || type == DIRICHLET_Y) continue;

        int shift = (type == NEUMANN_X) ? 0 : 1;

        for (int iEdge = 0; iEdge < domain->nElem; iEdge++) {
            int iElemEdge = domain->elem[iEdge];
            int mapEdge[nLocalEdge], mapU[nLocalEdge];
            double xEdge[nLocalEdge], yEdge[nLocalEdge], phiEdge[nLocalEdge];
            for (int j = 0; j < nLocalEdge; j++) {
                 mapEdge[j] = theEdges->elem[iElemEdge * nLocalEdge + j];
                 mapU[j] = 2 * mapEdge[j] + shift;
                 xEdge[j] = theNodes->X[mapEdge[j]];
                 yEdge[j] = theNodes->Y[mapEdge[j]];
            }

            double dx = xEdge[1] - xEdge[0];
            double dy = yEdge[1] - yEdge[0];
            double length = sqrt(dx*dx + dy*dy);
            double jacEdge = length / 2.0;

            for (int iInteg = 0; iInteg < theRuleEdge->n; iInteg++) {
                 double xsiEdge = theRuleEdge->xsi[iInteg];
                 double weightEdge = theRuleEdge->weight[iInteg];
                 femDiscretePhi(theSpaceEdge, xsiEdge, phiEdge);

                 for (int i = 0; i < nLocalEdge; i++) {
                     theBandSystem->B[mapU[i]] += phiEdge[i] * value * jacEdge * weightEdge;
                 }
            }
        }
    }
    // ============================================================


    // Allouer la copie (en supposant Option B)
    printf("Copying original band system for forces calculation...\n");
    theBandSystem_copy = femBandSystemCreate(size, bandWidth_for_alloc);
    if (!theBandSystem_copy) { /* Erreur */ }

    // Copier TOUT le bloc mémoire (B + données de A)
    // Taille du bloc = size * (band + 1) où band = bandWidth_for_alloc
    memcpy(theBandSystem_copy->B, theBandSystem->B, sizeof(double) * size * (bandWidth_for_alloc + 1));



    // ===== ÉTAPE 6 : Application des CL Dirichlet (BANDE - Manuel) =====
    printf("Applying Boundary Conditions (Band - Manual)...\n");
    int size_system = theBandSystem->size;
    // band_system est déjà défini comme bandWidth_for_alloc ? Utilisons ce nom.
    int band_system_width = theBandSystem->band; // = bandWidth_for_alloc

    for (int i = 0; i < theProblem->nBoundaryConditions; i++) {
        femBoundaryCondition *theCondition = theProblem->conditions[i];
        femDomain *theDomain = theCondition->domain;
        femBoundaryType type = theCondition->type; // Ajouté pour vérifier le type

        // Ne traiter que les conditions de Dirichlet ici
        if (type != DIRICHLET_X && type != DIRICHLET_Y) continue;

        for (int j = 0; j < theDomain->nElem; j++) {
            int iEdge = theDomain->elem[j];
            int nEdgeNodes = theGeometry->theEdges->nLocalNode;
            for (int k = 0; k < nEdgeNodes; k++) {
                 int node = theGeometry->theEdges->elem[iEdge * nEdgeNodes + k];
                 int constrainedDOF = -1;
                 double value = theCondition->value;

                 if (theCondition->type == DIRICHLET_X) constrainedDOF = node * 2 + 0;
                 if (theCondition->type == DIRICHLET_Y) constrainedDOF = node * 2 + 1;

                 if (constrainedDOF != -1) {
                     if (constrainedDOF < 0 || constrainedDOF >= size_system) {
                          femError("Constrained DOF out of bounds",__LINE__,file); continue;
                     }

                     // Appliquer contrainte (logique de femFullSystemConstrain)
                     // 1. Ajuster B
                     for (int row = 0; row < size_system; ++row) {
                         double A_ik = femBandSystemGet(theBandSystem, row, constrainedDOF);
                         theBandSystem->B[row] -= value * A_ik;
                     }
                     // 2. Zéro Colonne
                     for (int row = 0; row < size_system; ++row) {
                          femBandSystemSet(theBandSystem, row, constrainedDOF, 0.0);
                     }
                     // 3. Zéro Ligne (partie sup.)
                      for (int col = constrainedDOF; col < size_system; ++col) {
                          // Vérifie si dans la bande avant de mettre à zéro
                          if (col < constrainedDOF + band_system_width) {
                              femBandSystemSet(theBandSystem, constrainedDOF, col, 0.0);
                          }
                     }
                     // 4. Diagonale = 1
                     femBandSystemSet(theBandSystem, constrainedDOF, constrainedDOF, 1.0);
                     // 5. Set B
                     theBandSystem->B[constrainedDOF] = value;
                 }
            }
        }
    }
    // =====================================================================

    // ===== ÉTAPE 7 : Résoudre le système BANDE =====
    // ATTENTION : CETTE FONCTION MANQUE DANS fem.c !
    // Tu dois l'obtenir ou l'implémenter.
    printf("Solving Band System...\n");
    femBandSystemSolve(theBandSystem);
    // femError("FONCTION MANQUANTE : femBandSystemSolve n'est pas définie dans fem.c !",__LINE__,file);
    // ==============================================


    // ===== ÉTAPE 8 : Préparer la Solution =====f
    printf("Copying solution...\n");
    double *solution = (double*) malloc(sizeof(double) * size_system);
     if (!solution) {
        femError("Failed to allocate memory for solution copy", __LINE__, file);
        femBandSystemFree(theBandSystem); // Libérer avant de quitter
        return NULL;
    }
    // Copie B (qui contient la solution après solve) dans solution
    memcpy(solution, theBandSystem->B, size_system * sizeof(double));
    // Copie aussi dans theProblem->soluce si nécessaire pour la suite (ex: forces)
    memcpy(theProblem->soluce, theBandSystem->B, size_system * sizeof(double));
    // ========================================

    // ===== ÉTAPE 9 : Libérer la mémoire (BANDE) =====
    printf("Freeing Band System...\n");
    femBandSystemFree(theBandSystem);
    // ============================================

    // ===== ÉTAPE 10 : Retourner la Solution =====
    printf("Band solver workflow finished.\n");
    // Retourne la copie de la solution (ou theProblem->soluce si modifié)
    return solution; // L'appelant devra free(solution)
    // Ou return theProblem->soluce; si c'est la convention
    // =========================================
}


double * femElasticityForces(femProblem *theProblem){

    double *soluce = theProblem->soluce;
    int size = 2 * theProblem->geometry->theNodes->nNodes; // Obtenir size

    // Vérifier si la copie existe (créée dans femElasticitySolve)
    if (theBandSystem_copy == NULL) {
        femError("Original band system copy (theBandSystem_copy) not available in femElasticityForces", __LINE__, file);
        // Retourner des zéros ou NULL ?
         if (theProblem->residuals == NULL) {
             theProblem->residuals = (double*) calloc(size, sizeof(double)); // Alloue et met à zéro
         } else {
             memset(theProblem->residuals, 0, size * sizeof(double)); // Met à zéro
         }
        return theProblem->residuals;
    }

    double *residuals = theProblem->residuals;
    if(residuals==NULL){
        residuals = (double *) calloc(size, sizeof(double)); // calloc met à zéro
        theProblem->residuals = residuals;
         if(residuals==NULL) {
             femError("Failed to allocate residuals in femElasticityForces", __LINE__, file);
             return NULL;
         }
    }

    // Récupérer A et B originaux depuis la copie
    double **A_orig = theBandSystem_copy->A; // Attention à l'accès via Get
    double *B_orig = theBandSystem_copy->B;
    int band_orig = theBandSystem_copy->band;


    printf("Calculating residuals R = A_orig * U - B_orig...\n");
    // Calculer R = A_original * U - B_original
    for (int i = 0; i < size; i++) {
        double AxU_i = 0.0;
        // Calculer la i-ème composante de A_original * U
        // Boucle sur les colonnes j dans la bande de la ligne i
        for (int j = i; j < size && j < i + band_orig; j++) {
            // Utiliser Get sur la COPIE pour lire A_original(i, j)
            AxU_i += femBandSystemGet(theBandSystem_copy, i, j) * soluce[j];
        }
         // Si la matrice n'est PAS symétrique, il faut aussi la partie sous-diagonale
         // (qui n'est pas stockée par Get/Set tel quel).
         // Si on suppose la symétrie A(j,i) = A(i,j)
         for (int j = i - (band_orig - 1) ; j < i; j++) {
            if (j >= 0) { // Assurer j n'est pas négatif
                // A(i,j) = A(j,i) par symétrie. A(j,i) est dans la bande sup de la ligne j.
                 AxU_i += femBandSystemGet(theBandSystem_copy, j, i) * soluce[j];
            } 
         }


        residuals[i] = AxU_i - B_orig[i];
    }

    // On pourrait libérer la copie ici si elle n'est plus utile
    // femBandSystemFree(theBandSystem_copy);
    // theBandSystem_copy = NULL;
    // Mais si on relance le solveur/forces, il faudrait la refaire.
    // La gestion mémoire des copies globales peut être délicate.

    return theProblem->residuals;
}