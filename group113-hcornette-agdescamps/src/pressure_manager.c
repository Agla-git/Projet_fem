#include "pressure_manager.h"
#include "fem.h"

void applyPressureChanges(femProblem* theProblem, femGeo* theGeometry, double* pressures_intrados, double* pressures_extrados, int numPressuresIntrados, int numPressuresExtrados, int* currentPressureIntradosIndex, int* currentPressureExtradosIndex) {
    double pressureChangeTime = 2.0;  // Temps avant de changer la pression
    double lastPressureChangeTime = 0.0;  // Dernière fois où la pression a été changée
    double t = glfwGetTime();  // Récupère le temps actuel

    if (t - lastPressureChangeTime >= pressureChangeTime) {
        // Change la pression après 'pressureChangeTime' secondes
        *currentPressureIntradosIndex = (*currentPressureIntradosIndex + 1) % numPressuresIntrados;
        *currentPressureExtradosIndex = (*currentPressureExtradosIndex + 1) % numPressuresExtrados;
        lastPressureChangeTime = t;  // Réinitialiser le temps
    }

    // Appliquer la nouvelle pression dans la condition de Neumann
    for (int i = 0; i < theGeometry->nDomains; i++) {
        if (i % 2 == 0) {  // Intrados
            femElasticityAddBoundaryCondition(theProblem, theGeometry->theDomains[i]->name, NEUMANN_Y, pressures_intrados[*currentPressureIntradosIndex]);
        } else {  // Extrados
            femElasticityAddBoundaryCondition(theProblem, theGeometry->theDomains[i]->name, NEUMANN_Y, pressures_extrados[*currentPressureExtradosIndex]);
        }
    }
}
