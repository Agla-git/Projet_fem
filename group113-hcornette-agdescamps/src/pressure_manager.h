#ifndef PRESSURE_MANAGER_H
#define PRESSURE_MANAGER_H

#include "glfem.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gmshc.h"
#include <complex.h>

// Déclaration de la fonction pour appliquer les changements de pression
void applyPressureChanges(femProblem* theProblem, femGeo* theGeometry, double* pressures_intrados, double* pressures_extrados, int numPressuresIntrados, int numPressuresExtrados, int* currentPressureIntradosIndex, int* currentPressureExtradosIndex);

#endif
