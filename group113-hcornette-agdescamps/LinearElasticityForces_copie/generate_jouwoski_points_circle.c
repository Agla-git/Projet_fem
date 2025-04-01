#include <stdio.h>
#include <math.h>

#define N 100   // Nombre de points
#define R 1.2  // Rayon du cercle

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void generate_csv(const char *filename) {
    FILE *file = fopen(filename, "w");  // Ouvrir le fichier en écriture
    if (!file) {
        printf("Erreur d'ouverture du fichier !\n");
        return;
    }

    fprintf(file, "Re,Im\n");  // Écriture de l'en-tête

    for (int i = 0; i < N; i++) {
        //double theta = -M_PI / 2 + M_PI * i / (N);  // Angles entre -π/2 et π/2 ==> je pense que c'est pour créer l'intrados.
        //double theta = M_PI / 2 - M_PI * i / (N);  // Angles entre π/2 et -π/2
        double theta = -M_PI + 2 * M_PI * i / (N);  // Angles entre -π et π
        double x = R * cos(theta)- 0.04;  // Décalage pour simuler la cambrure
        double y = R * sin(theta) + 0.04;  // Coordonnée y

        fprintf(file, "%f,%f\n", x, y);  // Écriture des points dans le fichier
    }

    
    // Générer l'extrados (de π/2 à -π/2, dans l'autre sens)
    // N-1 pour ne pas avoir 2 fois le même point (on l'a déjà dans l'intrados)
    /*
    for (int i = N; i >= 1; i--) {
        double theta = -M_PI / 2 + M_PI * i / N;
        double x = -R * cos(theta)+ 0.04;
        double y = R * sin(theta);  // On inverse y pour créer l'extrados
        fprintf(file, "%f,%f\n", x, y);
    }
    */

    fclose(file);
    printf("Fichier '%s' généré avec succès !\n", filename);
}

int main() {
    generate_csv("jouwoski_base_points.csv");  // Générer le fichier CSV
    return 0;
}
 