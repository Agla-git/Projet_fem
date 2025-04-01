import pandas as pd
import matplotlib.pyplot as plt


# Charger le fichier CSV contenant les points
data = pd.read_csv('jouwoski_base_points.csv')

# Extraire les coordonnées
x = data['Re']
y = data['Im']

# Créer le graphique
plt.figure(figsize=(16, 8))
plt.plot(x, y, 'o', label='Points du profil', color='blue')

# Ajouter un cercle de rayon 1 pour référence
circle = plt.Circle((0, 0), 1, color='red', fill=False, linestyle='--')
plt.gca().add_artist(circle)

# Configurer les limites des axes pour inclure tout le cercle
# Vous pouvez ajuster ces valeurs pour avoir plus d'espace autour du cercle si nécessaire
plt.xlim(-1.5, 1.5)  # Plage de l'axe X de -1.5 à 1.5
plt.ylim(-1.5, 1.5)  # Plage de l'axe Y de -1.5 à 1.5

# Configurer le graphique
plt.title('Profil de la géométrie')
plt.xlabel('Axe Re')
plt.ylabel('Axe Im')
plt.axis('equal')  # Assurer que l'échelle soit la même sur les deux axes
plt.grid(True)
plt.legend()

# Afficher le graphique
plt.show()
