
###################
#MICKAEL 19/10/2024
###################

#On va élaborer un dictionnaire de correspondances des flags sam établir une correspondance pour pouvoir les renvoyers dans
#un fichier de sortie :

# Dans les notations : Read : Lecture single-end ou pair-end
#                      Segment : un fragment d'ADN sens ou antisens
# Dictionnaires de correspondances bits : Flags SAM

Binary_sam = { 1: "Read apparié.",
             2: "Segment apparié correctement selon les critères de l'aligneur.",
             4: "Segment particulier non aligné.",
             8: "Segment complémentaire aligné sur le brin négatif.",
            16: "Segment aligné sur le brin négatif.",
            32: "L'autre read est aligné en réverse sur le brin positif.",
            64: "Il s'agit du premier read d'une paire sur le brin positif (5'->3').",
           128: "Il s'agit du second read d'une paire sur le brin negatif (5' -> 3').",
           256: "Alignement secondaire (non spécifique, alignement multiple).",
           512: "Read qui n'a pas passé les filtres de qualité.",
          1024: "Duplication du à la PCR ou au processus optique.",
          2048: "Alignement supplémentaire (non spécifique, alignement multiple).", }


# nos commentaires concaténés (du dictionnaire) en fonction de la valeur de son flag,
# Retourne une liste de commentaires des flags présents dans le dictionnaire Binary_sam,
# en examinant chaque paire clé-valeur (bit, commentaire) on inclut le commentaire dans la liste si
# la valeur_flag à le bit correspondant activté (sur 1).
# On utilise l'opérateur & (bitwise AND) pour vérifier si le bit est activé.

def Decodage_flags(valeur_du_flag):
    # Parcourt les flags et imprime chaque commentaire avec une puce
    for bit, commentaire in Binary_sam.items():
        if valeur_du_flag & bit:
            print(f"- {commentaire}")  # Affiche chaque commentaire avec une puce

test = 99
print(Decodage_flags(test))


























































































²
