###################
# MICKAEL 19/10/2024
###################

# Dictionnaire qui stocke des bits comme clés et des commentaires comme valeurs
d_Binary_sam = {
    1: "Read apparié.",
    2: "Segment apparié correctement selon les critères de l'aligneur.",
    4: "Segment particulier non aligné.",
    8: "Segment complémentaire aligné sur le brin négatif.",
    16: "Segment aligné sur le brin négatif.",
    32: "L'autre read est aligné en réverse sur le brin positif.",
    64: "Il s'agit du premier read d'une paire sur le brin positif (5'->3').",
    128: "Il s'agit du second read d'une paire sur le brin négatif (5' -> 3').",
    256: "Alignement secondaire (non spécifique, alignement multiple).",
    512: "Read qui n'a pas passé les filtres de qualité.",
    1024: "Duplication due à la PCR ou au processus optique.",
    2048: "Alignement supplémentaire (non spécifique, alignement multiple).",
}

# Fonction pour décoder la valeur d'un flag en affichant les commentaires correspondants
def decodage_flags(valeur_du_flag):
    l_synthese = []

    # Chaque bit représente un flag spécifique et son commentaire associé
    for i_bit, s_commentaire in d_Binary_sam.items():
        if valeur_du_flag & i_bit: # En gros ici on vérification l'activation ou non de chaque bit pour la valeur du flag
            l_synthese.append(f"- {s_commentaire}")

    return "\n".join(l_synthese) # On concaténe avec un RC les commentaires activés.

i_test = 98
print(decodage_flags(i_test))

