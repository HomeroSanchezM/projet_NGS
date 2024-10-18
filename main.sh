#!/bin/bash

###########################################################################################################################################################################################################
#                                                                          COQUERELLE MICKAEL - HOMERO SANCHEZ 18.10.2024                                                                                #      
#			        		               --  SCORE : SAM Characterization and Observational Report for Evaluations  --                                                              
###########################################################################################################################################################################################################

# CONVENTIONS - NOMENCLATURE :
# on choisit dans nos regles de nommages d'identifier chaque variable avec comme première lettre son type de donnée pour stocker l'information le long du script et se répérer.
# bBooleen, dDictionnaire, sString, cCharact, iInteger, fFloat ...

#-------------------------------------- #
# Contrôles des conformités sur l'input #
#-------------------------------------- #

Rep="" # J'initialise mon Rep à null dans l'hypothèse où la conjonction de mes vérifications renverrait false.

# Vérification des conditions en une seule instruction if
if [ -z "$1" ] || [ ! -f "$1" ] || [ -d "$1" ] || [ ! -s "$1" ]; then

    # Montage du texte d'erreur dans l'hypothèse où au moins une de nos opérandes est fausse :
    [ -z "$1" ]   &&  Rep+="Erreur: Aucun fichier SAM spécifié.\n"
    [ ! -f "$1" ] &&  Rep+="Erreur: Le fichier n'existe pas.\n"
    [ -d "$1" ]   &&  Rep+="Erreur: Le fichier est un répertoire.\n"
    [ ! -s "$1" ] &&  Rep+="Erreur: Le fichier est vide.\n"
    
    # Affichage des messages d'erreurs
    echo -e "Erreur(s) détectée(s) :\n$Rep"
    exit 1
else
    # Si toutes les conditions sont satisfaites, on exécute notre Script Python

    debNGSpy=$(date +%s) # Correction : espace ajouté après `date`

    echo "Nos contrôles préliminaires sont terminés, élaboration du rapport : "
    
    # On lance NGS.py sur notre premier paramètre
    python3 NGS.py "$1"
    
    endNGSpy=$(date +%s) # Correction : espace ajouté après `date`
    
    # Affichage du temps de traitement
    echo -e "\nTemps de traitement : $((endNGSpy - debNGSpy)) secondes"
    
fi
