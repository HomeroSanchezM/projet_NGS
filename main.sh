#!/bin/bash

###########################################################################################################################################################################################################
#                                                                          COQUERELLE MICKAEL - HOMERO SANCHEZ 18.10.2024                                                                                #      
#			        		               --  SCORE : SAM Characterization and Observational Report for Evaluations  --                                                              
###########################################################################################################################################################################################################

# CONVENTIONS - NOMENCLATURE :
# on choisit dans nos regles de nommages d'identifier chaque variable avec comme première lettre son type de donnée pour stocker l'information le long du script et se répérer.
# bBooleen, dDictionnaire, sString, cCharact, iInteger, fFloat ...

#------------------------------------ #
# Contrôles de conformités de l'input #
#------------------------------------ #

Rep="" # J'initialise mon Rep à null dans l'hypothèse ou la conjonction de mes vérifications préliminaires renverrais false.

# Vérification des conditions en une seule instruction if
if [ -z "$1" ] || [ ! -f "$1" ] || [ -d "$1" ] || [ ! -s "$1" ]; then

    # Montage du texte d'erreur dans l'hypothèse ou au moins une de nos opérandes est fausse :
    [ -z "$1" ]   &&  Rep+="Erreur: Aucun fichier SAM spécifié.\n"
    [ ! -f "$1" ] &&  Rep+="Erreur: Le fichier n'existe pas.\n"
    [ -d "$1" ]   &&  Rep+="Erreur: Le fichier est un répertoire.\n"
    [ ! -s "$1" ] &&  Rep+="Erreur: Le fichier est vide.\n"
    
    # Affichage des messages d'erreurs
    echo -e "Erreur(s) détectée(s) :\n$Rep"
    exit 1
else
    # Si toutes les conditions sont satisfaites, on execute notre Script Python et ca bombarde ;) ...
    
    echo "Nos controles préliminaires sont terminés, élaboration du rapport en cours ..."
    
    python3 NGS.py "$1"
    
fi


