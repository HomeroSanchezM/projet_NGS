#!/bin/bash

###########################################################################################################################################################################################################
#                                                                          COQUERELLE MICKAEL - HOMERO SANCHEZ 18.10.2024                                                                                #      
#			        		               --  SCORE : SAM Characterization and Observational Report for Evaluations  --                                                              
###########################################################################################################################################################################################################

sep="--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"

Script1="NGS.py"

#-------------------------------------- #
# Contrôles de conformité sur l'input   #
#-------------------------------------- #

Rep="" # J'initialise mon Rep à null dans l'hypothèse où la conjonction de mes vérifications renverrait false.

# Vérification des conditions en une seule instruction if
if [ -z "$1" ] || [ ! -f "$1" ] || [ -d "$1" ] || [ ! -s "$1" ]; then
    [ -z "$1" ]   &&  Rep+="Erreur: Aucun fichier SAM spécifié.\n"
    [ ! -f "$1" ] &&  Rep+="Erreur: Le fichier n'existe pas.\n"
    [ -d "$1" ]   &&  Rep+="Erreur: Le fichier est un répertoire.\n"
    [ ! -s "$1" ] &&  Rep+="Erreur: Le fichier est vide.\n"
    
    echo -e "Erreur(s) détectée(s) :\n$Rep"
    exit 1
else
    # Si toutes les conditions sont satisfaites, on exécute notre Script Python
    debNGSpy=$(date +%s) 
   
    echo -e "\nNos contrôles préliminaires sont terminés, début du traitement : \n"

# ------------------------------------------------------------------------------------------------------------------------------------------------------#
# Création d'un Output format .csv avec une dénomination qui va gérer les éxecutions successives du programme pour éviter l'écrasement de l'output n-1  #
# ----------------------------------------------------------------------------------------------------------------------------------------------------- #
# Définir le nom du fichier de sortie avec la date :
#Output_CSV="Output_DATA_$(date +"%Y%m%d").csv"

# Créer le fichier avec la première ligne d'en-tête
#echo "Fichier d'extraction des métriques SAM : " > "$Output_CSV"

#-------------------------------------- #
# Gestion des options en arguments      #
#-------------------------------------- #
  # Extraction des options supplémentaires
    options="${@:2}" # Tous les arguments après le premier ($1 est le fichier SAM)

    # Gestion des options si on ne veut pas lancer l'intégralité de l'analyse : 
    if [ -z "$options" ]; then
        echo "Aucune option spécifique fournie, toutes les analyses seront effectuées par défaut."
        python3 "$Script1" "$1" --all
    else
        echo "Options spécifiées : $options"
        python3 "$Script1" "$1" $options
    fi
    

    if [ $? -ne 0 ]; then 
        echo "Erreur lors de l'extraction des flags et du comptage des occurrences, vérifier le script_1 ou si le nom du fichier n'est pas modifié."
        exit 1
    else 
        echo -e "\nExtractions des Flags et comptages des occurrences effectués avec succès\n"
    fi
    
    echo -e "\n$sep"
    # Affichage du temps de traitement
    endNGSpy=$(date +%s) 
    echo -e "\nTemps de traitement : $((endNGSpy - debNGSpy)) secondes"
   
fi

