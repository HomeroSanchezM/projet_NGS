#!/bin/bash

###########################################################################################################################################################################################################
#                                                                          COQUERELLE MICKAEL - HOMERO SANCHEZ 18.10.2024                                                                                #      
#			        		               --  SCORE : SAM Characterization and Observational Report for Evaluations  --                                                              
###########################################################################################################################################################################################################
sep="--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"

Script1="NGS.py"

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

    debNGSpy=$(date +%s) 
   
    echo -e "\nNos contrôles préliminaires sont terminés, création du rapport de sortie au format csv ... \n"

#-------------------------------------- #
# Génération de l'output au format .csv #
#-------------------------------------- #
  
    # Créer le fichier output.csv avec la première ligne d'en-tête
    echo "SCORE -- OUTPUT" > output.csv
    echo "FLAG_REPERE,NOMBRE_OCCURENCES,TRANSLATE_FLAG" >> output.csv
    echo "$sep"
    
    # On lance NGS.py  pour repérer les FLAGS présents et compter leurs occurences.
    python3 "$Script1" "$1"
    
    if [ $? -ne 0 ]; then #Si le code de retour de la dernière commande n'est pas égale à 0 
    echo "Erreur lors de l'extraction des flags et du comptage des occurences, vérifier le script_1 ou si le nom du fichier n'est pas modifié."
    exit 1
    else 
    echo -e "\nExtractions des Flags et comptages des occurences effectués avec succès\n"
    fi
    
    echo -e "\n$sep"
    # Affichage du temps de traitement
    endNGSpy=$(date +%s) 
    echo -e "\nTemps de traitement : $((endNGSpy - debNGSpy)) secondes"
    
fi

