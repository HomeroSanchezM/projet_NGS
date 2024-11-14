# Projet NGS : PROGRAMME D'ÉVALUATION DE FICHIER SAM
**COQUERELLE MICKAEL - HOMERO SANCHEZ**  
**Date** : 17.11.2024  
**SCORE** : SAM Characterization and Observational Report for Evaluations

---

## Table des matières

1. [Introduction](#introduction)
2. [Prérequis](#prérequis)
3. [Installation et configuration](#installation-et-configuration)
    1. [Créer un environnement virtuel](#1-créer-un-environnement-virtuel)
    2. [Installer Python](#2-installer-python)
    3. [Installer les dépendances avec `requirements.txt`](#3-installer-les-dépendances-avec-requirementstxt)
4. [Exécution des scripts](#exécution-des-scripts)
    1. [Fichier `main.sh`](#1-fichier-mainsh)
    2. [Fichier `NGS.py`](#2-fichier-ngspy)
    3. [Fichier `Translate_FLAG.py`](#3-fichier-translate_flagpy)
5. [Structure des fichiers](#structure-des-fichiers)
6. [Auteurs](#auteurs)

---

## Introduction

Ce projet a pour but d'analyser et de caractériser les fichiers **SAM (Sequence Alignment/Map)**. Il permet d'extraire par exemple des informations détaillées sur les flags présents dans un fichier SAM, de compter leur occurrence et de produire un rapport observant les différentes caractéristiques des données.

---

## Prérequis

Avant de pouvoir exécuter ce projet, vous devez disposer de certains outils et programmes installés sur votre machine.

- **Python 3.x** : Utilisé pour exécuter les scripts.
- **Pip** : Pour gérer les dépendances Python.
- **Système Linux ou similaire** : Ce projet est conçu pour fonctionner sous Linux.

---

## Installation et configuration

### 1. Créer un environnement virtuel

L'utilisation d'un environnement virtuel est essentielle pour isoler les dépendances spécifiques à ce projet. Cela vous permet de ne pas interférer avec d'autres projets ou les versions des bibliothèques installées sur votre système.

#### 1.1. Créer un environnement virtuel

Avant de créer un environnement virtuel, assurez-vous d'être dans le répertoire de votre projet. Si vous n'êtes pas dans le bon répertoire, utilisez la commande `cd` pour vous y rendre.

```bash
cd ~/projet

python3 -m venv SCORE_PROGRAMME

