# NGS Projet : SAM FILE EVALUATION PROGRAM
**COQUERELLE MICKAEL - HOMERO SANCHEZ**  
**Date** : 17.11.2024  
**SCORE** : SAM Characterization and Observational Report for Evaluations

---
 
## Table of content

1. [Introduction](#introduction)
2. [Prerequisites](#Prerequisites)
3. [Installation and configuration](#installation-and-configuration)
    1. [Import of the project](#1-Import-of-the-project)
    2. [installation of outbuildings with `requirements.txt`](#2-Instalation-of-outbuildings-with-requirement.txt)
    3. [Type of the input](#2-Type-of-the-imput)
4. [Script execution](#Script-execution)
    1. [Fichier `main.sh`](#1-fichier-mainsh)
    2. [Fichier `NGS.py`](#2-fichier-ngspy)
5. [Structure des fichiers](#structure-des-fichiers)
6. [Auteurs](#auteurs)

---

## Introduction

This project aims to analyze and characterize **SAM (Sequence Alignment/Map)** files. It allows you to extract, for example, detailed information on the flags present in a SAM file, to count their occurrence and to produce a report in the form of a table, observing the different characteristics of the data.
---

## Prerequisites

Before you can run this project, you need to have certain tools and programs installed on your machine.

- **Python 3.6** : Used to execute scripts.
- **Pip** : Used to manage Python dependencies.
- **Linux or similar** : This project is designed to run on Linux.
---

## Installation et configuration

### Import of the project 
To install all the files in our project, you need to enter the following line of code in your terminal: 
```
git clone https://github.com/HomeroSanchezM/projet_NGS.git 
```
### Installation of outbuildings  
The program uses the following packages:
- **panda**
- **pars**
- **os**
- **re**
  
In order to automatically download the packages we have included a requirement.txt file, in order to run it you must enter the following line of code in your terminal:  
```
pip install -r requirements.txt 
```
### Type of the input 

The project is adapted for the analisis of SAM files containing aligned DNA sequence reads. During the execution, the SAM file is the only Imput needed and have to be add to the same folder that contain the project's files. For more information ot the SAM format, please read the documentation **descriptive_sam.pdf** of the project. 

## Script execution



