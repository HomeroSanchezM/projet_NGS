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
4. [How to use the Script?](#how-to-use-the-script)
    1. [Complete execution](#1-Complete-execution)
    2. [Execution with options](#2-Execution-with-options)
5. [Output](#structure-des-fichiers)
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
- **pandas**
- **pars**
- **os**
- **re**
  
In order to automatically download the packages we have included a requirement.txt file, in order to run it you must enter the following line of code in your terminal:  
```
pip install -r requirements.txt 
```
### Type of the input 

The project is adapted for the analisis of SAM files containing aligned DNA sequence reads. During the execution, the SAM file is the only Imput needed and have to be add to the same folder that contain the project's files. For more information ot the SAM format, please read the documentation **descriptive_sam.pdf** of the project. 

## How to use the script? 

### Complete execution

Here you have an example of execution using a SAM file present in the project files that will execute the complete analysis of the SAM, you need to be in the same folder that the project files and enter the following line of code in your terminal:

```
./main.sh mapping.sam
```
this execution will analyse and print by default a complete analysis of the your SAM file, if you're only interested in doing a partial analysis, you can options.  

### execution with options

The program have differents options that allow the user to only extract the information needed, the options are:

- **--all** : analysis of all 
- **--cigar** : analysis of the cigar information 
- **--base** : percentage of different bases
- **--flag** : analysis of the flags
- **--pos** : analysis of the position of the fisrt AA
- **--ali** : Sequence of the alignement
- **--qual** : quality of the mapping

the list of options is avalable with the folowing line of code: 

```
./main.sh mapping.sam --help
```
## Output 


The execution of the program will print in the terminal tables that give a summary vieuw of the analysis.


Licence GPLv3

