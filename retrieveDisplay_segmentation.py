# -*- coding: utf-8 -*-
"""
Created on Thu Feb 05 11:07:47 2015

@ author (C) Cristina Gallego, University of Toronto
----------------------------------------------------------------------
 """
import os, os.path
import sys
from sys import argv, stderr, exit
import numpy as np
import dicom
import psycopg2
import pandas as pd

from sqlalchemy import Column, Integer, String
import datetime
from mybase import myengine
import mydatabase
from sqlalchemy.orm import sessionmaker
from sendNew2_mydatabase import *
from sendNew2_updatedatabase import *


if __name__ == '__main__':    
    # Get Root folder ( the directory of the script being run)
    path_rootFolder = os.path.dirname(os.path.abspath(__file__))
    print path_rootFolder
       
    # Open filename list
    file_masterlist = sys.argv[1]
    print "Processing masterlist: %s" % file_masterlist
    file_ids = open(file_masterlist,"r")
    file_ids.seek(0)
    line = file_ids.readline()
        
    while ( line ) : 
        # Get the line: 
        fileline = line.split()
        lesion_id = int(fileline[0] )
        BenignNMaligNAnt = fileline[1]
        SeriesID = fileline[2]
        T2SeriesID = fileline[3]
        StudyID = fileline[4]
        Lesionfile = fileline[5]
        DicomExamNumber = fileline[6]
        AccessionN = fileline[7]
        dateID = fileline[8]
        ismass = bool(fileline[9])
        isnonmass = bool(fileline[10])
        cond = fileline[15] 
        Diagnosis = fileline[16]     
        sideBreast = fileline[14]

        if(lesion_id>272):
            DicomExamNumber = AccessionN
            nameSegment = StudyID+'_'+AccessionN+'_'+str(lesion_id)+'.vtk'
        else:
            ## for old DicomExamNumber  
            nameSegment = StudyID+'_'+AccessionN+'_'+'lesion.vtk'
            
        img_folder ='Z:/Cristina/MassNonmass'+os.sep+cond[:-1]+'/'
        path_rootFolder= 'C:/Users/windows/Documents/repoCode-local/stage1features'
        pathSegment = 'C:\Users\windows\Documents\repoCode-local\stage1features\seg+T2'
 
        #############################
        ###### 1) Querying Research database for clinical, pathology, radiology data
        #############################
        SendNew2DB = SendNew()
        [img_folder, cond, BenignNMaligNAnt, Diagnosis, rowCase] = SendNew2DB.queryNewDatabase(StudyID, dateID, Diagnosis)        
        
        #############################
        ###### 2) Get T2 info
        #############################        
        Send2DB = SendNewUpdate()
        [path_T2Series, T2SeriesID] = Send2DB.processT2(T2SeriesID, img_folder, StudyID, DicomExamNumber)

        #############################                  
        ###### 3) Load segmentation and display
        #############################
        [series_path, phases_series, lesion3D, chgSeg] = Send2DB.checkSegment(path_rootFolder, cond, StudyID, DicomExamNumber, SeriesID, Lesionfile, T2SeriesID, path_T2Series, lesion_id, sideBreast)
        
        ## continue to next case
        line = file_ids.readline()
        print line
            
        
    file_ids.close()            

