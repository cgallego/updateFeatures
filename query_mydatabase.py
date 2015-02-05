# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 16:48:07 2014

@ author (C) Cristina Gallego, University of Toronto
"""
import sys, os
import string
import datetime
from numpy import *
import datetime
from mylocalbase import Base, myengine
import mylocaldatabase

from sqlalchemy.orm import sessionmaker

#!/usr/bin/env python
class QuerymyDatabase(object):
    """
    USAGE:
    =============
    localdata = QuerymyDatabase()
    """        
        
    def __init__(self): 
        """ initialize QueryDatabase """
        
                              
    def queryby_lesionid(self, lesion_id):
        """
        run : Query by Lesion_id on local database
        
        Inputs
        ======
        lesion_id : (int)    My CADStudy lesion_id        
        StudyID : (int)    CAD StudyID
        redateID : (int)  CAD StudyID Data of exam (format yyyy-mm-dd)
        
        Output
        ======
        """               
        # Create the database: the Session. 
        self.Session = sessionmaker()
        self.Session.configure(bind=myengine)  # once engine is available
        session = self.Session() #instantiate a Session
        
        # for cad_case in session.query(Cad_record).order_by(Cad_record.pt_id): 
        #     print cad_case.pt_id, cad_case.cad_pt_no_txt, cad_case.latest_mutation_status_int    
        for lesion in session.query(mylocaldatabase.Lesion_record).\
            filter(mylocaldatabase.Lesion_record.lesion_id == str(lesion_id)):
            # print results
            if not lesion:
                print "lesion is empty"

        return lesion 
                   

    def query_withT2(self, lesion_id):
        
         # Create the database: the Session. 
        self.Session = sessionmaker()
        self.Session.configure(bind=myengine)  # once engine is available
        session = self.Session() #instantiate a Session

        # for cad_case in session.query(Cad_record).order_by(Cad_record.pt_id): 
        #     print cad_case.pt_id, cad_case.cad_pt_no_txt, cad_case.latest_mutation_status_int    
        for lesion, T2record in session.query(mylocaldatabase.Lesion_record, mylocaldatabase.T2_features).\
            filter(mylocaldatabase.Lesion_record.lesion_id == mylocaldatabase.T2_features.lesion_id).\
            filter(mylocaldatabase.Lesion_record.lesion_id == str(lesion_id)):
            # print results
            if not lesion:
                print "lesion is empty"
            

        return lesion, T2record 
        
           
           

