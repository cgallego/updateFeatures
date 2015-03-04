# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 16:48:07 2014

@ author (C) Cristina Gallego, University of Toronto
"""
import sys, os
import string
import datetime
from numpy import *
import pandas as pd

import mylocaldatabase
from mylocalbase import Base, myengine

import database
from base import Base, engine

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
        T2record = []
        for T2record in session.query(mylocaldatabase.T2_features).\
            filter(mylocaldatabase.T2_features.lesion_id == str(lesion_id)):
            # print results
            if not T2record:
                print "lesion T2record is empty"
            

        return T2record 
        
           
    def queryBiomatrix(self, session, fStudyID, redateID):
        """
        run : Query by StudyID/AccesionN pair study to local folder NO GRAPICAL INTERFACE. default print to output
        
        Inputs
        ======
        fStudyID : (int)    CAD fStudyID
        redateID : (int)  CAD StudyID Data of exam (format yyyy-mm-dd)
        
        Output
        ======
        """                   
        datainfo = []; 
           
        for cad, exam, in session.query(database.Cad_record, database.Exam_record).\
                    filter(database.Cad_record.pt_id==database.Exam_record.pt_id).\
                    filter(database.Cad_record.cad_pt_no_txt == str(fStudyID)).\
                    filter(database.Exam_record.exam_dt_datetime == str(redateID)).all():
                        
            # print results
            if not cad:
                print "cad is empty"
            if not exam:
                print "exam is empty"          
        
            datainfo.append( [([cad.cad_pt_no_txt, cad.latest_mutation_status_int,
                                    exam.exam_dt_datetime, exam.mri_cad_status_txt, exam.comment_txt,
                                    exam.original_report_txt,
                                    exam.sty_indicator_rout_screening_obsp_yn,
                                    exam.sty_indicator_high_risk_yn,
                                    exam.sty_indicator_high_risk_brca_1_yn,  
                                    exam.sty_indicator_high_risk_brca_2_yn, 
                                    exam.sty_indicator_high_risk_brca_1_or_2_yn, 
                                    exam.sty_indicator_high_risk_at_yn, 
                                    exam.sty_indicator_high_risk_other_gene_yn, 
                                    exam.sty_indicator_high_risk_prior_high_risk_marker_yn, 
                                    exam.sty_indicator_high_risk_prior_personal_can_hist_yn, 
                                    exam.sty_indicator_high_risk_hist_of_mantle_rad_yn, 
                                    exam.sty_indicator_high_risk_fam_hist_yn, 
                                    exam.sty_indicator_add_eval_as_folup_yn, 
                                    exam.sty_indicator_folup_after_pre_exam_yn,
                                    exam.sty_indicator_pre_operative_extent_of_dis_yn,
                                    exam.sty_indicator_post_operative_margin_yn, 
                                    exam.sty_indicator_pre_neoadj_trtmnt_yn,
                                    exam.sty_indicator_prob_solv_diff_img_yn, 
                                    exam.sty_indicator_scar_vs_recurr_yn,
                                    exam.sty_indicator_folup_recommend_yn, 
                                    exam.sty_indicator_prior_2_prophy_mast_yn])] )
            
        ################### Send to table display  
        # add main CAD record table       
        colLabels = ("cad.cad_pt_no_txt", "cad.latest_mutation", "exam.exam_dt_datetime", "exam.mri_cad_status_txt", "exam.comment_txt",
                     "exam.original_report_txt",
                     "exam.sty_indicator_rout_screening_obsp_yn",
                     "exam.sty_indicator_high_risk_yn",
                     "exam.sty_indicator_high_risk_brca_1_yn",  
                     "exam.sty_indicator_high_risk_brca_2_yn", 
                     "exam.sty_indicator_high_risk_brca_1_or_2_yn", 
                     "exam.sty_indicator_high_risk_at_yn", 
                     "exam.sty_indicator_high_risk_other_gene_yn", 
                     "exam.sty_indicator_high_risk_prior_high_risk_marker_yn", 
                     "exam.sty_indicator_high_risk_prior_personal_can_hist_yn", 
                     "exam.sty_indicator_high_risk_hist_of_mantle_rad_yn", 
                     "exam.sty_indicator_high_risk_fam_hist_yn", 
                     "exam.sty_indicator_add_eval_as_folup_yn", 
                     "exam.sty_indicator_folup_after_pre_exam_yn",
                     "exam.sty_indicator_pre_operative_extent_of_dis_yn",
                     "exam.sty_indicator_post_operative_margin_yn", 
                     "exam.sty_indicator_pre_neoadj_trtmnt_yn",
                     "exam.sty_indicator_prob_solv_diff_img_yn", 
                     "exam.sty_indicator_scar_vs_recurr_yn",
                     "exam.sty_indicator_folup_recommend_yn", 
                     "exam.sty_indicator_prior_2_prophy_mast_yn")
                         
        # write output query to pandas frame.
        print len(datainfo)
        dinfo = pd.DataFrame(data=datainfo[0], columns=colLabels)
        print(dinfo['exam.original_report_txt'][0])
            
        return dinfo       

