# -*- coding: utf-8 -*-
"""
Check/Update features
USAGE:
from check_update_features import *
check = CheckUpdate()
check.check_update_selected('T2_features_mastertlist.txt', chgRad=1, chgSegfile=1, chgDyn=1, chgMor=1, chgTex=1, chgT2=1, chgStage1=1)
# only T2
check.check_update_selected('T2_features_mastertlist.txt', chgRad=0, chgSegfile=0, chgDyn=0, chgMor=0, chgTex=0, chgT2=1, chgStage1=0)

FLAGS:
chgRad = 1 To Add/change Radiology ----------| 
chgSegfile = 1  To change segmentation----------| if changing segmentation should change all features
chgDyn = 1  To change Dynamic features      |
chgMor = 1  To change Morphology features   |
chgTex = 1  To change Texture features   <--|
chgT2 = 1   To change T2 features
chgStage1 = 1   To change Stage1 features

Created on Wed Jan 28 16:40:49 2015

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
import pylab      
 
from sendNew2_updatedatabase import *
from retrieve_mydatabase import *
from newFeatures import *

from sqlalchemy.orm import sessionmaker
from mylocalbase import localengine
import mylocaldatabase

class CheckUpdate(object):
    """
    USAGE:
    =============
    check = CheckUpdate()
    """
    def __init__(self): 
        """ initialize database session """           
        #  create a top level Session configuration which can then be used throughout
        # Create the Session
        self.Session = sessionmaker()
        self.Session.configure(bind=localengine)  # once engine is available
        
        
    def __call__(self):       
        """ Turn Class into a callable object """
        CheckUpdate() 
        

    def check_update_selected(self, file_masterlist, chgRad, chgSegfile, chgDyn, chgMor, chgTex, chgT2, chgStage1):
           
        # Open filename list
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
            cond = fileline[15][:-1] 
            Diagnosis = fileline[16]     
            sideBreast = fileline[14]
            
            if(lesion_id>272):
                DicomExamNumber = AccessionN
                nameSegment = StudyID+'_'+AccessionN+'_'+str(lesion_id)+'.vtk'
            else:
                ## for old DicomExamNumber  
                nameSegment = StudyID+'_'+AccessionN+'_'+'lesion.vtk'
            
            #############################
            ###### 1) Querying Research database for clinical, pathology, radiology data
            #############################
            RetrieveDB = Retrieve()
            RetrieveDB.querymyDatabase(lesion_id)
            
            Send2DB = SendNewUpdate()
            # Query biomatrix for radiological complete info
            radioinfo = Send2DB.queryRadioData(StudyID, dateID)
            radioinfo = radioinfo.iloc[0]
            
            if chgRad:
                #############################
                ###### 4) Add Radiology data (raw reports and imaing reasons)
                #############################
                self.session = self.Session() #instantiate a Session
                # Send to database rad_records info 
                rad_records = mylocaldatabase.Radiology_record(lesion_id, radioinfo['cad.cad_pt_no_txt'], radioinfo['cad.latest_mutation'], radioinfo['exam.exam_dt_datetime'],
                        radioinfo['exam.mri_cad_status_txt'], radioinfo['exam.comment_txt'], 
                        str(radioinfo['exam.original_report_txt']),
                        radioinfo['exam.sty_indicator_rout_screening_obsp_yn'], 
                        radioinfo['exam.sty_indicator_high_risk_yn'], radioinfo['exam.sty_indicator_high_risk_brca_1_yn'], radioinfo['exam.sty_indicator_high_risk_brca_2_yn'], radioinfo['exam.sty_indicator_high_risk_brca_1_or_2_yn'], 
                        radioinfo['exam.sty_indicator_high_risk_at_yn'], radioinfo['exam.sty_indicator_high_risk_other_gene_yn'],
                        radioinfo['exam.sty_indicator_high_risk_prior_high_risk_marker_yn'], radioinfo['exam.sty_indicator_high_risk_prior_personal_can_hist_yn'], radioinfo['exam.sty_indicator_high_risk_hist_of_mantle_rad_yn'],
                        radioinfo['exam.sty_indicator_high_risk_fam_hist_yn'], radioinfo['exam.sty_indicator_add_eval_as_folup_yn'], radioinfo['exam.sty_indicator_folup_after_pre_exam_yn'], 
                        radioinfo['exam.sty_indicator_pre_operative_extent_of_dis_yn'], radioinfo['exam.sty_indicator_post_operative_margin_yn'], radioinfo['exam.sty_indicator_pre_neoadj_trtmnt_yn'],
                        radioinfo['exam.sty_indicator_prob_solv_diff_img_yn'], radioinfo['exam.sty_indicator_scar_vs_recurr_yn'], radioinfo['exam.sty_indicator_folup_recommend_yn'], 
                        radioinfo['exam.sty_indicator_prior_2_prophy_mast_yn'])
                                      
                self.session.add(rad_records)
                print "\n====================Sending rad_records to update db... " 
                
                # Finally send records to database
                try:
                    self.session.commit()
                except:
                    self.session.rollback()
                    raise
                finally:
                    self.session.close()
            
            
            #############################
            ###### 2) Get T2 info
            #############################        
            img_folder ='Z:/Cristina/MassNonmass'+os.sep+cond+'/'
            path_rootFolder= 'C:/Users/windows/Documents/repoCode-local/stage1features'
            pathSegment = 'C:\Users\windows\Documents\repoCode-local\stage1features\seg+T2'
            
            [path_T2Series, T2SeriesID] = Send2DB.processT2(T2SeriesID, img_folder, StudyID, DicomExamNumber)
            
            #############################                  
            ###### 3) Load segmentation and display
            #############################
            if chgSegfile:
                [series_path, phases_series, self.lesion3D, chgSeg] = Send2DB.checkSegment(path_rootFolder, cond, StudyID, DicomExamNumber, SeriesID, Lesionfile, T2SeriesID, path_T2Series, lesion_id, sideBreast)
                # SEgmentation details
                self.session = self.Session() #instantiate a Session
                # Send to database lesion info
                segment_records = mylocaldatabase.Segment_record(lesion_id, Send2DB.loadDisplay.lesion_bounds[0], Send2DB.loadDisplay.lesion_bounds[1], Send2DB.loadDisplay.lesion_bounds[2], Send2DB.loadDisplay.lesion_bounds[3], Send2DB.loadDisplay.lesion_bounds[4], Send2DB.loadDisplay.lesion_bounds[5],
                                                    Send2DB.loadDisplay.no_pts_segm, Send2DB.loadDisplay.VOI_vol, Send2DB.loadDisplay.VOI_surface, Send2DB.loadDisplay.VOI_efect_diameter, str(list(Send2DB.loadDisplay.lesion_centroid)), str(Send2DB.lesion_centroid_ijk))
                                                   
                self.session.add(segment_records)
                print "\n====================Sending segment_records to update db... " 
                
                # Finally send records to database
                try:
                    self.session.commit()
                except:
                    self.session.rollback()
                    raise
                finally:
                    self.session.close()                                                    
            else:
                [series_path, phases_series, self.lesion3D] = Send2DB.loadSegment(path_rootFolder, cond, StudyID, DicomExamNumber, SeriesID, Lesionfile, T2SeriesID, path_T2Series, lesion_id, sideBreast)
            

            if chgDyn:
                #############################
                ###### 4) Extract Dynamic features
                #############################
                [dyn_inside, dyn_contour] = Send2DB.extract_dyn(series_path, phases_series, self.lesion3D)
                 
                self.session = self.Session() #instantiate a Session
                # Send to database lesion info
                dyn_records = mylocaldatabase.Dynamic_features(lesion_id, dyn_inside['A.inside'],dyn_inside['alpha.inside'], dyn_inside['beta.inside'], dyn_inside['iAUC1.inside'], dyn_inside['Slope_ini.inside'], dyn_inside['Tpeak.inside'],
                         dyn_inside['Kpeak.inside'], dyn_inside['SER.inside'], dyn_inside['maxCr.inside'], dyn_inside['peakCr.inside'], dyn_inside['UptakeRate.inside'], dyn_inside['washoutRate.inside'], dyn_inside['maxVr.inside'],
                         dyn_inside['peakVr.inside'], dyn_inside['Vr_increasingRate.inside'], dyn_inside['Vr_decreasingRate.inside'], dyn_inside['Vr_post_1.inside'],
                         dyn_contour['A.contour'], dyn_contour['alpha.contour'], dyn_contour['beta.contour'], dyn_contour['iAUC1.contour'], dyn_contour['Slope_ini.contour'], dyn_contour['Tpeak.contour'],
                         dyn_contour['Kpeak.contour'], dyn_contour['SER.contour'], dyn_contour['maxCr.contour'], dyn_contour['peakCr.contour'], dyn_contour['UptakeRate.contour'], dyn_contour['washoutRate.contour'], dyn_contour['maxVr.contour'],
                         dyn_contour['peakVr.contour'], dyn_contour['Vr_increasingRate.contour'], dyn_contour['Vr_decreasingRate.contour'], dyn_contour['Vr_post_1.contour'])
                self.session.add(dyn_records)
                print "\n====================Sending dyn_records to update db... " 
                
                # Finally send records to database
                try:
                    self.session.commit()
                except:
                    self.session.rollback()
                    raise
                finally:
                    self.session.close()
        
            if chgMor:
                #############################
                ###### 5) Extract Morphology features
                #############################
                morphofeatures = Send2DB.extract_morph(series_path, phases_series, self.lesion3D)
                # Morphology
                self.session = self.Session() #instantiate a Session
                # Send to database lesion info
                morp_records = mylocaldatabase.Morpho_features(lesion_id, morphofeatures['min_F_r_i'], morphofeatures['max_F_r_i'], morphofeatures['mean_F_r_i'], morphofeatures['var_F_r_i'], morphofeatures['skew_F_r_i'], morphofeatures['kurt_F_r_i'],
                         morphofeatures['iMax_Variance_uptake'], morphofeatures['iiMin_change_Variance_uptake'], morphofeatures['iiiMax_Margin_Gradient'],  morphofeatures['k_Max_Margin_Grad'],
                         morphofeatures['ivVariance'], morphofeatures['circularity'], morphofeatures['irregularity'], morphofeatures['edge_sharp_mean'], morphofeatures['edge_sharp_std'], morphofeatures['max_RGH_mean'], morphofeatures['max_RGH_mean_k'], morphofeatures['max_RGH_var'], morphofeatures['max_RGH_var_k'])
                self.session.add(morp_records)
                print "\n====================Sending morp_records to update db... "
                
                # Finally send records to database
                try:
                    self.session.commit()
                except:
                    self.session.rollback()
                    raise
                finally:
                    self.session.close()
       
        
            if chgTex:            
                #############################        
                ###### 6) Extract Texture features
                #############################
                texturefeatures = Send2DB.extract_text(series_path, phases_series, self.lesion3D) 
                
                self.session = self.Session() #instantiate a Session
                # Send to database lesion info
                tex_records = mylocaldatabase.Texture_features(lesion_id, texturefeatures['texture_contrast_zero'], texturefeatures['texture_contrast_quarterRad'], texturefeatures['texture_contrast_halfRad'], texturefeatures['texture_contrast_threeQuaRad'], 
                          texturefeatures['texture_homogeneity_zero'], texturefeatures['texture_homogeneity_quarterRad'], texturefeatures['texture_homogeneity_halfRad'], texturefeatures['texture_homogeneity_threeQuaRad'], 
                          texturefeatures['texture_dissimilarity_zero'], texturefeatures['texture_dissimilarity_quarterRad'], texturefeatures['texture_dissimilarity_halfRad'], texturefeatures['texture_dissimilarity_threeQuaRad'], 
                          texturefeatures['texture_correlation_zero'], texturefeatures['texture_correlation_quarterRad'], texturefeatures['texture_correlation_halfRad'], texturefeatures['texture_correlation_threeQuaRad'], 
                          texturefeatures['texture_ASM_zero'], texturefeatures['texture_ASM_quarterRad'], texturefeatures['texture_ASM_halfRad'], texturefeatures['texture_ASM_threeQuaRad'], 
                          texturefeatures['texture_energy_zero'], texturefeatures['texture_energy_quarterRad'], texturefeatures['texture_energy_halfRad'], texturefeatures['texture_energy_threeQuaRad'])

                self.session.add(tex_records)
                print "\n====================Sending tex_records to update db... "
                
                # Finally send records to database
                try:
                    self.session.commit()
                except:
                    self.session.rollback()
                    raise
                finally:
                    self.session.close()
                    

            if chgT2: 
                #############################
                # 7) Extract Lesion and Muscle Major pectoralies signal                                   
                ############################# 
                if T2SeriesID != 'NONE':
                    # query existing T2 records      
                    RetrieveDB4T2 = Retrieve()
                    RetrieveDB4T2.querymyDatabase(lesion_id) 
 
                    if(RetrieveDB4T2.lesion.f_T2):                           
                        RetrieveDB.querymyDatabasefor_T2(lesion_id)     
                        line_muscleVOI = RetrieveDB.T2info.bounds_muscleSI
                        line_muscleVOI = line_muscleVOI.rstrip()
                        l = line_muscleVOI[line_muscleVOI.find('[')+1:line_muscleVOI.find(']')].split(",")
                        bounds_muscleSI = [float(l[0]), float(l[1]), float(l[2]), float(l[3]), float(l[4]), float(l[5]) ]
                        print "\n bounds_muscleSI from file:"
                        print bounds_muscleSI
                        [T2_muscleSI, muscle_scalar_range, bounds_muscleSI, T2_lesionSI, lesion_scalar_range, LMSIR, morphoT2features, textureT2features] = Send2DB.T2_loadupdate(T2SeriesID, path_T2Series, self.lesion3D, pathSegment, nameSegment, sideBreast, bounds_muscleSI)
                        #[T2_muscleSI, muscle_scalar_range, bounds_muscleSI, T2_lesionSI, lesion_scalar_range, LMSIR, morphoT2features, textureT2features] = Send2DB.T2_extract(T2SeriesID, path_T2Series, self.lesion3D, pathSegment, nameSegment, sideBreast)
                        # finally pick T2 BIRADS category
                        BIRADST2 = str(RetrieveDB.T2info.find_t2_signal_int)
                           
                    else:
                        # theres no prior record of T2
                        print "\n Needs to pick muscle seed:"
                        [T2_muscleSI, muscle_scalar_range, bounds_muscleSI, T2_lesionSI, lesion_scalar_range, LMSIR, morphoT2features, textureT2features] = Send2DB.T2_extract(T2SeriesID, path_T2Series, self.lesion3D, pathSegment, nameSegment, sideBreast)
                        # finally pick T2 BIRADS category
                        BIRADST2 = "check biomatrix"
                    
                    self.session = self.Session() #instantiate a Session
                    # Send to database lesion info
                    t2_records = mylocaldatabase.T2_features(lesion_id, BIRADST2, str(list(Send2DB.load.T2dims)), str(list(Send2DB.load.T2spacing)), str(Send2DB.load.T2fatsat), mean(T2_muscleSI), std(T2_muscleSI), str(muscle_scalar_range), str(bounds_muscleSI), mean(T2_lesionSI), std(T2_lesionSI), str(lesion_scalar_range), LMSIR, 
                                            morphoT2features['T2min_F_r_i'], morphoT2features['T2max_F_r_i'], morphoT2features['T2mean_F_r_i'], morphoT2features['T2var_F_r_i'], morphoT2features['T2skew_F_r_i'], morphoT2features['T2kurt_F_r_i'], morphoT2features['T2grad_margin'], morphoT2features['T2grad_margin_var'], morphoT2features['T2RGH_mean'], morphoT2features['T2RGH_var'], 
                                            textureT2features['T2texture_contrast_zero'], textureT2features['T2texture_contrast_quarterRad'], textureT2features['T2texture_contrast_halfRad'], textureT2features['T2texture_contrast_threeQuaRad'], 
                                            textureT2features['T2texture_homogeneity_zero'], textureT2features['T2texture_homogeneity_quarterRad'], textureT2features['T2texture_homogeneity_halfRad'], textureT2features['T2texture_homogeneity_threeQuaRad'], 
                                            textureT2features['T2texture_dissimilarity_zero'], textureT2features['T2texture_dissimilarity_quarterRad'], textureT2features['T2texture_dissimilarity_halfRad'], textureT2features['T2texture_dissimilarity_threeQuaRad'], 
                                            textureT2features['T2texture_correlation_zero'], textureT2features['T2texture_correlation_quarterRad'], textureT2features['T2texture_correlation_halfRad'], textureT2features['T2texture_correlation_threeQuaRad'], 
                                            textureT2features['T2texture_ASM_zero'], textureT2features['T2texture_ASM_quarterRad'], textureT2features['T2texture_ASM_halfRad'], textureT2features['T2texture_ASM_threeQuaRad'], 
                                            textureT2features['T2texture_energy_zero'], textureT2features['T2texture_energy_quarterRad'], textureT2features['T2texture_energy_halfRad'], textureT2features['T2texture_energy_threeQuaRad'])
                    self.session.add(t2_records)                
                    print "\n====================Sending t2_records to update db... "
                    
                    # Finally send records to database
                    try:
                        self.session.commit()
                    except:
                        self.session.rollback()
                        raise
                    finally:
                        self.session.close()
                                     

            
            
            if chgStage1:
                #############################
                ###### 8) Extract new features from each DCE-T1 and from T2 using segmented lesion
                #############################
                newfeatures = newFeatures(Send2DB.load, Send2DB.loadDisplay)
                [deltaS, t_delta, centerijk] = newfeatures.extract_MRIsamp(series_path, phases_series, self.lesion3D, T2SeriesID)
                
                # generate nodes from segmantation 
                [nnodes, curveT, earlySE, dce2SE, dce3SE, lateSE, ave_T2, prop] = newfeatures.generateNodesfromKmeans(deltaS['i0'], deltaS['j0'], deltaS['k0'], deltaS, centerijk, T2SeriesID)    
                [kmeansnodes, d_euclideanNodes] = prop
                
                # pass nodes to lesion graph
                G = newfeatures.createGraph(nnodes, curveT, prop)                   
            
                [degreeC, closenessC, betweennessC, no_triangles, no_con_comp] = newfeatures.analyzeGraph(G)        
                network_measures = [degreeC, closenessC, betweennessC, no_triangles, no_con_comp]
                        
                self.session = self.Session() #instantiate a Session
                # Send to database lesion info
                stage1_records = mylocaldatabase.Stage1_record(lesion_id, d_euclideanNodes, earlySE, dce2SE, dce3SE, lateSE, ave_T2, network_measures)
                self.session.add(stage1_records)
                print "\n====================Sending stage1_records to update db... "
                
                # Finally send records to database
                try:
                    self.session.commit()  
                except:
                    self.session.rollback()
                    raise
                finally:
                    self.session.close()
                            
            ## continue to next case
            line = file_ids.readline()
            print line
                        
        
        file_ids.close()            
        return

