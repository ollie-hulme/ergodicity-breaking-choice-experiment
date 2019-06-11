#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  2 10:07:24 2018

@author: davidm
"""
import numpy as np
import os, scipy.io
import pandas as pd

# Define data paths
datapath_Multi_txt = os.path.join(os.getcwd(),'data','TxtFiles_multiplicative')
datapath_ADD_txt = os.path.join(os.getcwd(),'data','TxtFiles_additive')
datapaths_txt = [datapath_Multi_txt,datapath_ADD_txt]
txtAppnd = ['_multi','_add']


# Retrieve all relevant data in a loop for each subject(=i) and for multiplicative/
# additive(=j) from txt files

redo = [100.0, 1.0]# The growth factors in the txt files are multiplied with a
# factor of 100 (e.g. 183.02 instead of 1.8302), this list will be used to correct
# for this

datadict = {} # Dictionary to be filled with variables along the script
for j in range(2):
    for i in range(19):
        subjID = str(1+i)
        f = open(os.path.join(datapaths_txt[j], subjID + '_' + '2.txt'),'rb') #2 means active phase
        df = pd.read_csv(f,sep='\t',header=0)
        
        '''
        append which session first
        '''
        if j == 0: #only want to add information once
            if i < 9: #First 9 subjects had multiplicative session on day 1
                datadict.setdefault('MultiplicativeSessionFirst',[]).append(1)
            else:
                datadict.setdefault('MultiplicativeSessionFirst',[]).append(0)
                
        '''
        Retrieve growth factor values and wealth and add to dictionary
        '''
        datadict.setdefault('Gam1_1'+txtAppnd[j],[]).append(np.array(df['Gam1_1'])/redo[j]) # setdefault 
        # function either creates a key-value pair for the dictionary if it 
        # doesn't exist yet and appends if it does exist.
        datadict.setdefault('Gam1_2'+txtAppnd[j],[]).append(np.array(df['Gam1_2'])/redo[j])
        datadict.setdefault('Gam2_1'+txtAppnd[j],[]).append(np.array(df['Gam2_1'])/redo[j])
        datadict.setdefault('Gam2_2'+txtAppnd[j],[]).append(np.array(df['Gam2_2'])/redo[j])
        datadict.setdefault('Wealth'+txtAppnd[j],[]).append(df['earnings'][0]) #the wealth variable
        #contains the theoretical wealth progression were the gambles to be realized
        #- which they never were, not even for payout. Thus, this vector is irrelevant,
        #the only relevant value is the very first entry denoting the wealth the 
        #subject took from the passive to the active session
        
        '''
        Retrieve keypresses and recode into accept/reject left side gamble
        '''        
        datadict.setdefault('KP_Final'+txtAppnd[j],[]).append(np.array(df['KP_Final'])) #if subject
        #pressed more than once (i.e. switched between left and right gamble), we take the final keypress as decision        
        datadict.setdefault('KP_Late'+txtAppnd[j],[]).append(np.array(df['KPlate'])) #if subject
        #also pressed after allowed response window, we define that one as the 
        #chosen option (even though the final keypress within response window
        #was highlighted for the subject)
        
        # Code keypresses as binary choices (1=chosen gamble on the left, 0= chosen gamble on the right)
        kp_final = datadict['KP_Final'+txtAppnd[j]][i] #retrieve from dictionary, coded as 8's and 9's.
        #This codes the last of all keypresses made within response window.
        kp_late = datadict['KP_Late'+txtAppnd[j]][i] #retrieve from dictionary, coded as 8's and 9's.
        #This codes keypresses made after response window.
        sublst = []
        for k in range(len(kp_final)):
            if np.isnan(kp_late[k]):#if no late keypress
                if kp_final[k] == 9:
                    sublst.append(1)
                elif kp_final[k] == 8:
                    sublst.append(0)
                else:
                    sublst.append(float('nan'))
            else:
                if kp_late[k] == 9:
                    sublst.append(1)
                elif kp_late[k] == 8:
                    sublst.append(0)
                else:
                    print('found unexpected value')
        
        datadict.setdefault('Choice'+txtAppnd[j],[]).append(sublst)
        
        '''
        Retrieve and calculate RTs - some trials have more than one RT (multiple
        button presses), separated by comma. They are only read as strings and 
        can therefore not be converted directly to floats as done for the other 
        variables above
        '''
        datadict.setdefault('RT_stringFormat'+txtAppnd[j],[]).append(df['RT']) #this still 
        #contains RTs as strings
        RT = datadict['RT_stringFormat'+txtAppnd[j]][i]
        datadict.setdefault('RTlate_stringFormat'+txtAppnd[j],[]).append(df['RTlate'])
        RTlate = datadict['RTlate_stringFormat'+txtAppnd[j]][i]
        
        sublst = []
        first = []
        last = []
        for k in range(len(RT)):
            # Split RTs by comma
            rts = str(RT[k]).split(",") #returns RTs in a list for each trial, separated by comma.
            #Always adds empty cell as final entry
            if len(rts)>1:
                for l in range(len(rts)-1):
                    try:
                        rts[l] = float(rts[l])
                    except:
                        rts[l] = float('nan')
            else:
                rts[0] = float('nan')
                    
            rts_late = str(RTlate[k]).split(",") #no empty cell as final entry here
            for l in range(len(rts_late)):
                try:
                    rts_late[l] = float(rts_late[l])
                except:
                    rts_late[l] = float('nan') 
            
            if np.isnan(rts_late[0]): #if no late keypress
                first.append(rts[0]) #append first reaction time
                try: #if there are more than two entries (more than one RT and 
                #the empty cell), add the final one as last RT
                    last.append(rts[-2:-1][0])
                except: #else the first one also is the last one
                    last.append(rts[0])
            else: #if there is a late keypress
                if np.isnan(rts[0]): #if there also is a non-late keypress
                    first.append(rts[0])
                    last.append(rts_late[-1:][0]) 
                else:
                    first.append(rts_late[0]) #append first reaction time
                    last.append(rts_late[-1:][0])#last RT is last rts_late entry
                    #(which oftentimes only has one entry and thus is the same
                    #as the first)
        
        
        datadict.setdefault('FBOnset'+txtAppnd[j],[]).append(pd.to_numeric(df['FBOnset'],downcast ='float'))
        onset = datadict['FBOnset'+txtAppnd[j]][i]
        firstRT = [a - float(b) for a,b in zip(first, onset)]#RTs were not recorded as 
        #times from decision onset but since experiment onset, thus FBOnset has to be subtracted
        lastRT = [a - float(b) for a,b in zip(last, onset)]
        datadict.setdefault('RT_firstPress'+txtAppnd[j],[]).append(firstRT)
        datadict.setdefault('RT_lastPress'+txtAppnd[j],[]).append(lastRT)
        
        '''
        Calculate expected values/time-average growth for the gambles
        '''
        if j == 0: #multiplicative session
            LinU_Gam1_1_multi = (datadict['Gam1_1'+txtAppnd[j]][i] * datadict['Wealth'+txtAppnd[j]][i])-datadict['Wealth'+txtAppnd[j]][i]
            LinU_Gam1_2_multi = (datadict['Gam1_2'+txtAppnd[j]][i] * datadict['Wealth'+txtAppnd[j]][i])-datadict['Wealth'+txtAppnd[j]][i]
            LinU_Gam2_1_multi = (datadict['Gam2_1'+txtAppnd[j]][i] * datadict['Wealth'+txtAppnd[j]][i])-datadict['Wealth'+txtAppnd[j]][i]
            LinU_Gam2_2_multi = (datadict['Gam2_2'+txtAppnd[j]][i] * datadict['Wealth'+txtAppnd[j]][i])-datadict['Wealth'+txtAppnd[j]][i]
            
            LogU_Gam1_1_multi = np.log(datadict['Gam1_1'+txtAppnd[j]][i])
            LogU_Gam1_2_multi = np.log(datadict['Gam1_2'+txtAppnd[j]][i])
            LogU_Gam2_1_multi = np.log(datadict['Gam2_1'+txtAppnd[j]][i])
            LogU_Gam2_2_multi = np.log(datadict['Gam2_2'+txtAppnd[j]][i])
            
            EU_Lin_Gam1_multi = (LinU_Gam1_1_multi + LinU_Gam1_2_multi)/2.0
            EU_Lin_Gam2_multi = (LinU_Gam2_1_multi + LinU_Gam2_2_multi)/2.0
            delta_EU_Lin_multi = EU_Lin_Gam1_multi-EU_Lin_Gam2_multi
            EU_Lin_LeftBetter_multi = (np.array(EU_Lin_Gam1_multi) - np.array(EU_Lin_Gam2_multi))>0.0001 #some of the TA differences
            #are very close to zero (to about 5 digits behind the comma), and that 
            #probably only due to rounding, in reality they are zero. So in that case
            #TA is not really contradicting EV, we exclude them with this definition
            EU_Lin_LeftBetter_multi = EU_Lin_LeftBetter_multi*1 #convert from boolean to 1/0
            
            EU_Log_Gam1_multi = (np.log(datadict['Gam1_1'+txtAppnd[j]][i]) + np.log(datadict['Gam1_2'+txtAppnd[j]][i]))/2.0
            EU_Log_Gam2_multi = (np.log(datadict['Gam2_1'+txtAppnd[j]][i]) + np.log(datadict['Gam2_2'+txtAppnd[j]][i]))/2.0
            delta_EU_Log_multi = EU_Log_Gam1_multi - EU_Log_Gam2_multi
            EU_Log_LeftBetter_multi = (np.array(EU_Log_Gam1_multi) - np.array(EU_Log_Gam2_multi))>0.0001
            EU_Log_LeftBetter_multi = EU_Log_LeftBetter_multi*1
            
            datadict.setdefault('LinU_Gam1_1'+txtAppnd[j],[]).append(np.array(LinU_Gam1_1_multi))
            datadict.setdefault('LinU_Gam1_2'+txtAppnd[j],[]).append(np.array(LinU_Gam1_2_multi)) 
            datadict.setdefault('LinU_Gam2_1'+txtAppnd[j],[]).append(np.array(LinU_Gam2_1_multi))
            datadict.setdefault('LinU_Gam2_2'+txtAppnd[j],[]).append(np.array(LinU_Gam2_2_multi))
            
            datadict.setdefault('LogU_Gam1_1'+txtAppnd[j],[]).append(np.array(LogU_Gam1_1_multi))
            datadict.setdefault('LogU_Gam1_2'+txtAppnd[j],[]).append(np.array(LogU_Gam1_2_multi)) 
            datadict.setdefault('LogU_Gam2_1'+txtAppnd[j],[]).append(np.array(LogU_Gam2_1_multi))
            datadict.setdefault('LogU_Gam2_2'+txtAppnd[j],[]).append(np.array(LogU_Gam2_2_multi))
            
            datadict.setdefault('EU_Lin_Gam1'+txtAppnd[j],[]).append(np.array(EU_Lin_Gam1_multi))
            datadict.setdefault('EU_Lin_Gam2'+txtAppnd[j],[]).append(np.array(EU_Lin_Gam2_multi)) 
            
            datadict.setdefault('EU_Log_Gam1'+txtAppnd[j],[]).append(np.array(EU_Log_Gam1_multi))
            datadict.setdefault('EU_Log_Gam2'+txtAppnd[j],[]).append(np.array(EU_Log_Gam2_multi))
            
            datadict.setdefault('delta_EU_Lin'+txtAppnd[j],[]).append(np.array(delta_EU_Lin_multi))
            datadict.setdefault('delta_EU_Log'+txtAppnd[j],[]).append(np.array(delta_EU_Log_multi))
            
            datadict.setdefault('EU_Lin_LeftBetter'+txtAppnd[j],[]).append(np.array(EU_Lin_LeftBetter_multi))
            datadict.setdefault('EU_Log_LeftBetter'+txtAppnd[j],[]).append(np.array(EU_Log_LeftBetter_multi))
        else: #additive session
            LinU_Gam1_1_add = datadict['Gam1_1'+txtAppnd[j]][i]
            LinU_Gam1_2_add = datadict['Gam1_2'+txtAppnd[j]][i]
            LinU_Gam2_1_add = datadict['Gam2_1'+txtAppnd[j]][i]
            LinU_Gam2_2_add = datadict['Gam2_2'+txtAppnd[j]][i]
            
            LogU_Gam1_1_add = np.log((datadict['Gam1_1'+txtAppnd[j]][i] + datadict['Wealth'+txtAppnd[j]][i])/datadict['Wealth'+txtAppnd[j]][i])
            LogU_Gam1_2_add = np.log((datadict['Gam1_2'+txtAppnd[j]][i] + datadict['Wealth'+txtAppnd[j]][i])/datadict['Wealth'+txtAppnd[j]][i])
            LogU_Gam2_1_add = np.log((datadict['Gam2_1'+txtAppnd[j]][i] + datadict['Wealth'+txtAppnd[j]][i])/datadict['Wealth'+txtAppnd[j]][i])
            LogU_Gam2_2_add = np.log((datadict['Gam2_2'+txtAppnd[j]][i] + datadict['Wealth'+txtAppnd[j]][i])/datadict['Wealth'+txtAppnd[j]][i])
            
            EU_Lin_Gam1_add = (datadict['Gam1_1'+txtAppnd[j]][i] + datadict['Gam1_2'+txtAppnd[j]][i])/2.0
            EU_Lin_Gam2_add = (datadict['Gam2_1'+txtAppnd[j]][i] + datadict['Gam2_2'+txtAppnd[j]][i])/2.0
            delta_EU_Lin_add = EU_Lin_Gam1_add-EU_Lin_Gam2_add
            EU_Lin_LeftBetter_add = (np.array(EU_Lin_Gam1_add) - np.array(EU_Lin_Gam2_add))>0.0001
            EU_Lin_LeftBetter_add = EU_Lin_LeftBetter_add*1
            
            EU_Log_Gam1_add = 0.5*np.log((datadict['Gam1_1'+txtAppnd[j]][i] + datadict['Wealth'+txtAppnd[j]][i])/datadict['Wealth'+txtAppnd[j]][i]) \
                + 0.5*np.log((datadict['Gam1_2'+txtAppnd[j]][i]+datadict['Wealth'+txtAppnd[j]][i])/datadict['Wealth'+txtAppnd[j]][i])
            EU_Log_Gam2_add = 0.5*np.log((datadict['Gam2_1'+txtAppnd[j]][i]+datadict['Wealth'+txtAppnd[j]][i])/datadict['Wealth'+txtAppnd[j]][i]) \
                + 0.5*np.log((datadict['Gam2_2'+txtAppnd[j]][i]+datadict['Wealth'+txtAppnd[j]][i])/datadict['Wealth'+txtAppnd[j]][i])
            delta_EU_Log_add = EU_Log_Gam1_add - EU_Log_Gam2_add
            EU_Log_LeftBetter_add = (np.array(EU_Log_Gam1_add) - np.array(EU_Log_Gam2_add))>0.0001
            EU_Log_LeftBetter_add = EU_Log_LeftBetter_add*1
            
            datadict.setdefault('LinU_Gam1_1'+txtAppnd[j],[]).append(np.array(LinU_Gam1_1_add))
            datadict.setdefault('LinU_Gam1_2'+txtAppnd[j],[]).append(np.array(LinU_Gam1_2_add)) 
            datadict.setdefault('LinU_Gam2_1'+txtAppnd[j],[]).append(np.array(LinU_Gam2_1_add))
            datadict.setdefault('LinU_Gam2_2'+txtAppnd[j],[]).append(np.array(LinU_Gam2_2_add))
            
            datadict.setdefault('LogU_Gam1_1'+txtAppnd[j],[]).append(np.array(LogU_Gam1_1_add))
            datadict.setdefault('LogU_Gam1_2'+txtAppnd[j],[]).append(np.array(LogU_Gam1_2_add)) 
            datadict.setdefault('LogU_Gam2_1'+txtAppnd[j],[]).append(np.array(LogU_Gam2_1_add))
            datadict.setdefault('LogU_Gam2_2'+txtAppnd[j],[]).append(np.array(LogU_Gam2_2_add))
            
            datadict.setdefault('EU_Lin_Gam1'+txtAppnd[j],[]).append(np.array(EU_Lin_Gam1_add))
            datadict.setdefault('EU_Lin_Gam2'+txtAppnd[j],[]).append(np.array(EU_Lin_Gam2_add)) 
            
            datadict.setdefault('EU_Log_Gam1'+txtAppnd[j],[]).append(np.array(EU_Log_Gam1_add))
            datadict.setdefault('EU_Log_Gam2'+txtAppnd[j],[]).append(np.array(EU_Log_Gam2_add))
            
            datadict.setdefault('delta_EU_Lin'+txtAppnd[j],[]).append(np.array(delta_EU_Lin_add))
            datadict.setdefault('delta_EU_Log'+txtAppnd[j],[]).append(np.array(delta_EU_Log_add))
            
            datadict.setdefault('EU_Lin_LeftBetter'+txtAppnd[j],[]).append(np.array(EU_Lin_LeftBetter_add))
            datadict.setdefault('EU_Log_LeftBetter'+txtAppnd[j],[]).append(np.array(EU_Log_LeftBetter_add))
            
        '''
        Calculate theoretically decisive gambles (gambles where time average 
        growth mandates choice of one gamble, while expected value mandates
        choice of the other gamble). 
        '''
        ev = list(datadict['delta_EU_Lin'+txtAppnd[j]][i]) #EV_Diff = Gamble A(left) - Gamble B(right)
        ta = list(datadict['delta_EU_Log'+txtAppnd[j]][i])
        choice = datadict['Choice'+txtAppnd[j]][i]
        evta = list(zip(ev,ta))
 
        dec_gambs_ind = datadict['EU_Lin_LeftBetter'+txtAppnd[j]][i] + datadict['EU_Log_LeftBetter'+txtAppnd[j]][i] #zero if
        #both say right is better, 2 if both say left is better, 1 if they don't agree
        dec_gambs_ind = [1 if k == 1 else 0 for k in dec_gambs_ind]
#        print dec_gambs_ind
        
        dec_gambs = []
        choice_TA = [] #choices according to TA
        for k in range(len(evta)):
            if dec_gambs_ind[k] == 1: #if a decisive gamble
                dec_gambs.append(evta[k])
                if (evta[k][0]-evta[k][1]) > 0: #if EV_Diff is positive and 
                    #TA_Diff is negative (i.e. EV_Diff - TA_Diff > 0), then 
                    #choosing left would be to choose in accordance with TA
                    #(since EV_ and TA_Diff are defined as difference left -
                    #right)
                    if choice[k] == 1:
                        choice_TA.append(0)
                    else:
                        choice_TA.append(1)
                else:
                    if choice[k] == 1:
                        choice_TA.append(1)
                    else:
                        choice_TA.append(0)
        datadict.setdefault('Decisive_gambles_index'+txtAppnd[j],[]).append(dec_gambs_ind)
        datadict.setdefault('Decisive_gambles'+txtAppnd[j],[]).append(dec_gambs)
        datadict.setdefault('DecGambs_choice_TA'+txtAppnd[j],[]).append(choice_TA)
            
        
        '''
        Calculate which of the two options (left or right) is best (wealth 
        optimising) and whether it was chosen
        '''
        diffList = ['delta_EU_Lin','delta_EU_Log']
        for t in range(2):
            bestChoice = []
            for x in datadict[diffList[t]+txtAppnd[j]][i]:
                if x < 0: #negative difference -> right better option
                    bestChoice.append(0)
                else:
                    bestChoice.append(1)
            bestChosen = [1 if k==l else 0 for k,l in zip(bestChoice, datadict['Choice'+txtAppnd[j]][i])]
            datadict.setdefault('BestChoice_Left1Right0_' + diffList[t] + txtAppnd[j],[]).append(bestChoice)
            datadict.setdefault('BestChosen_' + diffList[t] + txtAppnd[j],[]).append(bestChosen)
            
        '''
        Define non-brainer trials and record decisions
        '''
        allGams = np.vstack((datadict['Gam1_1'+txtAppnd[j]][i],datadict['Gam1_2'+txtAppnd[j]][i],
                                                   datadict['Gam2_1'+txtAppnd[j]][i],datadict['Gam2_2'+txtAppnd[j]][i]))
        noBrainChoice = []
        noBrainers =[]
        for k in range(np.shape(allGams)[1]):
            if len(np.unique(allGams[:,k])) < 4: #if not all 4 gamble values
                #are different, then it's a no-brainer trial
                if j == 0:
                    noBrainers.append(np.round(np.log(allGams[:,k]),4))
                else:
                    noBrainers.append(allGams[:,k])
                noBrainChoice.append(choice[k])
        
        corrNoBrainChoice = []
        noBrainDiff = []
        for k in range(len(noBrainers)):
            #if mean left side lower than right and right chosen or if mean left
            #side higher than right side and left chosen, correct no brainer choice
            if (np.mean(noBrainers[k][:2]) < np.mean(noBrainers[k][2:]) and noBrainChoice[k] == 0) or \
            (np.mean(noBrainers[k][:2]) > np.mean(noBrainers[k][2:]) and noBrainChoice[k] == 1):
                corrNoBrainChoice.append(1)
            else:
                corrNoBrainChoice.append(0)
            unq,unq_cnt = np.unique(noBrainers[k],return_counts = True) #unq returns
            #only the unique values (e.g. [-428., -321., -214.] if no brainer
            #trial has values [-321. -428. -214. -428.]). unq_cnt hten returns how
            #often they appear (here: [2,1,1] since -428 appears twice)
            #Then the line below takes the difference in value between the two
            #values that only appear once
            noBrainDiff.append(np.round(noBrainers[k][noBrainers[k]!=unq[unq_cnt>1]][0]-noBrainers[k][noBrainers[k]!=unq[unq_cnt>1]][1],2))
        datadict.setdefault('NoBrainerGambles'+txtAppnd[j],[]).append(noBrainers)
        datadict.setdefault('NoBrainerChoices'+txtAppnd[j],[]).append(noBrainChoice)
        datadict.setdefault('NoBrainerChoiceCorrect'+txtAppnd[j],[]).append(corrNoBrainChoice)
        datadict.setdefault('NoBrainerValueDifference'+txtAppnd[j],[]).append(noBrainDiff) 


scipy.io.savemat(os.path.join(os.getcwd(),'data','allData.mat'),datadict,oned_as='row')
np.savez(os.path.join(os.getcwd(),'data','allData.npz'),datadict = datadict)