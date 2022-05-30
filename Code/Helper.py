from psychopy import visual , core , event , gui , sound 
from scipy.stats import logistic
import math
from psychopy.constants import (PLAYING, PAUSED)
from psychopy.hardware import keyboard
from Parameters import*
import pandas as pd
from datetime import datetime 

def getUserInfo():
    again = False
    DialogGui = gui.Dlg(title="User information MDT")
    DialogGui.addField(" Subject Name: ") #0
    DialogGui.addField(" Subject Surname: ") #1
    DialogGui.addField(" Subject Number: ") #2
    DialogGui.addField(" Age: ") #3
    DialogGui.addField('Gender', choices=['Male', 'Female']) #4
    DialogGui.addField('Handedness', choices=['Right', 'Left']) #5
    DialogGui.addField('Task Type', choices=['High Feedback', 'Low Feedback']) #6
    DialogGui.addField('Ages Type', choices=['Child', 'Adult']) #7
    DialogGui.show()
    
    if (DialogGui.data[6]=="High Feedback"):
        Type = "HF"
    if (DialogGui.data[6]=="Low Feedback"):
        Type = "LF"
    if ( not DialogGui.data[2].isdigit() or not DialogGui.data[3].isdigit() ):
        again = True
        DialogGui.data[2] , DialogGui.data[3] = 99999999 , 99999999
        
    if (DialogGui.data[0].isdigit() or DialogGui.data[1].isdigit()):
        again = True

    if (not DialogGui.OK):
        core.quit()
    listInfo = [int(DialogGui.data[2]) , DialogGui.data[0] , DialogGui.data[1] , int(DialogGui.data[3]) , DialogGui.data[4] , DialogGui.data[5] , Type , DialogGui.data[7]]
    if '' or 99999999 in listInfo :
        again = True
    else:
        again = False
    
    return [listInfo , again]

def UserInfo():
    usrInfo , again = getUserInfo()
    while(again):
        usrInfo , again = getUserInfo()
    return usrInfo
    
def Txt_loader(win , title , duration):
    textstimlike=visual.TextBox(
        window=win,
        text=str(title),
        font_size=18,
        font_color=[-1,-1,1],
        color_space='rgb',
        size=(1.8,.1),
        pos=(0.0,.5),
        units='norm')
    imgClk = core.Clock()
    while (imgClk.getTime() < duration):
        textstimlike.draw()
        win.flip() 
    win.flip()
    
def play_encode_num (number):
    if number > 20 and (number%10 != 0) :
        a1 = int(number / 10)
        a2 = number % 10
        playSound(None,"numbers/" + str(a1 * 10) + "_" )
        playSound(None,"numbers/" + str(int(a2)))
    else:
        playSound(None,"numbers/" + str(int(number))  )
        

    
    
def FeedBack(windows,feedb_pic_dir,nCorRes,exptMins,pointsPerMin):
    
    if (feedb_pic_dir[3:8] == "Adult"):
        positions = [(0.29,0.34),(0.12,0.34),(-0.42,0.34),(0.5,0.24)]
    else :
        positions = [(0.38,0.29),(0.18,0.29),(-0.42,0.29),(0.5,0.20)]
        
    
    windows.color = [1,1,1]
    windows.flip()
    
    img = visual.ImageStim(
        win = windows,
        image = "../Asset/Pics/" + feedb_pic_dir + ".PNG" 
        )
#    Left_img = visual.ImageStim(
#        win = windows,
#        image = "../Asset/Pics/" + "L1" + ".PNG" ,
#        pos=(-770 , 0)
#        )
#    Right_img = visual.ImageStim(
#        win = windows,
#        image = "../Asset/Pics/" + "R1" + ".PNG" ,
#        pos=(770 , 0)
#        )
    
    txt_nCorRes = visual.TextBox(
        window=windows,
        text=str(int(nCorRes)),
        font_size=30,
        font_color=[-1,-1,-1],
        color_space='rgb',
        size=(1.9,.3),
        pos=positions[0],
        grid_horz_justification='center',
        units='norm')
    txt_exptMins = visual.TextBox(
        window=windows,
        text=str(int(exptMins)),
        font_size=30,
        font_color=[-1,-1,-1],
        color_space='rgb',
        size=(1.9,.3),
        pos=positions[1],
        grid_horz_justification='center',
        units='norm')
    txt_pointsPerMin = visual.TextBox(
        window=windows,
        text=str(int(pointsPerMin)),
        font_size=30,
        font_color=[-1,-1,-1],
        color_space='rgb',
        size=(1.9,.3),
        pos=positions[2],
        grid_horz_justification='center',
        units='norm')
    txt_score = visual.TextBox(
        window=windows,
        text=str(math.floor(pointsPerMin/5)),
        font_size=30,
        font_color=[-1,-1,-1],
        color_space='rgb',
        size=(1.9,.3),
        pos=positions[3],
        grid_horz_justification='center',
        units='norm')
    img.draw()
#    Left_img.draw()
#    Right_img.draw()
    txt_nCorRes.draw()
    txt_exptMins.draw()
    txt_pointsPerMin.draw()
    txt_score.draw()
    windows.flip()
    
    #fixed 
#    playSound(None, "Fixed/P1_" + feedb_pic_dir[3:8]  )
#    play_encode_num (nCorRes)
#    playSound(None ,"Fixed/P2_" + feedb_pic_dir[3:8]  )
#    play_encode_num (exptMins)
#    playSound(None,"Fixed/P3_" + feedb_pic_dir[3:8]  )
#    play_encode_num (pointsPerMin)
#    playSound(None,"Fixed/P4_" + feedb_pic_dir[3:8]  )
#    play_encode_num (math.floor(pointsPerMin/5))
#    playSound(None,"Fixed/P5_" + feedb_pic_dir[3:8] )   
    
#    if(feedb_pic_dir[:2] == "LF"):
#        playSound( None,"Rest/" +  feedb_pic_dir[3:8]  )
#    else:
#        playSound(None ,"HF_" + feedb_pic_dir[3:8] + "/" +  feedb_pic_dir )
        
    windows.color = [-1,-1,-1]



    
def playFeedbackGif(windows , movieNamel , movieNamer ):

    
    movl = visual.MovieStim3(win = windows, filename = "../Asset/gifs/left/" + movieNamel+ ".gif", flipVert=False , colorSpace='rgb'  , pos=(-500 , 0))
    movr = visual.MovieStim3(win = windows, filename = "../Asset/gifs/right/" + movieNamer + ".gif", flipVert=False , colorSpace='rgb'  , pos=(500 , 0))
    print(movieNamel,movieNamer)
    print(movl.duration,movr.duration)
    if (movl.duration > movr.duration) :
        while(movl.status != visual.FINISHED):
            movl.draw()
            movr.draw()
            windows.flip()
            if(movr.status == visual.FINISHED):
                movr.seek(movr.duration - 0.2)
                movr.draw()
                movl.draw()
                windows.flip()
            if(movl.status == visual.FINISHED and movr.status == visual.FINISHED):
                movl.seek(movl.duration - 0.2)
                movr.seek(movr.duration - 0.2)
                movl.draw()
                movr.draw()
                windows.flip()
            
               
                break
#            if(movr.status == visual.FINISHED):
#                while(movl.status != visual.FINISHED):
#                    core.wait(0.5)
#                    movl.draw()
#                    windows.flip()
    else:
        while(movr.status != visual.FINISHED):
            movr.draw()
            movl.draw()
            windows.flip()
            if(movl.status == visual.FINISHED):
                movl.seek(movl.duration - 0.2)
                movr.draw()
                movl.draw()
                windows.flip()
            if(movl.status == visual.FINISHED and movr.status == visual.FINISHED):


                movl.seek(movl.duration - 0.2)
                movr.seek(movr.duration - 0.2)
                movl.draw()
                movr.draw()
                
                windows.flip()
               
                
                break
                
#            if(movl.status == visual.FINISHED):
#                while(movr.status != visual.FINISHED):
#                    core.wait(0.5)
#                    movr.draw()
#                    windows.flip()
    windows.flip()
    core.wait(2)

    
    

    
def playSound(duration,soundName):

    c4 = sound.Sound( "../Asset/voice/" + soundName + ".wav")
    if(duration == None ):
        duration = c4.duration
    c4.play()
#    c4.setSound(value = 500 , secs=duration)
    core.wait(duration)

    
def playVideo(win,movieName):
    mov = visual.MovieStim3(win, filename = "../Asset/Videos/" + movieName + ".mkv", flipVert=False , colorSpace='rgb')
#    mov = visual.MovieStim(win, filename = "../Asset/Videos/" + movieName + ".mp4", flipVert=False)
    while mov.status != visual.FINISHED:
        mov.draw()
        win.flip()
        intrupt = event.getKeys(['q'])
        if 'q' in intrupt:
            break


def Instruction(win,usrInfo):
#    if (usrInfo[7] == "Child"):
#        playVideo(win,'introsheeps')

    displayImg(win,'Intro' + usrInfo[7],0,instr=True,size=None,pos=None)
    event.waitKeys(keyList=['space'] , clearEvents = True)
    
#    playFeedbackVoice(win)
#    playFeedbackGif(win , 'lg28' , 'lg28')
#    core.wait(2)
#    core.quit()
    
def set_window(screen_size,color):

    if color == 'black' or color == 'Black':
        Cwin = [-1,-1,-1]
    elif color == 'white' or color == 'White':
        Cwin = [1,1,1]
    elif color == 'grey' or color == 'Grey' :
        Cwin = [0,0,0]
        
    win = visual.Window(
    monitor=mon,
    size = screen_size,
    fullscr = True,
    units = "pix" , 
    color = Cwin
    )
    
    square_1 = visual.Rect(
        win=win,
        units="pix",
        width=50,
        height=50,
        pos = (-180.0 , 0),
        fillColor=[1, 1, 1],
        lineColor=[1, 1, 1]
    )

    square_2 = visual.Rect(
        win=win,
        units="pix",
        width=50,
        height=50,
        pos = (+180.0 , 0),
        fillColor=[1, 1, 1],
        lineColor=[1, 1, 1]
    )

    
    return win , square_1 , square_2
    
def setDotStim (win,nDots,coherence,
            dotSize,
            speed,
            color,
            opacity):
    
    
    target = visual.DotStim(win, nDots=nDots, coherence=coherence, fieldPos=fieldPos, fieldSize=fieldSize, 
               fieldShape=fieldShape, dotSize=dotSize, dotLife=dotLife, dir=dir, speed=speed, color=color, 
                colorSpace=colorSpace, opacity=opacity, signalDots=signalDots, 
                noiseDots= noiseDots)
                
    return target
    
def displayImg(win ,imgName , duration , instr, size, pos ):
    if size == None and pos == None :
            img = visual.ImageStim(
        win = win,
        image = "../Asset/Pics/" + imgName + ".PNG" 
        )
    else:
        img = visual.ImageStim(
        win = win,
        image = "../Asset/Pics/" + imgName + ".PNG" ,
#        size=(size[0], size[1]),
        pos=(pos[0], pos[1]) 
        )
    if instr :
        img.draw()
        win.flip() 
    else :
        imgClk = core.Clock()
        while (imgClk.getTime() < duration):
            img.draw()
            win.flip()
        
        win.flip()
        
        
def fixation(win,duration):

    fixation_cross = visual.TextStim(win, text="+")        
    fixation_cross.draw()
    win.flip()
    core.wait(duration)
    win.flip()

    
def stimulus(win,dot_stim,square_1,square_2):
    
    square_1.autoDraw = True
    square_2.autoDraw = True
    dot_stim.dir = np.random.choice((180,0))
    kb = keyboard.Keyboard()
    kb.clock.reset()
    kb.clearEvents()
    timer = core.Clock()
    while(True):
        er_t = timer.getTime()
        keys = kb.getKeys(['z' , 'slash'] , waitRelease=True)
        dot_stim.draw()
        win.flip()
        er_t = timer.getTime() - er_t
        if 'z' in keys or 'slash' in keys:
            respons = keys
            break
    square_1.autoDraw = False
    square_2.autoDraw = False
    return keys[0] , keys[0].rt - er_t , respons

def Trial_feedback(win , key , Rt , dot_stim , usrInfo):
    corr = False
    corr_RT = 0
    icorr_RT = 0
    
    if(dot_stim.dir == 0):
        corrAns = 'Right'
    else:
        corrAns = 'Left'
    
    if( not Rt >= Too_late and not Rt <= Too_soon ):
        Speed = "Normal"
        if(dot_stim.dir == 0):
            if (key.name=='z'):
                playSound(0,"End_trial/Wrong")
                displayImg(win , 'Wrong_' + usrInfo[7] , Wrong_feedback_dur , False, None, None )
                icorr_RT = Rt
            else :
                playSound(0,"End_trial/Correct")
                displayImg(win , 'Correct_' + usrInfo[7] , Correct_feedback_dur , False, None, None )
                corr = True
                corr_RT = Rt
        else :
            if (key.name=='slash'):
                playSound(0,"End_trial/Wrong")
                displayImg(win , 'Wrong_' + usrInfo[7] , Wrong_feedback_dur , False, None, None )
                icorr_RT = Rt
            else :
                playSound(0,"End_trial/Correct")
                displayImg(win , 'Correct_' + usrInfo[7] , Correct_feedback_dur , False, None, None )
                corr = True
                corr_RT = Rt
    else:
        if (Rt >= Too_late):
            Speed = "Slow"
            playSound(0,"End_trial/Too_slow")
            displayImg(win , '2late_' + usrInfo[7] , Too_LateorSoon_dur , False, None, None )
            icorr_RT = Rt
        elif (Rt<= Too_soon):
            Speed = "Fast"
            playSound(0,"End_trial/Too_Fast")
            displayImg(win , '2soon_' + usrInfo[7] , Too_LateorSoon_dur , False, None, None )
            icorr_RT = Rt
            
    return corr , corr_RT , icorr_RT , Speed , corrAns
"""
EZ-Diffusion model
Input :
    Pc -> proportion of correct decisions.
    VRT -> variance of response times for correct decisions.
    MRT -> mean response time for correct decisions.
    s ->  scaling parameter = 0.1
Returns :
    drift rate. , boundary separation , nondecision time
    
"""

def EZ_Diffusion(Pc, VRT, MRT, S):
    
    if (VRT == 0) : 
        VRT = 0.01
    # The default value for the scaling parameter s equals 0.1
    s2 = S**2
    L = logistic.ppf(Pc)
    # This gives drift rate.
    x = L * (L*Pc**2 - L*Pc+Pc - 0.5)/VRT  
    v = np.sign(Pc - 0.5) * S * x**(1/4) 
    # This gives boundary separation
    a = s2 * logistic.ppf(Pc) / v
    # This gives nondecision time
    y = -v * a / s2 
    MDT = (a/(2*v)) * (1-np.exp(y))/(1+np.exp(y))
    Ter = MRT - MDT
    
    return v,a,Ter
    
def get_opt_acc(a,v):
    return 1/(1 + (np.exp(-a*v)))


def get_opt_speed(a,v,Ter):
    return (a/(2*v))  *  ( (1-np.exp(-v*a))/(1+np.exp(-v*a)) ) + Ter
    
def ITI (win , duration):

    win.flip()
    core.wait(duration)

    
def training_phase(win,dot_stim,usrInfo,square_1,square_2):
    #training phase 25 trial
    dot_stim.coherence = 0.8
    Ncor_a = 0
    qualified = False
    for j in range(Trainning_phase_blocks):
#        while(not qualified):
        for i in range(Trainning_phase_trials):
            # Fixation -> stimulus -> Feedback -> ITI
            fixation(win,Fixation_dur)
            key , Rt , _ = stimulus(win,dot_stim,square_1 , square_2)
            corr, _ , _ ,_ , _ = Trial_feedback(win , key , Rt , dot_stim , usrInfo)
            ITI(win , ITI_dur)
            
            if (corr) : 
                Ncor_a +=1
            elif (not corr and Ncor_a != correct_sequence):
                Ncor_a = 0
            if (Ncor_a >= correct_sequence):
                qualified = True
        if (not qualified):
            Finish(win)
            core.quit()

def test_phase(win,dot_stim,usrInfo,square_1,square_2):
    #test phase 40 trial
    dot_stim.coherence = coherence
    score_corr = 0
    score_time = 0
    score_list = []
    timeStatdic = {'block' : 0,
                   'trial_time' : 0,
                   'Ncor_a' : 0}
    timeStat = []
    corr_rt = []
    icorr_rt = []
    
    otherAoptions = []
    otherPcs = []
    otherVRTs = []
    otherMRTs = []
    numerator = []
    denominator = []
    otherRewardRates = []
    
    
    
    CorrAns = []
    RTime = []
    userSpeed = []
    Trial = []
    CatNum = []
    userAns = []
    Acc = []
    Mean_mrt = []
    Best_mrt = []
    Mean_pc = [] 
    Best_pc = []
    bestReward = [] 
    userReward = []
    timer = core.Clock()
    for block_c in range(1 , Test_phase_blocks+1):
        Ncor_a = 0
        Nicor_a = 0
        s_time = timer.getTime()
        for trial_c in range(1 , Test_phase_trials+1):
            CatNum.append(block_c)
            # Fixation -> stimulus -> Feedback -> ITI
            fixation(win,Fixation_dur)
            key , Rt , responskey = stimulus(win,dot_stim,square_1 , square_2)
            RTime.append(Rt)
            if 'z' in responskey:
                userAns.append('Left')
            else:
                userAns.append('Right')
                
            corr, corr_RT , icorr_RT ,Speed , corrAns = Trial_feedback(win , key , Rt , dot_stim , usrInfo)
            userSpeed.append(Speed)
            CorrAns.append(corrAns)
            
            if (corr_RT != 0) :
                corr_rt.append(corr_RT)
            if (icorr_RT != 0) :
                icorr_rt.append(icorr_RT)
                
            ITI(win  , ITI_dur) #interTrialInterval
            
            if (corr) : 
                Ncor_a +=1
                Acc.append(1)
            else :
                Nicor_a +=1
                Acc.append(0)
                
            Trial.append(trial_c)
        trial_time = timer.getTime() - s_time
        timeStatdic['block'] = block_c 
        timeStatdic['trial_time'] = trial_time
        timeStatdic['Ncor_a'] = Ncor_a
        timeStat.append(timeStatdic.copy())
        

        #less than 200 trials    
        if (block_c <= Test_nBlock_bound):
            trialsBack = block_c * Test_phase_trials
            for items in timeStat :
                score_corr+=items['Ncor_a']
                score_time+=items['trial_time']
            score_list.append((score_corr * 60 / score_time, score_time))
            
            if (len(corr_rt)!=0):
                VRT = np.asarray(corr_rt).var()
                MRT = np.asarray(corr_rt).mean()
            else:
                VRT = epsilon
                MRT = epsilon
            Pc =  score_corr / (block_c * Test_phase_trials)
            if (Pc == 1):
                Pc = 1 - (1/(2*block_c *Test_phase_trials))
            if (Pc == 0.5):
                Pc = 0.5 + (1/(2*block_c *Test_phase_trials))
            if (Pc == 0):
                Pc = 0 + (1/(2*block_c *Test_phase_trials))
                
            v,a,Ter = EZ_Diffusion(Pc, VRT, MRT, S=stoch_s)
            
            score_corr, score_time = 0 , 0
            score_list.clear()
            
        #each 200 trials    
        else :
            trialsBack = Test_nBlock_bound
            for i in range (block_c - Test_nBlock_bound , len(timeStat) ) :
                score_corr+=timeStat[i]['Ncor_a']
                score_time+=timeStat[i]['trial_time']
            score_list.append((score_corr * 60 / score_time , score_time)) #accuracy + timespent 
             
            if (len(corr_rt)!=0):
                VRT = np.asarray(corr_rt[-Test_nBlock_bound:]).var()
                MRT = np.asarray(corr_rt[-Test_nBlock_bound:]).mean()
            else:
                VRT = epsilon
                MRT = epsilon
            Pc =  score_corr / (Test_nBlock_bound)
            if (Pc == 1):
                Pc = 1 - (1/(2*Test_nBlock_bound))
            if (Pc == 0.5):
                Pc = 0.5 + (1/(2*Test_nBlock_bound))
            if (Pc == 0):
                Pc =  (1/(2*Test_nBlock_bound))
            v,a,Ter = EZ_Diffusion(Pc, VRT, MRT, S=stoch_s)
            
            score_corr, score_time = 0 , 0 
            score_list.clear()
        

        otherAoptions = np.arange (0.001,1,0.001)
            

        for i in range(len(otherAoptions)):
            otherPcs.append(1 / (1 + np.exp(((- otherAoptions[i] *  v) / np.power(stoch_s,2)))))
            otherVRTs.append((((otherAoptions[i]*np.power(stoch_s,2))/(2*np.power(v,3)))*(((2*((-v*otherAoptions[i])/np.power(stoch_s,2))*np.exp(((-v*otherAoptions[i])/np.power(stoch_s,2))))-(np.exp(2*((-v*otherAoptions[i])/np.power(stoch_s,2))))+1)/np.power((np.exp(((-v*otherAoptions[i])/np.power(stoch_s,2)))+1),2))))
            otherMRTs.append((((otherAoptions[i]/(2*v))*(1-np.exp(((-v*otherAoptions[i])/np.power(stoch_s,2))))/(1+np.exp(((-v*otherAoptions[i])/np.power(stoch_s,2))))) + Ter))
            numerator.append(otherPcs[i]*trialsBack)
            denominator.append((otherMRTs[i]+(ITI_dur)+(Correct_feedback_dur)+((1-otherPcs[i])*(Wrong_feedback_dur)))*trialsBack)
            otherRewardRates.append(numerator[i]/denominator[i])
            
        bestRewardRate = max(otherRewardRates)
        numerator.clear()
        denominator.clear()
#        otherPcs.clear()
#        otherVRTs.clear()
#        otherMRTs.clear()
#        otherRewardRates.clear()        
        RRnum = 0
        RRden = 0
        
        RRnum = Pc*trialsBack
        RRden = ((MRT+(ITI_dur)+(Correct_feedback_dur)+((1-Pc)*(Wrong_feedback_dur)))*trialsBack)
        rewardRate=RRnum/RRden
        
        print(" bestRewardRate " , bestRewardRate , "\n rewardRate " , rewardRate)
        
        #feedback
        for i in range (len(otherAoptions)) :
            if(bestRewardRate == otherRewardRates[i]):
                positionOfBestRR = i
        bestRR_Pc = otherPcs[positionOfBestRR]
        bestRR_MRT = otherMRTs[positionOfBestRR]
        bestRR_VRT = otherVRTs[positionOfBestRR]
        
        meanPc = Pc
        meanMRT = MRT
        meanVRT = VRT
        
        nCorrectResponses = np.round(meanPc*trialsBack)
        
        potentialNTrials = 0
        potentialTrialTime = 0
        
        potentialTrialTime = (otherMRTs[positionOfBestRR]+(ITI_dur)+(Correct_feedback_dur)+((1-otherPcs[positionOfBestRR])*(Wrong_feedback_dur)))
        
        actualPotentialTime = 0
        for i in np.arange(potentialTrialTime,60,potentialTrialTime):
            potentialNTrials += 1
            actualPotentialTime = i
        potentialNCorrectResponses = bestRR_Pc * potentialNTrials
        
        experimentMins = RRden/60
        potentialExperimentMins = actualPotentialTime/60
        
        print("meanMRT " , meanMRT , " bestRR_MRT ", bestRR_MRT)
        MRTdiff = meanMRT - bestRR_MRT
        MRTdiff = float(format(MRTdiff,".3f"))
        print("meanPc " , meanPc , " bestRR_Pc ", bestRR_Pc)
        Pcdiff = meanPc - bestRR_Pc
        for i in range (Test_phase_trials):
            Mean_mrt.append(meanMRT)
            Best_mrt.append(bestRR_MRT)
            Mean_pc.append(meanPc)
            Best_pc.append(bestRR_Pc)
            bestReward.append(bestRewardRate)
            userReward.append(rewardRate)
        Pcdiff = Pcdiff * 100
        Pcdiff = float(format(Pcdiff,".3f"))
        
        correctResponsesDiffPerMin = ((nCorrectResponses/experimentMins) - (potentialNCorrectResponses/potentialExperimentMins))
        
        print("MRTdiff , " , MRTdiff , " Pcdiff , " , Pcdiff)
        if (bestRewardRate >= rewardRate) :
            if(Pcdiff > 0 and MRTdiff > 0): #3
                type_acc = "More"
                type_speed = "Less"
            elif (Pcdiff < 0 and MRTdiff < 0): # 7
                MRTdiff = MRTdiff * -1
                Pcdiff = Pcdiff * -1
                type_acc = "Less"
                type_speed = "More"
            elif(Pcdiff > 0 and MRTdiff < 0): #1
                MRTdiff = MRTdiff * -1
                type_acc = "More"
                type_speed = "More"
            else :                  #(Pcdiff < 0 and MRTdiff > 0): #9
                Pcdiff = Pcdiff * -1
                type_acc = "Less"
                type_speed = "Less"
                
                #        elif(abs(Pcdiff) <= efficient_thrsh and MRTdiff < 0): #4
                #            MRTdiff = MRTdiff * -1
                #            Pcdiff = abs(Pcdiff)
                #            type_acc = "Suff"
                #            type_speed = "More"
                            
                #        elif(abs(Pcdiff) <= efficient_thrsh and MRTdiff > 0): #6
                #            Pcdiff = abs(Pcdiff)
                #            type_acc = "Suff"
                #            type_speed = "Less"
                #        elif(abs(Pcdiff) <= efficient_thrsh and abs(MRTdiff) <= efficient_thrsh): #5
                #            MRTdiff = abs(MRTdiff)
                #            Pcdiff = abs(Pcdiff)
                #            type_acc = "Suff"
                #            type_speed = "Suff"            
                #        elif(Pcdiff < 0 and abs(MRTdiff) <= efficient_thrsh): #8
                #            MRTdiff = abs(MRTdiff)
                #            Pcdiff = Pcdiff * -1
                #            type_acc = "Less"
                #            type_speed = "Suff"
                #        else: #2
                #            MRTdiff = abs(MRTdiff)
                #            type_acc = "More"
                #            type_speed = "Suff"
        else:
            type_acc = "Suff"
            type_speed = "Suff"
            
        correctResponsesDiffPerMin = correctResponsesDiffPerMin * -1
        MRTpercent = (MRTdiff/meanMRT)*100
        pointsPerMin = nCorrectResponses/experimentMins
        meanMRT = float(format(meanMRT,".3f"))               
        nCorrectResponses = float(format(nCorrectResponses,".2f")) #. 
        experimentMins = float(format(experimentMins,".3f")) #. 
        pointsPerMin = float(format(pointsPerMin,".3f")) #. 
        correctResponsesDiffPerMin = float(format(correctResponsesDiffPerMin,".3f"))
        meanMRT = float(format(meanMRT,".3f"))
        
        
        # load gifs
        if (usrInfo[6] == "HF"):
            feedb_pic_dir = usrInfo[6] + "_" + usrInfo[7] + "_a_" + type_acc + "_s_" + type_speed
        else:
            feedb_pic_dir = usrInfo[6] + "_" + usrInfo[7]
        
        print("dir : " , feedb_pic_dir)
        FeedBack(win,feedb_pic_dir,nCorrectResponses,experimentMins,pointsPerMin)
        event.waitKeys(keyList=['space'] , clearEvents = True)


    SaveDate(usrInfo , Trial , CatNum , CorrAns  , userAns , RTime ,
    Acc , userSpeed ,coherence , nDots,  Correct_feedback_dur , 
    ITI_dur  , Wrong_feedback_dur ,Mean_mrt , Best_mrt ,Mean_pc ,Best_pc , bestReward, userReward)
        
def Finish(win):
    displayImg(win,'thanks',0,instr=True,size=None,pos=None)
    playSound(None,"End_thanks/Thanks")
    
def Index_check(arr , maxim):
    for i in range(maxim-len(arr)) :
        arr.append(' ')  
        
def SaveDate(UsrInfo , Trial , CatNum , CorrAns  , userAns , RTime , Acc , 
userSpeed ,Coherence , numberofDots,  FeedbackTime , ITI  , ErrorTimeout ,
Mean_mrt , Best_mrt ,Mean_pc ,Best_pc , bestReward, userReward):
    Id = []
    Name = []
    LastName = []
    Age = []
    Gender = []
    DomHand = []
    Type = []
    Adult_Child = []
    coh = []
    dot_num = []
    iti = []
    FDtime = []
    Errtime = []
    
    usrNum, usrName , usrLastName , usrAge , usrGender ,usrHand , taskType , usrType = UsrInfo
    maxim = max(len(Trial),len(CatNum),len(CorrAns),len(Acc),len(RTime),len(userAns) , len(userSpeed) , 
    len(Mean_mrt) , len(Best_mrt) ,len(Mean_pc) ,len(Best_pc) , len(bestReward), len(userReward))
    for i in range (maxim):
        Id.append(usrNum)
        Name.append('') 
        LastName.append('') 
        Age.append('') 
        Gender.append('') 
        DomHand.append('')
        Type.append('')
        Adult_Child.append('')
        coh.append(Coherence)
        dot_num.append(numberofDots)
        iti.append(ITI)
        FDtime.append(FeedbackTime)
        Errtime.append(ErrorTimeout)
    
    Name[0] = usrName
    LastName[0] = usrLastName
    Age[0] = usrAge
    Gender[0] = usrGender
    DomHand[0] = usrHand
    Type[0] = taskType
    Adult_Child[0] = usrType
    if ( ( len(Trial),len(CatNum),len(CorrAns),len(Acc),len(RTime),len(userAns) , len(userSpeed) ,
    len(Mean_mrt) , len(Best_mrt) ,len(Mean_pc) ,len(Best_pc) , len(bestReward), len(userReward)/ 13) != maxim ):
        l = [Trial,CatNum , CorrAns , Acc ,RTime , userAns , userSpeed , Mean_mrt , Best_mrt ,Mean_pc ,Best_pc , bestReward, userReward]
        for ind_l in l:
            Index_check(ind_l,maxim)
    data_dict = {
    "Subject.num" : Id ,
    "Subject.name" : Name ,
    "Subject.surName" : LastName ,
    "Age" : Age ,
    "Gender" : Gender ,
    "Handedness" : DomHand,
    "Subject.Type" : Adult_Child , 
    "FeedbackType" : Type ,
    "Trial" : Trial ,
    "Category.number" :CatNum ,
    "Correct.answer" : CorrAns , 
    "Acuuracy" : Acc ,
    "R_time" : RTime ,
    "User_answer" : userAns ,
    "User_speed" : userSpeed ,
    "Coherence" : coh ,
    "numberofDots"  : dot_num ,
    "ErrorTimeout" : Errtime ,
    "FeedbackTime" : FDtime , 
    "ITI": iti ,
    "Mean_MRT" : Mean_mrt ,
    "Best_MRT" : Best_mrt ,
    "Mean_PC" : Mean_pc ,
    "Best_PC" : Best_pc ,
    "Best_Reward_rate" : bestReward ,
    "Reward_rate" : userReward

    
    }
    
    UserInfoDF = pd.DataFrame(data_dict,columns= ['Subject.num','Subject.name','Subject.surName','Age' , 'Gender' , 'Handedness' , 
    'Subject.Type' , 'FeedbackType' , 'Trial','Category.number'
    ,'Correct.answer','Acuuracy' ,'R_time' , 'User_answer','User_speed' ,
    "Coherence" , "numberofDots" ,  "ErrorTimeout" , "FeedbackTime" , "ITI" , 
    "Mean_MRT" ,"Best_MRT"  ,"Mean_PC"  ,"Best_PC"  ,"Best_Reward_rate"  ,"Reward_rate"
    
    ])
    dir_csv = "../Output_File/"
    UserInfoDF.to_csv( dir_csv + str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))+"__" + str(usrName) + '_' + str(usrNum) + '.csv' ,index=False,header=True , line_terminator='\r\n')
