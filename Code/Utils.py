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


def select_gif(diff , Type , status):


    if(status == "p"):
        if(abs(diff) - 1 < 0):
            gif_n = np.round(abs(diff ) * 18)
        else:
            gif_n = np.round(((abs(diff) - 1 )/2)*18 + 23)
        
    if(status == "s"):
        if(abs(diff) - 1 < 0):
            gif_n =abs(np.round((abs(diff))*(-18) + 41  ))
        else:
            gif_n =abs( 18 - int(np.round(((abs(diff) - 1 )/2)*18)))
    if(gif_n > 41 ):
        gif_n = 41
    if(gif_n < 0):
        gif_n = 0
    if(gif_n == 22):
        gif_n = 23
    if(gif_n == 19):
        gif_n = 18
        
    if(Type == "Suff"):
        if(status == "p"):
            if(abs(diff) - 1 < 0):
                gif_n = np.round(abs(diff ) * 1 + 19)
            else:
                gif_n = np.round(((abs(diff) - 1) / 2 )*1 + 22)
            
        if(status == "s"):
            if(abs(diff) - 1 < 0):
                gif_n = np.round(((abs(diff)))*(-1) + 22)
            else:
                gif_n =abs(19 -  np.round((abs(diff) - 1 )/2 * 19))
            
            
            
        if(gif_n > 22):
            gif_n = 22
        if(gif_n < 19):
            gif_n = 19
            
                
    gif_dir = Type +"_"+ str(gif_n)

    return int(gif_n)
    
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
    done = 0
    a2 = number
    a3 = number
    a4 = number
    if(number < 20):
        playSound(None,"numbers/o" + str(int(number)))
        done = 1
    else:
        if (number >= 1000 and (number%1000 != 0)):
            a1 = int(number / 1000) 
            a2 = number % 1000   
            a3=a2
            playSound(None,"numbers/o" + str(a1 * 1000) + "_" )
            done = 0
        elif(number >= 1000 and (number%1000 == 0)):
            a1 = int(number / 1000)
            playSound(None,"numbers/o" + str(a1 * 1000) )
            done = 1
        if(done != 1):
            if(a2 >= 100 and (number%100 != 0)):
                a1 = int(a2 / 100) #1 2 3 4 5 
                a3 = a2 % 100   #23  5
                playSound(None,"numbers/o" + str(a1 * 100) + "_" )
                done = 0
            
            elif(a2 >= 100 and (a2%100 == 0)):
                a1 = int(a2 / 100)
                playSound(None,"numbers/o" + str(a1 * 100) )
                done = 1
            if(done != 1):
                if(a3 >= 20 and (a3%10 != 0)):
                    a1 = int(a3 / 10) #1 2 3 4 5 
                    a4 = a3 % 10   #5
                    playSound(None,"numbers/o" + str(a1 * 10) + "_" )
                    playSound(None,"numbers/o" + str(int(a4)))  
                    done = 1
                elif(a3 >= 20 and (a3%10 == 0)):
                    a1 = int(a3 / 10)
                    playSound(None,"numbers/o" + str(a1 * 10) )
                    done = 1
                elif(a3 < 20):
                    playSound(None,"numbers/o" + str(int(a3)))
                    done = 1
                

    
    
def FeedBack(windows,feedb_pic_dir,nCorRes,exptMins,pointsPerMin ,lrfeedb):
    
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
    Left_img = visual.ImageStim(
        win = windows,
        image = "../Asset/Pics/left - speed result/" + lrfeedb[0] + ".PNG" ,
        pos=(-800 , 0)
        )
    Right_img = visual.ImageStim(
        win = windows,
        image = "../Asset/Pics/right -accuracy result/" + lrfeedb[1] + ".PNG" ,
        pos=(800 , 0)
        )

    
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
    if (feedb_pic_dir[:3]=="HF"):
        Left_img.draw()
        Right_img.draw()
    txt_nCorRes.draw()
    txt_exptMins.draw()
    txt_pointsPerMin.draw()
    txt_score.draw()
    windows.flip()
    
    #fixed 
    playSound(None, "Fixed/P1_" + feedb_pic_dir[3:8]  )
    play_encode_num (nCorRes)
    playSound(None ,"Fixed/P2_" + feedb_pic_dir[3:8]  )
    play_encode_num (exptMins)
    playSound(None,"Fixed/P3_" + feedb_pic_dir[3:8]  )
    play_encode_num (pointsPerMin)
    playSound(None,"Fixed/P4_" + feedb_pic_dir[3:8]  )
    play_encode_num (math.floor(pointsPerMin/5))
    playSound(None,"Fixed/P5_" + feedb_pic_dir[3:8] )   
    
    if(feedb_pic_dir[:2] == "LF"):
        playSound( None,"Rest/" +  feedb_pic_dir[3:8]  )
    else:
        playSound(None ,"HF_" + feedb_pic_dir[3:8] + "/" +  feedb_pic_dir )
        
    windows.color = [-1,-1,-1]


    
def playSound(duration,soundName):

    c4 = sound.Sound( "../Asset/voice/" + soundName + ".wav")
    if(duration == None ):
        duration = c4.duration
    c4.play()
    core.wait(duration)

    
def playVideo(win,movieName):
    mov = visual.MovieStim3(win, filename = "../Asset/Videos/" + movieName + ".mp4", flipVert=False , colorSpace='rgb',size=(widthPix,heightPix))
#    mov = visual.MovieStim(win, filename = "../Asset/Videos/" + movieName + ".mp4", flipVert=False)
    while mov.status != visual.FINISHED:
        mov.draw()
        win.flip()
        intrupt = event.getKeys(['q'])
        if 'q' in intrupt:
            break


def Instruction(win,usrInfo):
    if (usrInfo[7] == "Child"):
        playVideo(win,'introsheeps')
        
    win.color = [1,1,1]
    win.flip()
    displayImg(win,'Intro_' + usrInfo[7],0,instr=True,size=None,pos=None)
    win.color = [-1,-1,-1]
    event.waitKeys(keyList=['space'] , clearEvents = True)
    

    
def set_window(screen_size,color):

    if color == 'black' or color == 'Black':
        Cwin = [-1,-1,-1]
    elif color == 'white' or color == 'White':
        Cwin = [1,1,1]
    elif color == 'grey' or color == 'Grey' :
        Cwin = [0,0,0]
        
    win = visual.Window(
#    monitor=mon,
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
        keys = kb.getKeys(['z' , 'slash' , 'q'] , waitRelease=False)
        dot_stim.draw()
        win.flip()
        er_t = timer.getTime() - er_t
        if 'z' in keys or 'slash' in keys:
            respons = keys
            break
        if 'q' in keys :
            core.quit()
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
            icorr_RT = Rt + Too_soon
            
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
    
    
    
    CorrAns_csv = []
    RTime_csv = []
    userSpeed_csv = []
    Trial_csv = []
    CatNum_csv = []
    userAns_csv = []
    Acc_csv = []
    Mean_mrt_csv = []
    Best_mrt_csv = []
    Mean_pc_csv = [] 
    Best_pc_csv = []
    bestReward_csv = [] 
    userReward_csv = []
    
    timer = core.Clock()
    for block_c in range(1 , Test_phase_blocks+1):
        Ncor_a = 0
        Nicor_a = 0
        s_time = timer.getTime()
        if (usrInfo[7] == "Child"):
            displayImg(win ,"sheep stages/Stage_" + str(block_c) , 4 , False, None, None )
        for trial_c in range(1 , Test_phase_trials+1):
            CatNum_csv.append(block_c)
            # Fixation -> stimulus -> Feedback -> ITI
            fixation(win,Fixation_dur)
            key , Rt , responskey = stimulus(win,dot_stim,square_1 , square_2)
            RTime_csv.append(Rt)
            if 'z' in responskey:
                userAns_csv.append('Left')
            else:
                userAns_csv.append('Right')
                
            corr, corr_RT , icorr_RT ,Speed , corrAns = Trial_feedback(win , key , Rt , dot_stim , usrInfo)
            userSpeed_csv.append(Speed)
            CorrAns_csv.append(corrAns)
            
            if (corr_RT != 0) :
                corr_rt.append(corr_RT)
            if (icorr_RT != 0) :
                corr_rt.append(icorr_RT)
                
            ITI(win  , ITI_dur) #interTrialInterval
            
            if (corr) : 
                Ncor_a +=1
                Acc_csv.append(1)
            else :
                Nicor_a +=1
                Acc_csv.append(0)
                
            Trial_csv.append(trial_c)
        trial_time = timer.getTime() - s_time
        timeStatdic['block'] = block_c 
        timeStatdic['trial_time'] = trial_time
        timeStatdic['Ncor_a'] = Ncor_a
        timeStat.append(timeStatdic.copy())
        

        #less than 200 trials    
        if (block_c <= Test_nBlock_bound ):
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
            trialsBack = Test_nBlock_bound * Test_phase_trials
            for i in range (block_c - Test_nBlock_bound , len(timeStat) ) :
                score_corr+=timeStat[i]['Ncor_a']
                score_time+=timeStat[i]['trial_time']
            score_list.append((score_corr * 60 / score_time , score_time)) #accuracy + timespent 
             
            if (len(corr_rt)!=0):
                VRT = np.asarray(corr_rt[-trialsBack:]).var()
                MRT = np.asarray(corr_rt[-trialsBack:]).mean()
            else:
                VRT = epsilon
                MRT = epsilon
            Pc =  score_corr / (trialsBack)
            if (Pc == 1):
                Pc = 1 - (1/(2*trialsBack))
            if (Pc == 0.5):
                Pc = 0.5 + (1/(2*trialsBack))
            if (Pc == 0):
                Pc =  (1/(2*trialsBack))
            v,a,Ter = EZ_Diffusion(Pc, VRT, MRT, S=stoch_s)
            
            score_corr, score_time = 0 , 0 
            score_list.clear()
        
        timeStatdic = { 'otherRewardRates' : 0,
                        'otherMRTs' : 0,
                        'otherPcs' : 0}
#        otherAoptions = np.arange (0.001,1,0.001)
        for i in np.arange (0.001,1,0.001):
            otherAoptions.append(i)
            
        otherPcs.clear()
        otherVRTs.clear()
        otherMRTs.clear()
        otherRewardRates.clear()
        numerator.clear()
        denominator.clear()
      
        for i in range(len(otherAoptions)):
            otherPcs.append(get_opt_acc(otherAoptions[i],v))
            otherVRTs.append((((otherAoptions[i]*np.power(stoch_s,2))/(2*np.power(v,3)))*(((2*((-v*otherAoptions[i])/np.power(stoch_s,2))*np.exp(((-v*otherAoptions[i])/np.power(stoch_s,2))))-(np.exp(2*((-v*otherAoptions[i])/np.power(stoch_s,2))))+1)/np.power((np.exp(((-v*otherAoptions[i])/np.power(stoch_s,2)))+1),2))))
            otherMRTs.append(get_opt_speed(otherAoptions[i],v,Ter))                  
            numerator.append(otherPcs[i]*trialsBack)
            denominator.append((otherMRTs[i]+(ITI_dur)+(Correct_feedback_dur)+((1-otherPcs[i])*(Wrong_feedback_dur)))*trialsBack)
            otherRewardRates.append(numerator[i]/denominator[i])
            
        bestRewardRate = max(otherRewardRates)

        
        RRnum = 0
        RRden = 0
        
        RRnum = Pc*trialsBack
        RRden = ((MRT+(ITI_dur)+(Correct_feedback_dur)+((1-Pc)*(Wrong_feedback_dur)))*trialsBack)
        rewardRate=RRnum/RRden
        print("RRden      " , RRden)
     
        
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
    

        MRTdiff = meanMRT - bestRR_MRT
        MRTdiff = float(format(MRTdiff,".3f"))

        Pcdiff = meanPc - bestRR_Pc

        for i in range (Test_phase_trials):
            if (i == Test_phase_trials-1):
                Mean_mrt_csv.append(meanMRT)
                Best_mrt_csv.append(bestRR_MRT)
                Mean_pc_csv.append(meanPc)
                Best_pc_csv.append(bestRR_Pc)
                bestReward_csv.append(bestRewardRate)
                userReward_csv.append(rewardRate)
            else :
                Mean_mrt_csv.append("-")
                Best_mrt_csv.append("-")
                Mean_pc_csv.append("-")
                Best_pc_csv.append("-")
                bestReward_csv.append("-")
                userReward_csv.append("-")
        Pcdiff = Pcdiff * 100
        Pcdiff = float(format(Pcdiff,".3f"))
        
        correctResponsesDiffPerMin = ((nCorrectResponses/experimentMins) - (potentialNCorrectResponses/potentialExperimentMins))
     
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
        
        
        feedback_gif_pc = "rf" + str(select_gif((meanPc/bestRR_Pc),type_acc , "p") )
        feedback_gif_sp  = "lf" + str(select_gif((meanMRT/bestRR_MRT),type_speed , "s"))
        
        correctResponsesDiffPerMin = correctResponsesDiffPerMin * -1
        MRTpercent = (MRTdiff/meanMRT)*100
        pointsPerMin = nCorrectResponses/experimentMins
        meanMRT = float(format(meanMRT,".0f"))               
        nCorrectResponses = float(format(nCorrectResponses,".0f")) #. 
        experimentMins = float(format(experimentMins,".0f")) #. 
        pointsPerMin = float(format(pointsPerMin,".0f")) #. 
        correctResponsesDiffPerMin = float(format(correctResponsesDiffPerMin,".0f"))
        meanMRT = float(format(meanMRT,".0f"))
        
        
        # load gifs
        if (usrInfo[6] == "HF"):
            feedb_pic_dir = usrInfo[6] + "_" + usrInfo[7] + "_a_" + type_acc + "_s_" + type_speed
        else:
            feedb_pic_dir = usrInfo[6] + "_" + usrInfo[7]
       
        FeedBack(win,feedb_pic_dir,nCorrectResponses,experimentMins,pointsPerMin , (feedback_gif_sp,feedback_gif_pc))
        event.waitKeys(keyList=['space'] , clearEvents = True)


    SaveDate(usrInfo , Trial_csv , CatNum_csv , CorrAns_csv  , userAns_csv , RTime_csv ,
    Acc_csv , userSpeed_csv ,coherence , nDots,  Correct_feedback_dur , 
    ITI_dur  , Wrong_feedback_dur ,Mean_mrt_csv , Best_mrt_csv ,Mean_pc_csv ,Best_pc_csv , bestReward_csv, userReward_csv)
        
def Finish(win):
    displayImg(win,'thanks',0,instr=True,size=None,pos=None)
    playSound(None,"End_thanks/Thanks")
    
def Index_check(arr , maxim):
    for i in range(maxim-len(arr)) :
        arr.append(' ')  
        
def SaveDate(UsrInfo , Trial_csv , CatNum_csv , CorrAns_csv  , userAns_csv , RTime_csv , Acc_csv , 
userSpeed_csv ,Coherence , numberofDots,  FeedbackTime , ITI  , ErrorTimeout ,
Mean_mrt_csv , Best_mrt_csv ,Mean_pc_csv ,Best_pc_csv , bestReward_csv, userReward_csv):
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
    maxim = max(len(Trial_csv),len(CatNum_csv),len(CorrAns_csv),len(Acc_csv),len(RTime_csv),len(userAns_csv) , len(userSpeed_csv) , 
    len(Mean_mrt_csv) , len(Best_mrt_csv) ,len(Mean_pc_csv) ,len(Best_pc_csv) , len(bestReward_csv), len(userReward_csv))
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
    if ( ( len(Trial_csv),len(CatNum_csv),len(CorrAns_csv),len(Acc_csv),len(RTime_csv),len(userAns_csv) , len(userSpeed_csv) ,
    len(Mean_mrt_csv) , len(Best_mrt_csv) ,len(Mean_pc_csv) ,len(Best_pc_csv) , len(bestReward_csv), len(userReward_csv)/ 13) != maxim ):
        l = [Trial_csv,CatNum_csv , CorrAns_csv , Acc_csv ,RTime_csv , userAns_csv , userSpeed_csv , Mean_mrt_csv , Best_mrt_csv ,Mean_pc_csv ,Best_pc_csv , bestReward_csv, userReward_csv]
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
    "Trial_csv" : Trial_csv ,
    "Category.number" :CatNum_csv ,
    "Correct.answer" : CorrAns_csv , 
    "Acuuracy" : Acc_csv ,
    "R_time" : RTime_csv ,
    "User_answer" : userAns_csv ,
    "User_speed" : userSpeed_csv ,
    "Coherence" : coh ,
    "numberofDots"  : dot_num ,
    "ErrorTimeout" : Errtime ,
    "FeedbackTime" : FDtime , 
    "ITI": iti ,
    "Mean_MRT" : Mean_mrt_csv ,
    "Best_MRT" : Best_mrt_csv ,
    "Mean_PC" : Mean_pc_csv ,
    "Best_PC" : Best_pc_csv ,
    "Best_Reward_rate" : bestReward_csv ,
    "Reward_rate" : userReward_csv

    
    }
    
    UserInfoDF = pd.DataFrame(data_dict,columns= ['Subject.num','Subject.name','Subject.surName','Age' , 'Gender' , 'Handedness' , 
    'Subject.Type' , 'FeedbackType' , 'Trial_csv','Category.number'
    ,'Correct.answer','Acuuracy' ,'R_time' , 'User_answer','User_speed' ,
    "Coherence" , "numberofDots" ,  "ErrorTimeout" , "FeedbackTime" , "ITI" , 
    "Mean_MRT" ,"Best_MRT"  ,"Mean_PC"  ,"Best_PC"  ,"Best_Reward_rate"  ,"Reward_rate"
    
    ])
    dir_csv = "../Output_File/"
    UserInfoDF.to_csv( dir_csv + str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S"))+"__" + str(usrName) + '_' + str(usrNum) + '.csv' ,index=False,header=True , line_terminator='\r\n')
    
