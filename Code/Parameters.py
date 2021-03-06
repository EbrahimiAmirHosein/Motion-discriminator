from psychopy import visual, monitors
import random 
import numpy as np
 
#extra
epsilon = 1e-9

#delays
Fixation_dur = random.uniform(0.5,0.75) #500ms - 700ms fixation duration
ITI_dur = random.uniform(0.4,0.6) #500ms - 700ms  ITI duration
Wrong_feedback_dur = 0.8 
Correct_feedback_dur = 0.5 
Too_LateorSoon_dur = 2.5
Too_late = 5
Too_soon = 0.3

#training phase parameters
Trainning_phase_blocks = 1
Trainning_phase_trials = 25
correct_sequence = 4

#test phase parameters:
Test_phase_blocks = 24 #24
Test_phase_trials = 4 #40
Test_nBlock_bound = 5 #200 5 * 40 = 200
stoch_s = 0.1


#dots stimuli parameters
nDots = 200
coherence = 0.3
dotSize=4.0
speed=5
opacity=1

fieldPos =(0.0, 0.0)
fieldSize = (300.0, 300.0)
color=(255.0,255.0,255.0)

fieldShape='circle'
colorSpace='rgb'
signalDots='same'
noiseDots= 'direction'

dotLife=-1
dir = np.random.choice((180,0))


# Display parameters
widthPix = 1600 # screen width in px
heightPix = 900 # screen height in px
monitorwidth = 40 # monitor width in cm
viewdist = 60. # viewing distance in cm
monitorname = 'SonyG500'
scrn = 0 # 0 to use main screen, 1 to use external screen
mon = monitors.Monitor(monitorname, width=monitorwidth, distance=viewdist)
mon.setSizePix((widthPix, heightPix))
