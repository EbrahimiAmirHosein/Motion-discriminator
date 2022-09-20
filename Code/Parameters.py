from psychopy import visual, monitors
import random
import numpy as np

# extra
LIMIT_DETAILS = 5
epsilon = 1e-9

#durarions
Err_feedback_dur = 0.3
Correct_feedback_dur = 0.5
Too_LateorSoon_dur = 2.5
Too_late = 4
Too_soon = 0.3
beep_delay = 60

# training phase parameters
Trainning_phase_blocks = 1
Trainning_phase_trials = 10
correct_sequence = 4

# test phase parameters:
Test_phase_blocks = 24  # 24
Test_phase_trials = 40  # 40
Test_nBlock_bound = 5   # 5 -> 5 * 40 = 200 
stoch_s = 0.1
STAGE_DUR = 1

mid_rest_time = int(Test_phase_blocks/2) #not been used


# dots stimuli parameters
nDots = 100
coherence = 0.4
dotSize = 6.0
speed = 4
opacity = 1

fieldPos = (0.0, 0.0)
fieldSize = (300.0, 300.0)
color = (255.0, 255.0, 255.0)

fieldShape = 'circle'
colorSpace = 'rgb'
signalDots = 'same'
noiseDots = 'direction'

dotLife = 60
dir = np.random.choice((0, 45))

# Display parameters
widthPix = 1920  # screen width in px
heightPix = 1080  # screen height in px
monitorwidth = 40 # monitor width in cm
viewdist = 60. # viewing distance in cm
monitorname = 'SonyG500'
scrn = 0 # 0 to use main screen, 1 to use external screen
mon = monitors.Monitor(monitorname, width=monitorwidth, distance=viewdist)
mon.setSizePix((widthPix, heightPix))

parent_dir = "../Output_File/"


