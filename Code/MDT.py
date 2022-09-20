from psychopy import visual , core , event 
from Parameters import *
from Utils import *

usrInfo , Fname = UserInfo() #getting userInfo
win , square_1  , square_2 = set_window(screen_size = [widthPix, heightPix] , color = 'black') #preparing window
Instruction(win,usrInfo) # wpop-up instructions
dot_stim = setDotStim(win,nDots,coherence,dotSize,speed,color,opacity) #setting up dot stimuli
training_phase(win,dot_stim,usrInfo,square_1,square_2,Fname) # training phase
test_phase(win,dot_stim,usrInfo,square_1,square_2,Fname) # test phase
Finish(win)
