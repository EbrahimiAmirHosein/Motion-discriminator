from psychopy import visual , core , event 
from Parameters import*
from Utils import *

 
usrInfo = UserInfo()
win , square_1 , square_2 = set_window(screen_size = [widthPix, heightPix] , color = 'black')
Instruction(win,usrInfo)
dot_stim = setDotStim(win,nDots,coherence,dotSize,speed,color,opacity)
training_phase(win,dot_stim,usrInfo,square_1,square_2)
test_phase(win,dot_stim,usrInfo,square_1,square_2)
Finish(win)
