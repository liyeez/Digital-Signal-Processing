# Digital-Signal-Processing
DSP Project: Human Ear Organ (Basilar Membrane) model

# This project will focus on modelling the spectral analyses carried out by the human cochlea.

# BM_Passive: 
mechanical model of Basilar Membrane, given a sound input (eg. sine signal) outputs displacement of points on the basilar membrane 

# MagResp: 
function that computes the magnitude response given poles and zeros of transfer function

# fft_gamma: 
provided with a specific frequency, computes the gamma formula output and adjusts gain according to design

# Filter_building: 
Main function that dissects an input sound signal into multiple desired peak frequency of membrane and proceeds to build a filter bank that models the behaviour of basilar membrane relying on output from BM_Passive

# outer_middle_filter:
Models the behaviour of the outer/middle of ear that filters raw sound into signal ready to be dissect by inner ear 
