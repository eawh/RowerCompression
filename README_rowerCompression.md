{\rtf1\ansi\ansicpg1252\cocoartf1671\cocoasubrtf600
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 Scripts for the rower compression process in both the ellipse and crescent shaped trajectory.\
\
The code for the ellipse trajectory is DrivingPotential_CompressEllipse_Quadri_clean.m\
This code is currently configured to run the left flagellum of the tracked quadriflagellate.\
It should work for any force structure created using cilia tracking code, see 10.5281/zenodo.3518403\
It is worth checking the distance from the trap is consistent with the assumed type of the driving potential, ie an attractive trap should always have the distance decreasing and a repulsive trap has the distance increasing. \
\
The code for the crescent trajectory is DrivingPotential_CompressCrescent_Mouse_clean.m\
This code is configured for the mouse brain cilium from Pelliciotta2020 in PNAS. In this particular case the cilium is not perpendicular cell wall, and instead the power stroke is aligned with the x-direction from the tracking code. This script can be configured for any crescent trajectory associated with a force structure from the same tracking code. If a cell wall can be found, it is better to use the commented section which uses the parallel and perpendicular components as measured relative to the wall. This is just before STEP 1a.\
There are several other points where the code should be checked when changing the input.\
At STEP 2 check that a complete number of cycles is input for the fit\
At STEP 3 check the write indices are describing the endpoints and classifying the inner and outer parts of the trajectory\
At STEP 8 double check that the distance is consistent with the style of trap, i.e. same as ellipse the distance should be decreasing if the trap is attractive.}