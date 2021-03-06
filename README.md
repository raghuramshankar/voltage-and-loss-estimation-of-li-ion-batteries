# Voltage and loss estimation for Li-ion batteries during drive cycles with various parameter identification methods
Repository for my project as part of the course EEK150 "Electric Power Engineering Design Project" at Chalmers.

A research paper was drafted and is now provisionally accepted at [EPE 2021](http://www.epe2021.com/).

# Overview
There exist two fundamentally different kinds of battery cell models.  One is based on the underlying physics of the operation of the Li-ion cell and capturing the cell dynamics (physics based modelling). Another method approximates the behaviour of the cell’s voltage in response to different input currents(equivalent circuit modelling). Although not as accurate, the equivalent circuit modelling is preferred for applications in battery management systems (BMSs) due to its largely reduced computational complexity compared to the physics based models, which contains much more degrees of freedom. The method to find the RC link parameters for Li-ion batteries from time domain data such as a pulse discharge test. 

The contribution of this work is to quantify the battery RC link model accuracy with different parameterextraction methods. A comparative study is presented where four C-rates are used in the pulse discharge to characterize the cell. The model is validated under transient driving cycles and a physics based model is used as the reference. The comparison is performed on two aspects i.e. terminal voltage estimation as well as loss prediction.

# RC Parameters

## R0:
![r0](https://github.com/raghuramshankar/voltage-and-loss-estimation-for-li-ion-batteries/blob/main/images/r0.jpg)

## R1 & R2:
![r1](https://github.com/raghuramshankar/voltage-and-loss-estimation-for-li-ion-batteries/blob/main/images/r1.jpg) ![r2](https://github.com/raghuramshankar/voltage-and-loss-estimation-for-li-ion-batteries/blob/main/images/r2.jpg)

## C1 & C2:
![c1](https://github.com/raghuramshankar/voltage-and-loss-estimation-for-li-ion-batteries/blob/main/images/c1.jpg) ![c2](https://github.com/raghuramshankar/voltage-and-loss-estimation-for-li-ion-batteries/blob/main/images/c2.jpg)

# Loss Profiles

## WLTC Drive Cycle:
![volt-1](https://github.com/raghuramshankar/voltage-and-loss-estimation-for-li-ion-batteries/blob/main/images/volt-1.jpg) ![volt-2](https://github.com/raghuramshankar/voltage-and-loss-estimation-for-li-ion-batteries/blob/main/images/volt-2.jpg)
