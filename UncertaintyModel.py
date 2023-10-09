# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
from numpy import sqrt, sin, cos
from scipy.integrate import quad

"""For parts of this program, a cartesian coordinate system is used in three dimensions, from a side perspective of the submarine
x - in line with the plane, right and left when looking at the page
y - out of the plane, going into and out of the page
z - vertical, up and down when looking at the page

Orientation is considered with 0 in all values being the submarine lying along the z axis facing upright with the sail on the negative x axis side. 
In-plane rotation will be pitch up and down, positive and negative respectively
Out-of-plane rotation will be yaw right and left, positive and negative respectively
Longitudinal rotation will be roll right and left, positive and negative respectively
"""
# Initial parameters
t = 0 # s - time in seconds
t_prime = 100.0 # s - time in seconds used for integration
def func_a_t(t): #calculates acceleration with respect to time
    return 0 # m/s^2 - acceleration with respect to time in line with the orientation of the submarine
theta_in = 0.0 # rad - initial orientation in the in plane or pitch orientation
theta_out = 0.0 # rad - initial orientation in the out of plane or yaw orientation
theta_long = 0.0 # rad - initial orientation in the longitudinal or roll orientation
x = 0 # m - initial position in the x axis
y = 0 # m - initial position in the y axis
z = 0 # m - initial position in the z axis
v = 10.0 # m/s - initial velocity in line with the accelerometers

#Variables for uncertainty
x_0 = 0 # m - initial position uncertainty in the x axis
y_0 = 0 # m - initial position uncertainty in the y axis
z_0 = 0 # m - initial position uncertainty in the z axis
v_x0 = 0 # m/s - initial velocity uncertainty in the x axis
v_y0 = 0 # m/s - initial velocity uncertainty in the y axis
v_z0 = 0 # m/s - initial velocity uncertainty in the z axis
theta_in0 = 0.0 # rad - initial orientation uncertainty in the in-plane or pitch orientation
theta_out0 = 0.0 # rad - initial orientation uncertainty in the out-of-plane or yaw orientation
theta_long0 = 0.0 # rad - initial orientation uncertainty in the longitudinal or roll orientation
alpha_1 = 1*(10**(-9)) # rad/s - rate of gyroscopic drift
alpha_2 = 5*(10**(-7)) # rad/ms - acceleration dependent gyroscopic drift
sigma_a0 = .01 # m/s^2 - constant accelerometer bias
alpha_a = 1 + (5*(10**(-4))) # unitless - nonlinearity of the accelerometer


#Error equations

##Acceleration
def sigma_a(t_prime): #calculates error due to acceleration. Should only output error from constant bias as acceleration is set to 0
    a_t = func_a_t(t_prime)
    return sigma_a0 + alpha_a*a_t


##Orientation
"""these three equations calculate the rate of the accumulation of orientation uncertainty in each of the three axes"""
def rate_sigma_thetaIn(t_prime): # rate of error accumulation in the plane
    a_t = func_a_t(t_prime)
    return alpha_1 + alpha_2*a_t
def rate_sigma_thetaOut(t_prime): # rate of error accumulation out of the plane
    a_t = func_a_t(t_prime)
    return alpha_1 + alpha_2*a_t
def rate_sigma_thetaLong(t_prime): # rate of longitudinal error accumulation
    a_t = func_a_t(t_prime)
    return alpha_1 + alpha_2*a_t


"""these three equations then integrate the three above equations to the total orientation error as a function of time"""
def func_sigma_thetaIn(t_prime):
    return quad(rate_sigma_thetaIn, 0, t_prime)[0] #total error in the plane
def func_sigma_thetaOut(t_prime):
    return quad(rate_sigma_thetaOut, 0, t_prime)[0] #total error out of the plane
def func_sigma_thetaLong(t_prime):
    return quad(rate_sigma_thetaLong, 0, t_prime)[0] #total longitudinal error


##Velocity
"""These equations calculate the uncertainty in velocity, then compute the velocity by integrating the acceleration and adding uncertainty
This is then converted from the velocity in line with the sub to the velocity in the cartesian plane"""
def func_sigma_v(t_prime):
    return quad(sigma_a, 0, t_prime)[0] #uncertainty in velocity
def func_v(t_prime):
    sigma_v = func_sigma_v(t_prime)
    return quad(func_a_t, 0, t_prime)[0] + sigma_v #velocity when accounting for uncertainty

def v_x(t_prime):
    v = func_v(t_prime)
    return v*float(sin(theta_in))*float(cos(theta_long)) #velocity in the x axis
def v_y(t_prime):
    v = func_v(t_prime)
    return v*float(sin(theta_in))*float(sin(theta_long)) #velocity in the y axis
def v_z(t_prime):
    v = func_v(t_prime)
    return v*float(cos(theta_in))*float(cos(theta_out)) #velocity in the z axis

"""These equations then calculate the velocity uncertainty in the cartesian plane"""
def sigma_vX(t_prime):
    sigma_thetaIn = func_sigma_thetaIn(t_prime)
    sigma_v = func_sigma_v(t_prime)
    return sqrt((sigma_v**2)*(sin(theta_in)**2)+(v**2)*(sigma_thetaIn**2)*(cos(theta_in)**2)) #velocity uncertainty in the x axis
def sigma_vY(t_prime):
    sigma_thetaIn = func_sigma_thetaIn(t_prime)
    return sqrt((v**2)*(sigma_thetaIn**2)*(sin(theta_in)**2)) #velocity uncertainty in the y axis
def sigma_vZ(t_prime):
    sigma_thetaIn = func_sigma_thetaIn(t_prime)
    sigma_v = func_sigma_v(t_prime)
    return sqrt((sigma_v**2)*(cos(theta_in)**2)+(v**2)*(sigma_thetaIn**2)*(sin(theta_in)**2)) #velocity uncertainty in the z axis


##Position
def func_sigma_X(t_prime, t_0 = 0):
    return quad(sigma_vX, t_0, t_prime)[0] #position uncertainty in the x axis
def func_sigma_Y(t_prime, t_0 = 0):
    return quad(sigma_vY, t_0, t_prime)[0] #position uncertainty in the y axis
def func_sigma_Z(t_prime, t_0 = 0):
    return quad(sigma_vZ, t_0, t_prime)[0] #position uncertainty in the z axis


#Outputs
##Plots X uncertainty
x = np.linspace(0, 140.4676, 1250000) #creates the x axis of the plot going from 0 to 140 days (roughly 5 months), in 1.25 million 10 second steps
yvals_cycle = [] #the array for the error data points in a single cycle of error growth
yvals = [] #the array where all 5 error cycles will be stored
tempvar = 0 #the var used to store each error data point before getting added into the array
above1 = False #used for the print statements that give the milestone times for each error size
above10 = False
above100 = False

for i in range(250000): #computes a single cycle of error buildup between recalibrations of the INS system
#the increment i represents 10 second intervals
    tempvar = tempvar + func_sigma_X((i*10)+10,i*10)
    yvals_cycle.append(tempvar)
    if not above1:
        if tempvar > 1: #checks the time at which error goes above 1 m
            print(f"{i*10} seconds, {i/360} hours, {i/(360*24)} days underwater before X error hits 1 m")
            above1 = True
    if not above10:
        if tempvar > 10: #checks the time at which error goes above 10 m
            print(f"{i*10} seconds, {i/360} hours, {i/(360*24)} days underwater before X error hits 10 m")
            above10 = True
    if not above100:
        if tempvar > 100: #checks the time at which error goes above 100 m
            print(f"{i*10} seconds, {i/360} hours, {i/(360*24)} days underwater before X error hits 100 m")
            above100 = True
for i in range(5): #plots the single cycle 5 times to represent the longest theoretical patrol and the error values throughout
    for i in range(250000): #loops through the single cycle for each data point and adds it to the full array
        yvals.append(yvals_cycle[i]) 

#Generating the error plot for the X uncertainty
fig, ax = plt.subplots(figsize=(10, 5), dpi = (100), layout='constrained')
ax.set_xlabel('Time (days)')
ax.set_ylabel('Error (m)')
ax.set_title('Positional error in X direction')
plt.plot(x,yvals)
plt.show()

#Generating the histogram plot based on the X uncertainty
fig, ax = plt.subplots(figsize=(10, 5), dpi = (100), layout='constrained')
ax.set_xlabel('Error (m)')
ax.set_ylabel('% of time above a certain error')
ax.set_title('Histogram of X uncertainty')
plt.hist(yvals, histtype = "stepfilled", cumulative = -1, bins = 50, weights=np.ones(len(yvals)) / len(yvals))
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.show()

##Plots Y uncertainty
x = np.linspace(0, 140.4676, 1250)
yvals_cycle = []
yvals = []
tempvar = 0
above1 = False
above10 = False
above100 = False

for i in range(250):  #computes a single cycle of error buildup between recalibrations of the INS system
#the increment i is in steps of 10000 seconds because the error is always 0 so this speeds up calculation
    tempvar = tempvar + func_sigma_Y((i)+1,i)
    yvals_cycle.append(tempvar)
    if not above1:
        if tempvar > 1: #checks the time at which error goes above 1 m
            print(f"{i*10} seconds, {i/360} hours, {i/(360*24)} days underwater before X error hits 1 m")
            above1 = True
    if not above10:
        if tempvar > 10: #checks the time at which error goes above 10 m
            print(f"{i*10} seconds, {i/360} hours, {i/(360*24)} days underwater before X error hits 10 m")
            above10 = True
    if not above100:
        if tempvar > 100: #checks the time at which error goes above 100 m
            print(f"{i*10} seconds, {i/360} hours, {i/(360*24)} days underwater before X error hits 100 m")
            above100 = True
for i in range(5): #plots the single cycle 5 times to represent the longest theoretical patrol and the error values throughout
    for i in range(250): #loops through the single cycle for each data point and adds it to the full array
        yvals.append(yvals_cycle[i]) 
        
fig, ax = plt.subplots(figsize=(5, 5), dpi=(100), layout='constrained')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Error (m)')
ax.set_title('Positional error in Y direction')
plt.plot(x,yvals)
plt.show()

fig, ax = plt.subplots(figsize=(5, 5), dpi=(100), layout='constrained')
ax.set_xlabel('Error (m)')
ax.set_ylabel('% of time above a certain error')
ax.set_title('Histogram of Y uncertainty')
plt.hist(yvals, histtype = "stepfilled", cumulative = -1, bins = 50, weights=np.ones(len(yvals)) / len(yvals))
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.show()


##Plots Z uncertainty
x = np.linspace(0, 140.4676, 1250000) #creates the x axis of the plot going from 0 to 140 days (roughly 5 months), in 1.25 million 10 second steps
yvals_cycle = [] #the array for the error data points in a single cycle of error growth
yvals = [] #the array where all 5 error cycles will be stored
tempvar = 0 #the var used to store each error data point before getting added into the array
above1 = False #used for the print statements that give the milestone times for each error size
above10 = False
above100 = False

for i in range(250000): #computes a single cycle of error buildup between recalibrations of the INS system
#the increment i represents 10 second intervals
    tempvar = tempvar + func_sigma_Z((i*10)+10,i*10) #adds to the previously calculated value the integral from the new time step to the next time step. 
    yvals_cycle.append(tempvar)
    if not above1:
        if tempvar > 1: #checks the time at which error goes above 1 m
            print(f"{i*10} seconds, {i/360} hours, {i/(360*24)} days underwater before Z error hits 1 m")
            above1 = True
    if not above10:
        if tempvar > 10: #checks the time at which error goes above 10 m
            print(f"{i*10} seconds, {i/360} hours, {i/(360*24)} days underwater before Z error hits 10 m")
            above10 = True
    if not above100:
        if tempvar > 100: #checks the time at which error goes above 100 m
            print(f"{i*10} seconds, {i/360} hours, {i/(360*24)} days underwater before Z error hits 100 m")
            above100 = True
for i in range(5): #plots the single cycle 5 times to represent the longest theoretical patrol and the error values throughout
    for i in range(250000): #loops through the single cycle for each data point and adds it to the full array
        yvals.append(yvals_cycle[i]) 
        
#Generating the error plot for the Z uncertainty
fig, ax = plt.subplots(figsize=(10, 5), dpi = (100), layout='constrained') #sets the plot visual size and resolution for the error plot
ax.set_xlabel('Time (days)') #x axis label
ax.set_ylabel('Error (m)') #y axis label
ax.set_title('Positional error in Z direction') #plot title
plt.plot(x,yvals) #plots the arrays of time steps to the x axis and error data to the y axis
plt.show()

#Generating the histogram plot based on the Z uncertainty
fig, ax = plt.subplots(figsize=(5, 5), dpi=(100), layout='constrained')  #sets the plot visual size and resolution for the histogram
ax.set_xlabel('Error (m)') #x axis label
ax.set_ylabel('% of time above a certain error') #y axis label
ax.set_title('Histogram of Z uncertainty') #plot title
plt.hist(yvals, histtype = "stepfilled", cumulative = -1, bins = 50, weights=np.ones(len(yvals)) / len(yvals)) #creates the histogram with a cumulative scale showing the percent of error data above a certain error scale
plt.gca().yaxis.set_major_formatter(PercentFormatter(1)) #formats the y axis data from raw (0.50) to percentages (50%)
plt.show()


#Write outline for writeup in latex
