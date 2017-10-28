# -*- coding: utf-8 -*-
# 
# 
# BEGIN PROGRAM
#
# Get necessary packages
#
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist as AA
from scipy.optimize import curve_fit
import copy
#
# Set error flags
#
final_abs_flag = "false"
different_k_flag = "false" # Not yet implemented
#
# INPUT DATA AND CONVERT TO NUMPY ARRAYS 'time' and 'abs'
#
print()
print("Data will be written to a file having the same name as the data file, ")
print("but with _Analysis.txt at then end.")
name = raw_input("Enter the filename of the data file (it should end in .txt): ")
data_file = open(name, 'r')
data = data_file.readline()
newdata = []
newdata = data.split()
size = len(newdata)
for i in range(size):
        newdata[i] = float(newdata[i])
num_data = size/2
timedata = np.zeros(num_data)
absdata = np.zeros(num_data)            
for i in range(num_data):
        timedata[i]=newdata[2*i]
        absdata[i]=newdata[2*i+1]
# 
# Write the data out to a file (for testing purposes). I will probably remove
# this from the final version
#
outputfilename = str(name[:-4] + "_Analysis.txt")
with open(outputfilename, 'a') as outputfile:
        outputfile.write('Input File: '+name+'\n')
        outputfile.write('Output File: '+outputfilename+'\n')
#
print()
print("If the number of fits selected is 1, then the program will just fit")
print("the original data set and won't do an extrapolation.")
steps = raw_input("Enter the number of extrapolation fits (1-10): ")
steps = int(steps)
if steps < 1:
    steps = 1 
if steps > 10:
    steps = 10
with open(outputfilename, 'a') as outputfile:
        outputfile.write('The number of fits selected is: '+str(steps)+'\n')
#
print()
print("If the enzyme concentration is not known, entering 1 will mean that")
print("the reported kcat/Km is actually k(obs).")
enzconc = raw_input("Enter the enzyme concentration in M (e.g 3.2E-9): ")
enzconc = float(enzconc)
with open(outputfilename, 'a') as outputfile:
        outputfile.write('The enzyme concentration is: '+str(enzconc)+'\n')
#
# Setting up the arrays to hold the results from the two different fits
#
kcat_Km = np.zeros(steps)
kcat_Km2 = np.zeros(steps)
substrate = np.zeros(steps)
stddev = np.zeros(steps)
stddev2 = np.zeros(steps)
starting_abs = np.zeros(steps)
#
for loopcount in range(steps):
#
# SET UP THE CURRENT TIME AND ABS DATA
#
    if loopcount == 0:
        current_startabs = absdata[0]
    else:
        current_startabs = absdata[0] + loopcount*(deltaA_overall)/10 # This may not correspond to 10%rxn
    startpt_counter = 0
    while absdata[startpt_counter] < current_startabs:
        startpt_counter = startpt_counter + 1
    current_num_data = num_data-startpt_counter
    print("startpt_counter = ", startpt_counter)
    print("current_num_data = ", current_num_data)
    abs = np.zeros(current_num_data)
    time = np.zeros(current_num_data)
    for i in range(startpt_counter, num_data):
        abs[i-startpt_counter] = absdata[i]
        time[i-startpt_counter] = timedata[i] # I guess it isn't necessary to reset the data to start at time zero?
    with open(outputfilename, 'a') as outputfile:
        outputfile.write('\n')
        outputfile.write("Fit # "+str(loopcount+1)+'\n')
        outputfile.write("The current number of data points is: ")
        outputfile.write(str(current_num_data)+'\n')
#
# The following was used for trouble shooting
#
#       for i in range(current_num_data):
#           outputfile.write(str(time[i]))
#           outputfile.write('\t')  # insert tab between data
#           outputfile.write(str(abs[i]))
#           outputfile.write('\n')  # end the line
#
# ESTIMATING PARAMETERS FOR CURVE FIT - Mthod 1
#
# LINEAR FIT TO BEGINNING 5% TO ESTIMATE vi
# This is a bit crude, but should be good enough
#
    if loopcount == 0:
        timetwo = []
        abstwo = []
        initfitpts = int(current_num_data*0.05)
        for i in range(initfitpts):
           timetwo.append(time[i])
           abstwo.append(abs[i])
        initialfit = np.polyfit(timetwo, abstwo, 1) #Fitting to a line and getting the slope
        vi = initialfit[0]*1.2
#   
# For all but the first cycle, it will just be using the vi from the previous run
#
# A STUPID, BUT GOOD ENOUGH ESTIMATE FOR k, assuming data span about 5 half-lives
# For all but the first cycle, it will use the previous k
#
    if loopcount == 0:
        k = 5./(time[current_num_data-1]-time[0]) # I know this is stupid, but it works.
    init_abs = abs[0]  # This seems obvious
    initguess = np.zeros(3)
    initguess[0]=abs[0]
    initguess[1]=vi
    initguess[2]=k  #This could be changed to k2 if using second routine
#
# FIT TO CURVE - Method 1
# Here's where the magic happens!
#
    def func(time, init_abs, vi, k):
        return init_abs+(vi)*(1-np.exp(-1.*k*time))/k
    popt, pcov = curve_fit(func , time, abs, p0=initguess)
#
    pcov1=copy.copy(pcov) # Make a copy to have available for printing out StdDev on plot
    init_abs = popt[0]
    vi = popt[1]
    k = popt[2]
    kcat_Km[loopcount] = k/enzconc
#
#  Initially assume that the deltaA_overall corresponds to the average of the
#  last 5% of datapoints - the starting absorbance. This is used to estimate  
#  where the 10% reaction points are for the extrapolation.
#
    if loopcount == 0:
        endstart = int(0.95*(num_data-1))
        counter = 0
        sum = 0
        for i in range(endstart,num_data-1):
            counter = counter + 1
            sum = sum + abs[i]
        deltaA_overall = sum/counter-abs[0]
        print("deltaA_overall = ", deltaA_overall)
    substrate[loopcount] = 1-(abs[0]-absdata[0])/deltaA_overall
    starting_abs[loopcount] = abs[0]
    stddev[loopcount] = kcat_Km[loopcount]*np.sqrt(pcov[2,2])/k
    print(substrate[loopcount], kcat_Km[loopcount], stddev[loopcount])

    print()
    print('[Substrate] = ', substrate[loopcount])
    print('Initial Abs = ', init_abs, ' +/- ', np.sqrt(pcov[0,0]))
    print('Initial Rate = ', vi, ' +/- ', np.sqrt(pcov[1,1]))
    print('Rate constant = ', k, ' +/- ', np.sqrt(pcov[2,2]))
    print('[Enzyme] = ', enzconc)
    print('kcat/Km = ', kcat_Km[loopcount], ' +/- ', kcat_Km[loopcount]*np.sqrt(pcov[2,2])/k)
#
    with open(outputfilename, 'a') as outputfile:
        outputfile.write('[Relative S] = '+str(substrate[loopcount])+'\n')
        outputfile.write("Initial Abs = "+str(init_abs)+' +/- '+str(np.sqrt(pcov[0,0]))+'\n')
        outputfile.write('Initial Rate = '+str(vi)+' +/- '+str(np.sqrt(pcov[1,1]))+'\n')
        outputfile.write('Rate constant = '+str(k)+' +/- '+str(np.sqrt(pcov[2,2]))+'\n')
        outputfile.write('[Enzyme] = '+str(enzconc)+'\n')
        outputfile.write('kcat/Km = '+str(kcat_Km[loopcount])+' +/- '+str(kcat_Km[loopcount]*np.sqrt(pcov[2,2])/k)+'/s/M \n')
#
#
# ESTIMATING PARAMETERS FOR CURVE FIT - Method 2
#
    if loopcount == 0:
        deltaAf = abs[current_num_data-1] - abs[0] #after the first time around just use the previous estimate
    initguess2 = np.zeros(3)
    initguess2[0]=init_abs #Using the init_abs from the Method 1 curve fit should be a good guess
    initguess2[1]=deltaAf #Just use the value determined last tiome
    initguess2[2]=k # Using the k from the Method 1 curve fit should be a good guess!
#
# FIT TO CURVE - Method 2
#
    def func(time, init_abs2, deltaAf, k2):
        return init_abs2+(deltaAf)*(1-np.exp(-1.*k2*time))

    popt, pcov = curve_fit(func , time, abs, p0=initguess2)

    init_abs2 = popt[0]
    deltaAf = popt[1]
    k2 = popt[2]
    kcat_Km2[loopcount] = k2/enzconc
#    
#   substrate[loopcount] is the same as already determined
#    
    stddev2[loopcount] = kcat_Km2[loopcount]*np.sqrt(pcov[2,2])/k2
#
    print(substrate[loopcount], kcat_Km2[loopcount], stddev2[loopcount])
    print()
    print('Initial Abs = ', init_abs2, ' +/- ', np.sqrt(pcov[0,0]))
    print('Delta Abs = ', deltaAf, ' +/- ', np.sqrt(pcov[1,1]))
    print('Final Abs = ', init_abs2+deltaAf)
    print('Rate constant = ', k2, ' +/- ', np.sqrt(pcov[2,2]))
    print('kcat/Km = ', kcat_Km2[loopcount], ' +/- ', kcat_Km2[loopcount]*np.sqrt(pcov[2,2])/k2)
#
    with open(outputfilename, 'a') as outputfile:
        outputfile.write('\n')
        outputfile.write('[Relative S] = '+str(substrate[loopcount])+'\n')
        outputfile.write("Initial Abs = "+str(init_abs2)+' +/- '+str(np.sqrt(pcov[0,0]))+'\n')
        outputfile.write('Delta Abs = '+str(deltaAf)+' +/- '+str(np.sqrt(pcov[1,1]))+'\n')
        outputfile.write('Final Abs = '+str(init_abs2+deltaAf)+'\n')
        outputfile.write('Rate constant = '+str(k2)+' +/- '+str(np.sqrt(pcov[2,2]))+'\n')
        outputfile.write('[Enzyme] = '+str(enzconc)+'\n')
        outputfile.write('kcat/Km = '+str(kcat_Km2[loopcount])+' +/- '+str(kcat_Km2[loopcount]*np.sqrt(pcov[2,2])/k)+'/s/M \n')
#        
# SET UP PLOTS OF DATA
#
# NOT EXACTLY SURE ABOUT THE NEXT 4 LINES, But putting data on plot doesn't work without ax defined.
#
    f = plt.figure()
    f.subplots_adjust(right=0.85)
    ax = AA.Subplot(f,1,1,1)
    f.add_subplot(ax)
#
    plot_title = name[:-4] + str(loopcount+1) + "  Relative [S] = " + str(substrate[loopcount])
    x_axislabel = "Time (sec)"
    y_axislabel = "Absorbance"
    plt.title(plot_title)
    plt.xlabel(x_axislabel)
    plt.ylabel(y_axislabel)
    plt.axis([time[0],time[current_num_data-1],min(abs),min(abs)+1.1*(max(abs)-min(abs))])
    plt.plot(time, abs, color='orange') #The actual data for this fit
#
#  The calculated fits superimposed
#
    pcalc = np.zeros(current_num_data)
    pcalc = init_abs+(vi)*(1-np.exp(-1.*k*time))/k #Method 1 fit
    plt.plot(time,pcalc,'k', linewidth=1.0)
    pcalc = init_abs2+(deltaAf)*(1-np.exp(-1.*k2*time)) #Method 2 fit
    plt.plot(time,pcalc,'b', linewidth=1.0)
#
# PRINT DATA ON THE PLOT from Fit Method 1
#
    show_data = "Initial Abs = {:.3e} +/- {:.2e} AU".format(init_abs, np.sqrt(pcov1[0,0]))
    plt.text(0.6, 0.6, show_data, ha='center', va='center', transform=ax.transAxes)
    show_data = "Initial Rate = {:.3e}  +/- {:.2e} AU/sec".format(vi, np.sqrt(pcov1[1,1]))
    plt.text(0.6, 0.56, show_data, ha='center', va='center', transform=ax.transAxes)
    show_data = "Rate constant, k = {:.3e} +/- {:.2e} /sec".format(k, np.sqrt(pcov1[2,2]))
    plt.text(0.6, 0.52, show_data, ha='center', va='center', transform=ax.transAxes)
    show_data = "[Enzyme] = {:.2e} M" .format(enzconc)
    plt.text(0.6, 0.48, show_data, ha='center', va='center', transform=ax.transAxes)
    show_data = "kcat/Km = {:.3e} +/- {:.2e} /M/sec".format(kcat_Km[loopcount], kcat_Km[loopcount]*np.sqrt(pcov1[2,2])/k)
    plt.text(0.6, 0.44, show_data, ha='center', va='center', transform=ax.transAxes)
#
# PRINT DATA ON THE PLOT -- From Fit Method 2
#
    show_data = "Initial Abs = {:.3e} +/- {:.2e} AU".format(init_abs2, np.sqrt(pcov[0,0]))
    plt.text(0.6, 0.36, show_data, ha='center', va='center', transform=ax.transAxes)
    show_data = "DeltaAbs = {:.3e}  +/- {:.2e} AU".format(deltaAf, np.sqrt(pcov[1,1]))
    plt.text(0.6, 0.32, show_data, ha='center', va='center', transform=ax.transAxes)
    show_data = "Final Abs = {:.3e} AU".format(init_abs2+deltaAf)
    plt.text(0.6, 0.28, show_data, ha='center', va='center', transform=ax.transAxes)
    show_data = "Rate constant, k = {:.3e} +/- {:.2e} /sec".format(k2, np.sqrt(pcov[2,2]))
    plt.text(0.6, 0.24, show_data, ha='center', va='center', transform=ax.transAxes)
    show_data = "kcat/Km = {:.3e} +/- {:.2e} /M/sec".format(kcat_Km2[loopcount], kcat_Km2[loopcount]*np.sqrt(pcov[2,2])/k)
    plt.text(0.6, 0.20, show_data, ha='center', va='center', transform=ax.transAxes)
#
# SAVE THE PLOT AS A PDF FILE
#
    plfilenm = name[:-4] +"plot" + str(loopcount+1) + ".pdf"
    plt.savefig(plfilenm, format='pdf')
    print("Saved as ", plfilenm)
    plt.show()  # May want to disable this.....
    plt.close(plfilenm)
    with open(outputfilename, 'a') as outputfile:
        outputfile.write('Plot saved as '+plfilenm+'\n')
    with open(outputfilename, 'a') as outputfile:
        outputfile.write('\n')
#
# AFTER GOING THROUGH AS MANY CYCLES AS DESIRED, WE NOW DROP DOWN TO DO THE EXTRAPOLATION
# OR JUST DROP OUT IF ONLY ONE CYCLE
#
#
# Check to see that the estimate of the Final Abs from the last fit, matches
# the final absorbance calculated from the last fit. If these differ by more 
# than 1% give a warning.
#
loopcount = loopcount+1
if loopcount > 1:
    print("Final abs observed = ", deltaA_overall+absdata[0], " Final abs calcd = ", init_abs2+deltaAf)
    print("Calcd final abs/observed final abs = ", (init_abs2+deltaAf)/(deltaA_overall+absdata[0]))
    with open(outputfilename, 'a') as outputfile:
        outputfile.write('\n')
        outputfile.write('Calcd final abs/observed final abs = '+str((init_abs2+deltaAf)/(deltaA_overall+absdata[0]))+'\n')
    if (init_abs2+deltaAf)/(deltaA_overall+absdata[0]) > 1.001:
        print("WARNING: the final absorbance calculated differs from the final absorbance observed")
        final_abs_flag = "true"
        with open(outputfilename, 'a') as outputfile:
            outputfile.write('WARNING: the final absorbance calculated differs from the final absorbance observed'+'\n')
    if (init_abs2+deltaAf)/(deltaA_overall+absdata[0]) < 0.999:
        print("WARNING: the final absorbance calculated differs from the final absorbance observed")
        final_abs_flag = "true"
        with open(outputfilename, 'a') as outputfile:
            outputfile.write('WARNING: the final absorbance calculated differs from the final absorbance observed'+'\n')
    nothing = raw_input("Just pausing so you can see the ratio, push return to continue")
    if final_abs_flag == "true":
        extrap_plot = 2
    else:
        extrap_plot = 1 
    for k in range(extrap_plot):
        xvalues = []
        yvalues = []
        svalues = []
        for i in range(loopcount):
            xvalues.append(substrate[i])
            yvalues.append(kcat_Km[i])
            svalues.append(stddev[i])
            print(xvalues[i], yvalues[i])
        extrap_kcat_km = np.polyfit(xvalues, yvalues, 1)
        print('Extrapolated kcat/Km without weighting = ', extrap_kcat_km[1])
#       
#
# Do a weighted linearfit of relative[S] vs 1/(kcat/Km) using modified linearfit
# program
#
# Fit data to a line:
#  y=mx+b
#
# Adapted from Numerical Recipes: The Art of Scientific Computing
# Press, W.H.; Flannery, B.P.; Teukolsky, S.A.; Vetterling, W.T.
# (Cambridge University Press, 1986) Chapt. 14.2 "Fitting Data to a
# Straight Line pp 504-509. "Subroutine FIT" 
#
# If the Std. Dev. are unknown, then set them all to 1. (or setting all to 0 or 
# negative numbers will also work similarly.)
#  
# 
        sx=0
        sy=0
        st2=0
        b=0
        x = np.zeros(loopcount)
        y = np.zeros(loopcount)
        s = np.zeros(loopcount)
        w = np.zeros(loopcount)
#
#
# INPUT DATA and MAKE INITIAL ESTIMATES
#
#
        ss = 0.
        for i in range(loopcount):
            x[i] = substrate[i]
            y[i] = 1.0/(kcat_Km[i]/1000)
            s[i] = stddev[i]/kcat_Km[i]*y[i]
            if s[i] > 0:
                w[i] = 1./s[i]**2
            else:
                w[i] = 1.
            print("{:0.4e} {:0.4e} {:0.4e} {:0.4e}".format(x[i], y[i], s[i], w[i]))
            ss = ss+w[i]
            sx = sx+x[i]*w[i]
            sy = sy+y[i]*w[i]
#
        sxoss = sx/ss
#
        for i in range(loopcount):
            t = (x[i]-sxoss)/(1./np.sqrt(w[i]))
            st2 = st2 + t*t
            b = b + t*y[i]/(1./np.sqrt(w[i]))
#
# SOLVE FOR SLOPE(b), INTERCEPT(a), AND SIGMAS (siga, sigb)
#
        b = b/st2
        a = (sy-sx*b)/ss
        siga = np.sqrt((1.+sx*sx/(ss*st2))/ss)
        sigb = np.sqrt(1./st2)
#
#
# OUTPUT RESULTS
#   
# Convert kcat/Km to /s/mM from /s/M
#
#    print
        print("Slope = {:0.3e} +/- {:0.1e} /s/mM".format(b, sigb))
        print("Intercept = {:0.3e} +/- {:0.1e} s•mM".format(a, siga))
        print("1/Intercept = kcat/Km = {:0.3e} +/- {:0.1e} /s/mM".format(1.0/a, 1.0/a*siga/a))
        print("Slope/intercept = {:0.3e}".format(b/a))
#
        with open(outputfilename, 'a') as outputfile:
            outputfile.write('\n')
            if extrap_plot == 1:
                outputfile.write('Extrapolating based on last 5% of data = Final Absorbance'+'\n')
            else:
                outputfile.write('Extrapolating based on fit Final Absorbance'+'\n')
            outputfile.write('Slope = '+str(b)+' +/- '+str(sigb)+'\n')
            outputfile.write('Intercept = '+str(a)+' +/- '+str(siga)+'\n')
            outputfile.write('1/Intercept = kcat/Km '+str(1.0/a)+' +/- '+str(1.0/a*siga/a)+' /mM/s'+'\n')
            outputfile.write('Slope/Intercept = '+str(b/a)+'\n')
        
#
# PLOT DATA
#
        fig = plt.figure()
        ax = AA.Subplot(fig,1,1,1)
        fig.add_subplot(ax)
#
        
        if k == 1:
            plot_title  = name[:-4]+" -- Extrapolation plot, based on last 5% of data = final absorbance"
        else:
            plot_title  = name[:-4]+" -- Extrapolation plot, based on fit final absorbance"
        x_axislabel = "Relative [Substrate]"
        y_axislabel = "Km/kcat (s*mM)"
#
        plt.plot(x, y, 'ko')  # Putting the data points on the plot
        plt.title(plot_title)
        plt.xlabel(x_axislabel)
        plt.ylabel(y_axislabel)
        plt.axis([0,1.1*max(x),0.9*min(y),1.1*max(y)])
        ty = []
        tx = []
        tx.append(0.)
        ty.append(a)
        tx.append(1.1*max(x))
        ty.append(a+1.1*max(x)*b)
        plt.plot(tx,ty, 'k', linewidth=1.5) # Putting the line on the plot

#
# Include the Kinetic Data
#

        show_data = "Slope = {:0.3e} +/- {:0.1e} /s/mM".format(b, sigb)
        plt.text(0.5, 0.92, show_data, ha='center', va='center', transform=ax.transAxes)
        show_data = "y-intercept = Km/kcat = {:0.3e} +/- {:0.1e} s mM".format(a, siga)
        plt.text(0.5, 0.88, show_data, ha='center', va='center', transform=ax.transAxes)
        show_data = "1/intercept = kcat/Km = {:0.3e} +/- {:0.1e} /s/mM".format(1.0/a, 1.0/a*siga/a)
        plt.text(0.5, 0.84, show_data, ha='center', va='center', transform=ax.transAxes)
        show_data = "Slope/intercept = {:0.3e}".format(b/a)
        plt.text(0.5, 0.80, show_data, ha='center', va='center', transform=ax.transAxes)
#
# SAVE THE PLOT AS A PDF FILE
#
        plfilenm = name[:-4]+"_ExtrapPlot"+str(k)+".pdf"
        plt.savefig(plfilenm, format='pdf')
        plt.show()
        plt.close(plfilenm)
#
# If final absorbances didn't match reset [S] values to
# what would be obtained using the calculated final aborbance 
#
        if final_abs_flag == "true":
            for i in range(loopcount):
                substrate[i] = 1-(starting_abs[i]-absdata[0])/(init_abs2+deltaAf)
#         
print("PROGRAM COMPLETE")
