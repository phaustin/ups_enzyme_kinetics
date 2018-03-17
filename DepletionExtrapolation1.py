"""
  to run the script in iteractive mode:

     python DepletionExtrapoloation1.py

  to run the script in batch mode with a default config_file

     python DepletionExtrapoloation1.py  -r

  to run the the script in batch mode with a specified config_file

     python DepletionExtrapoloation1.py  -c path_to_config_file
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist as AA
from scipy.optimize import curve_fit
import copy
from builtins import input
from builtins import str
from builtins import range
import csv
import pdb
from collections import namedtuple
import json
import argparse
import textwrap
import sys
from pathlib import Path
from textwrap import dedent
from ups_enzyme.input_output import (initial_setup,make_plots,final_plot,
                                     get_config_dir, get_output_dir, get_plot_dir)

def make_parser():
    linebreaks = argparse.RawTextHelpFormatter
    descrip = textwrap.dedent(globals()['__doc__'])
    parser = argparse.ArgumentParser(formatter_class=linebreaks,
                                     description=descrip)
    parser.add_argument('-c','--config_file',type=str,help='optional name of json config file')
    parser.add_argument('-r','--run_default',action='store_true',
                           dest='default_run',help='run the default case and create a config file')
    return parser


def main(the_args=None):
    """
    the_args: optional -- if missing then args will be taken from command line
    """
    parser = make_parser()
    args = parser.parse_args(the_args)
    if args.config_file:
        # 
        # config_file argument found, run in batch mode
        #
        print('batch run with config_file={}'.format(args.config_file))
        do_default=False
        interactive=False
        config_file=args.config_file
    else:
        #
        #  if no config_file but default_run requested, will create default config
        #
        if args.default_run:
            interactive=False
            do_default = True
            config_file='default_config.json'
            print(('batch run requested, will create default '
                   'config file '.format(config_file)))
        else:
            #
            # no default request, no config file, run interactive
            #
            text = f"""
                starting interactive run
                to see the help message, do
                python {sys.argv[0]} -h
            """
            print(dedent(text))
            config_file=None
            do_default=False
            interactive=True
    #pdb.set_trace()
    name, steps, enzconc = initial_setup(config_file,do_default)
    input_file = get_config_dir() / Path(name)
    if not input_file.is_file():
        error = (f"cannot find input file, looking for\n"
                 f"{str(input_file)}")
        raise FileNotFoundError(error)
    final_abs_flag = "false"
    different_k_flag = "false"  # Not yet implemented
    timedata = []
    absdata = []
    with open(input_file, 'rU') as data_file:
        items = csv.reader(data_file, delimiter='\t')
        for pair in items:
            convert_pair = [float(value) for value in pair]
            timedata.append(convert_pair[0])
            absdata.append(convert_pair[1])
    timedata = np.array(timedata)
    absdata = np.array(absdata)
    num_data = len(timedata)
    #
    # Write the data out to a file (for testing purposes). I will probably remove
    # this from the final version
    #
    output_filename = get_output_dir() /Path(f'{input_file.stem}_Analysis.txt')
    with open(output_filename, 'a') as outputfile:
        outputfile.write(f'Input File: {str(input_file)}\n')
        outputfile.write(f'Output File: {str(output_filename)}\n')
    #
    if steps < 1:
        steps = 1
    if steps > 10:
        steps = 10
    with open(output_filename, 'a') as outputfile:
        outputfile.write('The number of fits selected is: ' + str(steps) +
                         '\n')
    #
    with open(output_filename, 'a') as outputfile:
        outputfile.write('The enzyme concentration is: ' + str(enzconc) + '\n')
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
    loop_keys=['abs','deltaAf','enzconc','init_abs','init_abs2',
               'current_num_data','k','k2','kcat_Km','kcat_Km2',
               'loopcount','name','pcov','pcov1','substrate','time','vi','input_file']
    loop_tuple= namedtuple('loop_tuple', sorted(loop_keys))
    final_keys=['a','b','kval','name','siga','sigb','x','y','input_file']
    final_tuple=namedtuple('final_tuple',sorted(final_keys))
    
    for loopcount in range(steps):
        #
        # SET UP THE CURRENT TIME AND ABS DATA
        #
        if loopcount == 0:
            current_startabs = absdata[0]
        else:
            current_startabs = absdata[0] + loopcount * (
                deltaA_overall) / 10  # This may not correspond to 10%rxn
        startpt_counter = 0
        while absdata[startpt_counter] < current_startabs:
            startpt_counter = startpt_counter + 1
        current_num_data = num_data - startpt_counter
        print("startpt_counter = ", startpt_counter)
        print("current_num_data = ", current_num_data)
        abs = np.zeros(current_num_data)
        time = np.zeros(current_num_data)
        for i in range(startpt_counter, num_data):
            abs[i - startpt_counter] = absdata[i]
            time[i - startpt_counter] = timedata[
                i]  # I guess it isn't necessary to reset the data to start at time zero?
        with open(output_filename, 'a') as outputfile:
            outputfile.write('\n')
            outputfile.write("Fit # " + str(loopcount + 1) + '\n')
            outputfile.write("The current number of data points is: ")
            outputfile.write(str(current_num_data) + '\n')
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
            initfitpts = int(current_num_data * 0.05)
            #pdb.set_trace()
            for i in range(initfitpts):
                timetwo.append(time[i])
                abstwo.append(abs[i])
            initialfit = np.polyfit(
                timetwo, abstwo, 1)  #Fitting to a line and getting the slope
            vi = initialfit[0] * 1.2
    #
    # For all but the first cycle, it will just be using the vi from the previous run
    #
    # A STUPID, BUT GOOD ENOUGH ESTIMATE FOR k, assuming data span about 5 half-lives
    # For all but the first cycle, it will use the previous k
    #
        if loopcount == 0:
            k = 5. / (time[current_num_data - 1] - time[0]
                      )  # I know this is stupid, but it works.
        init_abs = abs[0]  # This seems obvious
        initguess = np.zeros(3)
        initguess[0] = abs[0]
        initguess[1] = vi
        initguess[2] = k  #This could be changed to k2 if using second routine
        #
        # FIT TO CURVE - Method 1
        # Here's where the magic happens!
        #
        def func(time, init_abs, vi, k):
            return init_abs + (vi) * (1 - np.exp(-1. * k * time)) / k
        popt, pcov = curve_fit(func, time, abs, p0=initguess)
        #
        pcov1 = copy.copy(
            pcov
        )  # Make a copy to have available for printing out StdDev on plot
        init_abs = popt[0]
        vi = popt[1]
        k = popt[2]
        kcat_Km[loopcount] = k / enzconc
        #
        #  Initially assume that the deltaA_overall corresponds to the average of the
        #  last 5% of datapoints - the starting absorbance. This is used to estimate
        #  where the 10% reaction points are for the extrapolation.
        #
        if loopcount == 0:
            endstart = int(0.95 * (num_data - 1))
            counter = 0
            sum = 0
            for i in range(endstart, num_data - 1):
                counter = counter + 1
                sum = sum + abs[i]
            deltaA_overall = sum / counter - abs[0]
            print("deltaA_overall = ", deltaA_overall)
        substrate[loopcount] = 1 - (abs[0] - absdata[0]) / deltaA_overall
        starting_abs[loopcount] = abs[0]
        stddev[loopcount] = kcat_Km[loopcount] * np.sqrt(pcov[2, 2]) / k
        print(substrate[loopcount], kcat_Km[loopcount], stddev[loopcount])
        print()
        print('[Substrate] = ', substrate[loopcount])
        print('Initial Abs = ', init_abs, ' +/- ', np.sqrt(pcov[0, 0]))
        print('Initial Rate = ', vi, ' +/- ', np.sqrt(pcov[1, 1]))
        print('Rate constant = ', k, ' +/- ', np.sqrt(pcov[2, 2]))
        print('[Enzyme] = ', enzconc)
        print('kcat/Km = ', kcat_Km[loopcount], ' +/- ',
              kcat_Km[loopcount] * np.sqrt(pcov[2, 2]) / k)
        #
        with open(output_filename, 'a') as outputfile:
            outputfile.write('[Relative S] = ' + str(substrate[loopcount]) +
                             '\n')
            outputfile.write("Initial Abs = " + str(init_abs) + ' +/- ' +
                             str(np.sqrt(pcov[0, 0])) + '\n')
            outputfile.write('Initial Rate = ' + str(vi) + ' +/- ' +
                             str(np.sqrt(pcov[1, 1])) + '\n')
            outputfile.write('Rate constant = ' + str(k) + ' +/- ' +
                             str(np.sqrt(pcov[2, 2])) + '\n')
            outputfile.write('[Enzyme] = ' + str(enzconc) + '\n')
            outputfile.write(
                'kcat/Km = ' + str(kcat_Km[loopcount]) + ' +/- ' + str(
                    kcat_Km[loopcount] * np.sqrt(pcov[2, 2]) / k) + '/s/M \n')
    #
    #
    # ESTIMATING PARAMETERS FOR CURVE FIT - Method 2
    #
        if loopcount == 0:
            deltaAf = abs[current_num_data -
                          1] - abs[0]  #after the first time around just use the previous estimate
        initguess2 = np.zeros(3)
        initguess2[
            0] = init_abs  #Using the init_abs from the Method 1 curve fit should be a good guess
        initguess2[1] = deltaAf  #Just use the value determined last tiome
        initguess2[
            2] = k  # Using the k from the Method 1 curve fit should be a good guess!
        #
        # FIT TO CURVE - Method 2
        #
        def func(time, init_abs2, deltaAf, k2):
            return init_abs2 + (deltaAf) * (1 - np.exp(-1. * k2 * time))
        popt, pcov = curve_fit(func, time, abs, p0=initguess2)
        init_abs2 = popt[0]
        deltaAf = popt[1]
        k2 = popt[2]
        kcat_Km2[loopcount] = k2 / enzconc
        #
        #   substrate[loopcount] is the same as already determined
        #
        stddev2[loopcount] = kcat_Km2[loopcount] * np.sqrt(pcov[2, 2]) / k2
        #
        print(substrate[loopcount], kcat_Km2[loopcount], stddev2[loopcount])
        print()
        print('Initial Abs = ', init_abs2, ' +/- ', np.sqrt(pcov[0, 0]))
        print('Delta Abs = ', deltaAf, ' +/- ', np.sqrt(pcov[1, 1]))
        print('Final Abs = ', init_abs2 + deltaAf)
        print('Rate constant = ', k2, ' +/- ', np.sqrt(pcov[2, 2]))
        print('kcat/Km = ', kcat_Km2[loopcount], ' +/- ',
              kcat_Km2[loopcount] * np.sqrt(pcov[2, 2]) / k2)
        #
        with open(output_filename, 'a') as outputfile:
            outputfile.write('\n')
            outputfile.write('[Relative S] = ' + str(substrate[loopcount]) +
                             '\n')
            outputfile.write("Initial Abs = " + str(init_abs2) + ' +/- ' +
                             str(np.sqrt(pcov[0, 0])) + '\n')
            outputfile.write('Delta Abs = ' + str(deltaAf) + ' +/- ' +
                             str(np.sqrt(pcov[1, 1])) + '\n')
            outputfile.write('Final Abs = ' + str(init_abs2 + deltaAf) + '\n')
            outputfile.write('Rate constant = ' + str(k2) + ' +/- ' +
                             str(np.sqrt(pcov[2, 2])) + '\n')
            outputfile.write('[Enzyme] = ' + str(enzconc) + '\n')
            outputfile.write(
                'kcat/Km = ' + str(kcat_Km2[loopcount]) + ' +/- ' + str(
                    kcat_Km2[loopcount] * np.sqrt(pcov[2, 2]) / k) + '/s/M \n')
        keep_dict={}
        plot_dir = get_plot_dir()
        current_dict=locals()
        for key_val in loop_keys:
            keep_dict[key_val]=current_dict[key_val]
        plot_data=loop_tuple(**keep_dict)
        plfilenm=make_plots(plot_data)
        with open(output_filename, 'a') as outputfile:
            outputfile.write(f'Plot saved as {plfilenm}\n\n')
    #
    # AFTER GOING THROUGH AS MANY CYCLES AS DESIRED, WE NOW DROP DOWN TO DO THE EXTRAPOLATION
    # OR JUST DROP OUT IF ONLY ONE CYCLE
    #
    #
    # Check to see that the estimate of the Final Abs from the last fit, matches
    # the final absorbance calculated from the last fit. If these differ by more
    # than 1% give a warning.
    #
    loopcount = loopcount + 1
    if loopcount > 1:
        print("Final abs observed = ", deltaA_overall + absdata[0],
              " Final abs calcd = ", init_abs2 + deltaAf)
        print("Calcd final abs/observed final abs = ",
              (init_abs2 + deltaAf) / (deltaA_overall + absdata[0]))
        with open(output_filename, 'a') as outputfile:
            outputfile.write('\n')
            outputfile.write('Calcd final abs/observed final abs = ' + str(
                (init_abs2 + deltaAf) / (deltaA_overall + absdata[0])) + '\n')
        if (init_abs2 + deltaAf) / (deltaA_overall + absdata[0]) > 1.001:
            print(
                "WARNING: the final absorbance calculated differs from the final absorbance observed"
            )
            final_abs_flag = "true"
            with open(output_filename, 'a') as outputfile:
                outputfile.write(
                    'WARNING: the final absorbance calculated differs from the final absorbance observed'
                    + '\n')
        if (init_abs2 + deltaAf) / (deltaA_overall + absdata[0]) < 0.999:
            print(
                "WARNING: the final absorbance calculated differs from the final absorbance observed"
            )
            final_abs_flag = "true"
            with open(output_filename, 'a') as outputfile:
                outputfile.write(
                    'WARNING: the final absorbance calculated differs from the final absorbance observed'
                    + '\n')
        if interactive:
            nothing = input(
                "Just pausing so you can see the ratio, push return to continue")
        if final_abs_flag == "true":
            extrap_plot = 2
        else:
            extrap_plot = 1
        for kval in range(extrap_plot):
            xvalues = []
            yvalues = []
            svalues = []
            for i in range(loopcount):
                xvalues.append(substrate[i])
                yvalues.append(kcat_Km[i])
                svalues.append(stddev[i])
                print(xvalues[i], yvalues[i])
            extrap_kcat_km = np.polyfit(xvalues, yvalues, 1)
            print('Extrapolated kcat/Km without weighting = ',
                  extrap_kcat_km[1])
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
            sx = 0
            sy = 0
            st2 = 0
            b = 0
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
                y[i] = 1.0 / (kcat_Km[i] / 1000)
                s[i] = stddev[i] / kcat_Km[i] * y[i]
                if s[i] > 0:
                    w[i] = 1. / s[i]**2
                else:
                    w[i] = 1.
                print("{:0.4e} {:0.4e} {:0.4e} {:0.4e}".format(
                    x[i], y[i], s[i], w[i]))
                ss = ss + w[i]
                sx = sx + x[i] * w[i]
                sy = sy + y[i] * w[i]
    #
            sxoss = sx / ss
            #
            for i in range(loopcount):
                t = (x[i] - sxoss) / (1. / np.sqrt(w[i]))
                st2 = st2 + t * t
                b = b + t * y[i] / (1. / np.sqrt(w[i]))
    #
    # SOLVE FOR SLOPE(b), INTERCEPT(a), AND SIGMAS (siga, sigb)
    #
            b = b / st2
            a = (sy - sx * b) / ss
            siga = np.sqrt((1. + sx * sx / (ss * st2)) / ss)
            sigb = np.sqrt(1. / st2)
            #
            #
            # OUTPUT RESULTS
            #
            # Convert kcat/Km to /s/mM from /s/M
            #
            #    print
            print("Slope = {:0.3e} +/- {:0.1e} /s/mM".format(b, sigb))
            print("Intercept = {:0.3e} +/- {:0.1e} s•mM".format(a, siga))
            print("1/Intercept = kcat/Km = {:0.3e} +/- {:0.1e} /s/mM".format(
                1.0 / a, 1.0 / a * siga / a))
            print("Slope/intercept = {:0.3e}".format(b / a))
            #
            keep_dict={}
            current_dict=locals()
            for key_val in final_keys:
                keep_dict[key_val]=current_dict[key_val]
            plot_data=final_tuple(**keep_dict)
            filename=final_plot(plot_data)
            print('creating {}'.format(filename))
            with open(output_filename, 'a') as outputfile:
                outputfile.write('\n')
                if kval == 0:
                    outputfile.write(
                        'Extrapolating based on last 5% of data = Final Absorbance'
                        + '\n')
                else:
                    outputfile.write(
                        'Extrapolating based on fit Final Absorbance' + '\n')
                outputfile.write('Slope = ' + str(b) + ' +/- ' + str(sigb) +
                                 '\n')
                outputfile.write('Intercept = ' + str(a) + ' +/- ' + str(siga)
                                 + '\n')
                outputfile.write('1/Intercept = kcat/Km ' +
                                 str(1.0 / a) + ' +/- ' + str(
                                     1.0 / a * siga / a) + ' /mM/s' + '\n')
                outputfile.write('Slope/Intercept = ' + str(b / a) + '\n')
            
            #
            # If final absorbances didn't match reset [S] values to
            # what would be obtained using the calculated final aborbance
            #
            if final_abs_flag == "true":
                for i in range(loopcount):
                    substrate[i] = 1 - (starting_abs[i] - absdata[0]) / (init_abs2 + deltaAf - absdata[0])
    #
    print(f'final output written to: {str(output_filename)}')
    print("PROGRAM COMPLETE")

if __name__== "__main__":
    main()
