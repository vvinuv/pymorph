#++++++++++++++++++++++++++
#
# TITLE: cmp_results
#
# PURPOSE: Compare the Results of the
#          postage stamp tests to the
#          expected fits produced
#          in the Meert Catalog
#
# INPUTS: Reads the result.csv file from
#         the postage stamp results.
#
# OUTPUTS: Plots comparing some fit parameters
#
# PROGRAM CALLS: numpy, pylab
#
# BY: Alan Meert
#     Department of Physics and Astronomy
#     University of Pennsylvania
#
# DATE: 7 JAN 2014
#
#-----------------------------------

import numpy as np
import pylab as pl

def read_col_names(filename):
    """Reads the output names of fitted parameters from the result.csv file"""
    infile = open(filename)
    nameline = infile.readline()
    infile.close()

    nameline = nameline.strip() #remove linebreaks and whitespace
    cols = nameline.replace('#','').split(',')
    col_names = dict([('_'.join(a.split('_')[0:-1]), int(a.split('_')[-1])) for a in cols])
    
    return col_names

    
if __name__ == "__main__":

    ex_res_file = './expected_results/result.csv'
    new_res_file = './results/result.csv'
    cols_to_compare = ['C', 'n', 're_pix', 'z']

    plots_square = int(np.ceil(np.sqrt(len(cols_to_compare))))

    ex_res_names = read_col_names(ex_res_file)
    new_res_names = read_col_names(new_res_file)

    new_col_nums = [ new_res_names[a]-1 for a in  cols_to_compare]
    ex_col_nums = [ ex_res_names[a]-1 for a in  cols_to_compare]
    
    new_data = np.loadtxt(new_res_file, usecols = new_col_nums, unpack=True, delimiter = ',', skiprows = 1)
    ex_data = np.loadtxt(ex_res_file, usecols = ex_col_nums, unpack=True, delimiter = ',', skiprows = 1)

    
    for count, (plotname, ex_val, new_val) in enumerate(zip(cols_to_compare,ex_data,  new_data)):
        pl.subplot(plots_square, plots_square, count+1)
        pl.scatter(ex_val, new_val)
        pl.xlabel('expected %s' %plotname)
        pl.ylabel('new %s' %plotname)
        pl.title('%s comparison' %plotname)

    pl.subplots_adjust(hspace = 0.5, wspace = 0.5)
    pl.savefig('test_cmp.pdf')
    
