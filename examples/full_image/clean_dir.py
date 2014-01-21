#!/usr/bin/python

import os
import sys

to_remove = ['*.fits','*.html', '*.cat',
             '*.png','*.txt','*.pyc','*.in','*.out',
             '*.sex', '*.Shallow', '*.csv', '*.con',
             '*.log']

targetdir = './results/'
    
thisdir = os.getcwd()

os.chdir(targetdir)
for del_file in to_remove:
    os.system('rm %s' %del_file)

os.chdir(thisdir)

os.system("rm ./data/Wtest_237.fits ./data/Itest_237.fits") 
