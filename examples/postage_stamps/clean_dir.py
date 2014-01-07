#!/data2/home/ameert/python/bin/python2.5

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

