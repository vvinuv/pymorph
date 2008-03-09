#!/usr/bin/env python
import os
import compileall
compileall.compile_dir('.', force=True)
os.system('chmod 755 pymorph.pyc')
os.system('chmod 755 pymorph.py')
