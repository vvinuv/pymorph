#!/usr/bin/env python
import os
import compileall
compileall.compile_dir('/home/vinu/serial_pipeline/trunk/', force=True)
os.system('chmod 755 pymorph.pyc')
