#!/usr/bin/env python

import os, re, string, subprocess, sys, importlib
Import('env')
sys.path.append(os.getenv("MUSE_WORK_DIR")+'/site_scons')
#------------------------------------------------------------------------------
# print("Stntuple/SConscript:muse branch: PWD:"+os.getenv("PWD"))

# pass the package name as a parameter

x = subprocess.call(os.getenv("MUSE_WORK_DIR")+'/Stntuple/scripts/build_config_muse tvst02',shell=True)
# print("Stntuple/SConscript back from build_config_muse")

tvst02_env = env.Clone()
#------------------------------------------------------------------------------
# done
#------------------------------------------------------------------------------
exec(open(os.environ['MUSE_WORK_DIR']+"/site_scons/stntuple_site_init.py").read())

from stntuple_helper    import *

tvst02_env.Append(BUILDERS = {'StntupleCodegen'  : stntuple_codegen })
tvst02_env.Append(BUILDERS = {'StntupleRootCint' : stntuple_rootcint})

tvst02_env['CPPPATH' ].append(os.environ['MUSE_WORK_DIR']+'/include');

tvst02_env.Append(FORTRANPATH = [os.environ['MUSE_WORK_DIR']+'/include']);

# print(tvst02_env.Dump())

Export('tvst02_env')
Export('stntuple_helper')
