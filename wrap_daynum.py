#!/usr/bin/env python

import subprocess

def call(cmd):
  print(cmd)
  subprocess.call('%s' % (cmd), shell=True, executable='/bin/csh')

hours = range(0,60,1)

for h in hours:
  #cmd = 'ncl plot_raw_panel.ncl \'daynum=%02d\'' % (int(h))
  #cmd = 'ncl plot_data_panel.ncl \'daynum=%02d\'' % (int(h))
  #cmd = 'ncl plot_mask_panel.ncl \'daynum=%02d\'' % (int(h))
  #cmd = 'ncl plot_prime_panel.ncl \'daynum=%02d\'' % (int(h))
  cmd = 'ncl plot_st_input.ncl \'daynum=%02d\'' % (int(h))
  call(cmd)
