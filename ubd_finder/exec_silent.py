#!/usr/bin/python

import commands
import os

def exec_silent(command):
  os.system(command)
  return
  status, output = commands.getstatusoutput(command)

  if status != 0:
    print 'Command "%s" failed with status %s' % (command, status)
    print output

  
