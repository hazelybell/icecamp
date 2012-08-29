#!/usr/bin/python
# -*- coding: utf-8 -*-
#    Test Harness for l1l2.F90
#    Copyright (C) Joshua Charles Campbell 2012

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import fileinput
import re
from numpy import *
import matplotlib.pyplot

var = 0
x = []
y = []
v = []
t = []
reading = ''


for line in fileinput.input():
  line = line.rstrip()
  m = re.search('Variable\s+(\d+)', line, flags=re.IGNORECASE)
  trianglestart = re.search('Triangles', line, flags=re.IGNORECASE)
  threefloats = re.search('^\s*([\d.E+-]+)\s*([\d.E+-]+)\s*([\d.E+-]+)', line)
  threeints = re.search('^\s*(\d+)\s*(\d+)\s*(\d+)', line)
  endof = re.search('End', line, flags=re.IGNORECASE)
  #print line
  if (endof):
    print "END END"
    if (reading == 't'):
      t = asarray(t)
      #print t
    elif (reading == 'v'):
      x = asarray(x)
      y = asarray(y)
      v = asarray(v)
      #print x
      #print y
      #print v
      matplotlib.pyplot.figure(var)
      matplotlib.pyplot.tripcolor(x, y, t, v)
      matplotlib.pyplot.colorbar()
      if (var <= 5):
        contours = matplotlib.pyplot.tricontour(x, y, t, v, colors='k') 
        matplotlib.pyplot.clabel(contours, fmt='%1.0f')
    reading = ''
  elif (m):
    x = []
    y = []
    v = []
    var = int(m.group(1))
    print "Reading var %i" % var
    reading = 'v'
  elif (trianglestart):
    print "Reading Triangles"
    t = []
    reading = 't'
  elif (reading == 'v' and threefloats):
    x.append(float(threefloats.group(1)))
    y.append(float(threefloats.group(2)))
    v.append(float(threefloats.group(3)))
  elif (reading == 't' and threeints):
    t.append([int(threefloats.group(1)), int(threefloats.group(2)), int(threefloats.group(3))])
  else:
    print line

matplotlib.pyplot.show()