#!/usr/bin/python

# this code is written for python >= 3.6 and require matplotlib and numpy
#
# usage: python cadi24h.py sourcedir outputfile
# sourcedir = directory name of the CADI data (one day)
# outputfile = file name of the resulting plot
#
# outputfile should have extension .png to produce png-files or .pdf to produce pdf-files.
# supported formats eps, pdf, pgf, png, ps, raw, rgba, svg, svgz (depends on version of matplotlib)
#
# the code is written based on IDL code provided by Chris Meek
# IDL code originally written by Ian Grant and modified by the same and JWM
#

import glob, os
import sys
import struct
import datetime
from time import gmtime,strptime
import numpy as np
import matplotlib
import copy
#matplotlib.use('Agg')        # use the 'Agg' backend when $DISPLAY environment is not defined
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (10, 4)

max_ntimes = 256
max_ndopbins = 300000
nantennas = 4
dheight = 3.0
height=[]
ff=[]
timeaxis=[]

for filename in glob.glob(sys.argv[1]+"/*.md49"):     # change file extention of CADI data if necessary
   f = open(filename, "rb")
   try:
      # 1) read header information as described in the documentation p. 26-27
      site = f.read(3).decode("utf-8")
      ascii_datetime = f.read(22).decode("utf-8")
      filetype = f.read(1).decode("utf-8")
      nfreqs = struct.unpack("<H", f.read(2))[0]
      ndops = struct.unpack("<B", f.read(1))[0]
      minheight = struct.unpack("<H", f.read(2))[0]
      maxheight = struct.unpack("<H", f.read(2))[0]
      pps = struct.unpack("<B", f.read(1))[0]
      npulses_avgd = struct.unpack("<B", f.read(1))[0]
      base_thr100 = struct.unpack("<H", f.read(2))[0]
      noise_thr100 = struct.unpack("<H", f.read(2))[0]
      min_dop_forsave = struct.unpack("<B", f.read(1))[0]
      dtime = struct.unpack("<H", f.read(2))[0]
      gain_control = f.read(1).decode("utf-8")
      sig_process = f.read(1).decode("utf-8")
      noofreceivers = struct.unpack("<B", f.read(1))[0]
      spares = f.read(11).decode("utf-8")

      month = ascii_datetime[1:4]
      day = int(ascii_datetime[5:7])
      hour = int(ascii_datetime[8:10])
      minute = int(ascii_datetime[11:13])
      sec = int(ascii_datetime[14:16])
      year = int(ascii_datetime[17:21])

      month_number = strptime(month, '%b').tm_mon
      mydate = datetime.date(year, month_number, day)
      jd = mydate.toordinal() + 1721424.5
      jd0jd = datetime.date(1986, 1, 1)
      jd0 = jd0jd.toordinal() + 1721424.5

      time_header = (jd - jd0)*86400 + hour*3600 + minute*60 + sec
      time_hour = 3600 * (time_header/3600)

      # 2) read all frequencies used

      freqs = [struct.unpack("<f", f.read(4))[0] for i in range(nfreqs)]

      if filetype == 'I':
         max_nfrebins = nfreqs
      else:
         max_nfrebins = min(max_ntimes * nfreqs, max_ndopbins)

      # 3) read rawdata

      nheights = int(maxheight / dheight + 1)

      times = []
      frebins = []
      frebins_x = []
      frebins_gain_flag = []
      frebins_noise_flag = []
      frebins_noise_power10 = []
      time_min = 0
      time_sec = 0
      timex = -1
      freqx = nfreqs - 1
      dopbinx = -1
      frebinx = -1
      iq_bytes = np.zeros((noofreceivers,2))
      dopbin_x_timex = []
      dopbin_x_freqx = []
      dopbin_x_hflag = []
      dopbin_x_dop_flag = []
      dopbin_iq = []
      ht_limitmin = 90
      ht_limitmax = 500
      hflag = 0

      time_min = struct.unpack("<B", f.read(1))[0]
      while time_min != 255:
         time_sec = struct.unpack("<B", f.read(1))[0]
         flag = struct.unpack("<B", f.read(1))[0]  # gainflag
         timex += 1
         times.append(time_hour + 60 * time_min + time_sec)
         for freqx in range(nfreqs):
            noise_flag = struct.unpack("<B", f.read(1))[0] # noiseflag
            noise_power10 = struct.unpack("<H", f.read(2))[0]
            frebinx += 1
            frebins_gain_flag.append(flag)
            frebins_noise_flag.append(noise_flag)
            frebins_noise_power10.append(noise_power10)
            htlimitmin = ht_limitmin / 3
            htlimitmax = ht_limitmax / 3
            flag = struct.unpack("<B", f.read(1))[0]
            while flag < 224:
               ndops_oneh = struct.unpack("<B", f.read(1))[0]
               hflag = flag
               if ndops_oneh >= 128:
                  ndops_oneh = ndops_oneh - 128
                  hflag = hflag + 200
               for dopx in range(ndops_oneh):
                  dop_flag = struct.unpack("<B", f.read(1))[0]
                  for rec in range(noofreceivers):
                     iq_bytes[rec,0] = struct.unpack("<B", f.read(1))[0]
                     iq_bytes[rec,1] = struct.unpack("<B", f.read(1))[0]
                  dopbinx += 1
                  dopbin_iq.append(copy.deepcopy(iq_bytes))
                  dopbin_x_timex.append(timex)
                  dopbin_x_freqx.append(freqx)
                  dopbin_x_hflag.append(hflag)
                  dopbin_x_dop_flag.append(dop_flag)
               flag = struct.unpack("<B", f.read(1))[0] # next hflag/gainflag/FF
         time_min = flag # next record
   finally:
      f.close()

   # frebins_gain_flag[] is the gainflag for each frequency
   # frebins_noise_flag[] is the noiseflag for each frequency
   # frebins_noise_power10[] is the 10 x averagenoisepower for each frequency

   # dopbin_iq[] is an array of real and imaginary components for the receivers

   # frequency in MHz
   frequency=[]
   for i in range(len(dopbin_x_freqx)):
       frequency.append(freqs[dopbin_x_freqx[i]]/1000000.0)

   # --- noise reduction - remove lowest and highest frequency for each height
   for index in range(320):
      vh = 30 + index
      collection = [i for i, val in enumerate(dopbin_x_hflag) if val == vh]

      if len(collection) == 1:
         del dopbin_x_hflag[collection[0]]
         del frequency[collection[0]]

      if len(collection) == 2:
         del dopbin_x_hflag[collection[0]]
         del frequency[collection[0]]
         del dopbin_x_hflag[collection[1]-1]
         del frequency[collection[1]-1]

      minimum = 11.
      maximum = 1.
      if len(collection) > 2:
         for j in range(len(collection)):
            if frequency[collection[j]] < minimum:
               minindex = collection[j]
               minimum = frequency[collection[j]]
            if frequency[collection[j]] > maximum:
               maxindex = collection[j]
               maximum = frequency[collection[j]]

         if minindex < maxindex:
            del dopbin_x_hflag[minindex]
            del frequency[minindex]
            del dopbin_x_hflag[maxindex-1]
            del frequency[maxindex-1]
         else:
            del dopbin_x_hflag[maxindex]
            del frequency[maxindex]
            del dopbin_x_hflag[minindex-1]
            del frequency[minindex-1]
   # --- end noise reduction

   height.append(list(np.array(dopbin_x_hflag) * 3))
   hhh = list(np.array(dopbin_x_hflag) * 3)
   ff.append(frequency)
   ttt = [hour+minute/60.] * len(hhh)
   timeaxis.append([hour+(minute+7.5)/60.] * len(hhh))


x = list(matplotlib.cbook.flatten(timeaxis))
y = list(matplotlib.cbook.flatten(height))
col = list(matplotlib.cbook.flatten(ff))

title = 'CADI ' + ascii_datetime[1:8] + ascii_datetime[17:21]

fig, ax = plt.subplots()

cm = plt.cm.get_cmap('jet')
sc = plt.scatter(x, y, c=col, marker="_", s=22, cmap=cm)
ax.grid(True, which='both')
ax.set_xlim(0,24)                    # set X limits (0h to 24h)
ax.set_ylim(0,800)                   # set Y limits (min and max height in km)
ax.set_title(title)
ax.set_xlabel('UT TIME')
ax.set_ylabel('Height (km)')
cbar = plt.colorbar(sc)
cbar.set_label('Frequency (MHz)')

# toggle output to file, comment/uncomment these lines 
# if you don't want output to file: usage: python cadi24h.py sourcedir
savefilename = sys.argv[2]
fig.savefig(savefilename)

# toggle output to screen, comment/uncomment the next line (note you can have both!)
#plt.show()

