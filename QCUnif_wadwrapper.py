#!/usr/bin/env python
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# PyWAD is open-source software and consists of a set of plugins written in python for the WAD-Software medical physics quality control software. 
# The WAD Software can be found on https://github.com/wadqc
# 
# The pywad package includes plugins for the automated analysis of QC images for various imaging modalities. 
# PyWAD has been originaly initiated by Dennis Dickerscheid (AZN), Arnold Schilham (UMCU), Rob van Rooij (UMCU) and Tim de Wit (AMC) 
#
#
# Changelog:
#
#
# Description of this plugin:
# This plugin calculates the NEMA differential and integral uniformity of an image 
#


__version__ = '20180714'
__author__ = 'DD, tdw'


try:
    import pydicom as dicom
    from pydicom import tag
except ImportError:
    import dicom
    from  dicom import tag

import sys,os
import getopt
#from dicom import tag
import numpy as np
from numpy import random as rnd

import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

# this will fail unless wad_qc is already installed
from wad_qc.module import pyWADinput
from wad_qc.modulelibs import wadwrapper_lib
import NEMA_unif_lib as nemalib

import os



class uniformity_nm:
    def __init__(self, data, results, action):

        self.params = action.get('params', {})
        self.detector_identification = self.params.get('detector_identification',{})
        self.dome_corr = self.params.get('perform_dome_correction',0)
        self.dome_corr = self.dome_corr in ['yes','ja','1','true',1]

        print("Perform dome correction = {}".format(self.dome_corr))

        # data comes either as:
        #   - two series of each 1 instance (1 instance for each detector)
        #   - one series, one instance (multiframe: frame0=detector1, frame1=detector2)
        #   - one series, two instances (1 instance for each detector)
        serieslist = data.getAllSeries()
        det1=det2=None

        if len(serieslist)==0:            # #series=0
            print("Error: no data present")
            exit(False)
        elif len(serieslist)==1:          # #series=1
            if len(serieslist[0])==0:
                print("Error: no data present")
                exit(False)
            elif len(serieslist[0])==1:   # #series=1 #instances=1 so either multiframe or single detector
                dicomobject=serieslist[0][0]
                frames=int(dicomobject[tag.Tag("0028","0008")].value)
                tmppixelmap=self.read_dicom(dicomobject)

                if frames>1:              # multiframe
                    det1 = tmppixelmap[0]
                    det2 = tmppixelmap[1]
                else:                     # singe-frame single series (single detector)
                    # check which detector matches this single-frame instance
                    if not self.detector_identification:
                        print("no specification: assume det1=instance1")
                        det1=tmppixelmap[0]
                    else:
                        if self.match(dicomobject,"det1"):
                            det1 = tmppixelmap[0]
                        elif self.match(dicomobject,"det2"):
                            det2 = tmppixelmap[0]
                        else:
                            print("Data does not match search detector_identification criteria")

            elif len(serieslist[0])==2:  # #series=1 #instances=2
                # check if each instance is single-frame
                if int(serieslist[0][0][tag.Tag("0028","0008")].value)>1 or int(serieslist[0][1][tag.Tag("0028","0008")].value)>1:
                    print("Error: multi-frame not supported for series with multiple instances")
                    exit(False)

                # check which instance belongs to which detector
                if not self.detector_identification:
                    print("no specification: assume det1=instance1, det2=instance2")
                    det1=self.read_dicom(serieslist[0][0])[0]
                    det2=self.read_dicom(serieslist[0][1])[0]
                else:
                    for i in serieslist[0]:
                        if self.match(i,"det1"):
                            det1 = self.read_dicom(i)[0]
                        if self.match(i,"det2"):
                            det2 = self.read_dicom(i)[0]

            else:
                print("Error: more than 2 detectors currently not supported")
                exit(False)
        elif len(serieslist)==2:          # #series=2
            if all(len(s)==1 for s in serieslist):  # max 1 instance/series supported
                print("#series: {}".format(len(serieslist)))

                # check if each instance is single-frame
                if int(serieslist[0][0][tag.Tag("0028","0008")].value) > 1:
                    print("Error: multi-frame not supported for multiple series")
                    exit(False)
                # check which series belongs to which detector
                if not self.detector_identification:
                    print("no specification: assume det1=series1, det2=series2")
                    det1=self.read_dicom(serieslist[0][0])[0]
                    det2=self.read_dicom(serieslist[1][0])[0]
                else:
                    for s in serieslist:
                        if self.match(s[0],"det1"):
                            det1 = self.read_dicom(s[0])[0]
                        if self.match(s[0],"det2"):
                            det2 = self.read_dicom(s[0])[0]

            else:
                print("Error: only one instance per series supported for multiple series")
                exit(False)
        else:
            print("Error: more than 2 series per study currently not supported")
            exit(False)

        # at this point det1 and/or det2 are correctly identified, so perform calculation

        if det1 is not None:
            print("Performing uniformity calculation for detector 1...")
            self.Uniformity_main(det1,results,'det1',self.dome_corr)
        if det2 is not None:
            print("Performing uniformity calculation for detector 2...")
            self.Uniformity_main(det2,results,'det2',self.dome_corr)

        print("End of class init")


    def match(self,dicomobject,detectorname):
        # if no matching criteria are provided return False
        found_match = (len(self.detector_identification[detectorname])>0)
        for k,v in self.detector_identification[detectorname].items():
            tmp = k.split(',')
            try:
                tmptag = tag.Tag(tmp[0],tmp[1])
                found_match = found_match and dicomobject[tmptag].value==str(v)
            except KeyError:
                print("Dicom tag ({}) not present".format(tmptag))
                found_match = False
        return found_match


    def read_dicom(self,dicomobject):
       try:
            bits = int(dicomobject[tag.Tag("0028","0100")].value)
            frames = int(dicomobject[tag.Tag("0028","0008")].value) 
            rows = int(dicomobject[tag.Tag("0028","0010")].value)
            columns = int(dicomobject[tag.Tag("0028","0011")].value)

            count = rows*columns*frames

            if bits == 8:
                dtype = np.uint8
            elif bits == 16:
                dtype = np.uint16
            elif bits == 32:
                dtype = np.uint32
       except:
            print("Error reading dicom image")
            exit(False)

       return np.reshape(np.fromstring(dicomobject.PixelData,count=count,dtype=dtype),(frames,rows,columns))


    def check_dicomdict_studies(self,dicomdict):
        if len(dicomdict) > 1:
            print("dicomdict contains more than 1 study! aborting!")
            sys.exit()


    def savewadfig(self,impath,pixelmap):

        try:
            plt.imsave(impath,pixelmap)

        except:
            print("could not save image")

        return


    def Uniformity_main(self, pixelmap, results, label, dome_corr):

        if pixelmap is None:
            return False

        print('Starting Uniformity_main')
        detector = label
        print('Calculating results')
        print('pixelmap size: {}'.format(np.shape(pixelmap)))
        
        output = nemalib.calculate_nema_uniformity(pixelmap, (64,64), results, dome_corr)
        
        '''
          output[0] = DUx in UFOV
          output[1] = DUy in UFOV
          output[2] = coordinates max(DUx) in UFOV
          output[3] = coordinates max(DUy) in UFOV

          output[4] = DUx in CFOV
          output[5] = DUy in CFOV
          output[6] = coordinates max(DUx) in CFOV
          output[7] = coordinates max(DUy) in CFOV

          output[8] = IU in UFOV
          output[9] = IU in CFOV
          output[10] = ufov pixelmap (numpy masked array)
          output[11] = cfov pixelmap (numpy masked array)
        '''

        results.addFloat('{}_DU_x (UFOV)'.format(detector), output[0])
        results.addFloat('{}_DU_y (UFOV)'.format(detector), output[1])
        results.addFloat('{}_DU_x (CFOV)'.format(detector), output[4])
        results.addFloat('{}_DU_y (CFOV)'.format(detector), output[5])

        results.addFloat('{}_IU (UFOV)'.format(detector), output[8])
        results.addFloat('{}_IU (CFOV)'.format(detector), output[9])
        
        print("writing image results")

        try:
           filename = 'ufov_{}.png'.format(detector)
           nemalib.save_imgmap(output[10],output[2],output[3],filename)
           results.addObject('UFOV image {}'.format(detector),filename)
        except:
           print("Could not save {}".format(filename))

        try:
           filename = 'cfov_{}.png'.format(detector)
           nemalib.save_imgmap(output[11],output[6],output[7],filename)
           results.addObject('CFOV image {}'.format(detector),filename)
        except:
           print("Could not save {}".format(filename))



if __name__ == "__main__":
    data, results, config = pyWADinput()

    # read runtime parameters for module
    for name,action in config['actions'].items():
        #if name == 'acqdatetime':
        #    acqdatetime_series(data, results, action)
        if name == 'qc_series':
            uniformity_nm(data, results, action)

    results.write()
