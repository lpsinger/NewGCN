import numpy as np
import matplotlib as plt

from check_link import simple_get
from requests import get

import re

from IPython.display import clear_output

import csv
import json
from edit_csv import *

from bs4 import BeautifulSoup, Tag

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import Distance


class SaveGCN:

    with open("List_of_Tel.json") as f:
        _list_of_tel = json.load(f)
                     
    _month = ['Jan.', 'Feb.', 'Mar.', 'Apr.', 'May', 'Jun.', 'Jul.', 'Aug.', 'Sep.', 'Oct.', 'Nov.', 'Dec.']
    _month_f = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']

    def __init__(self, gcnNum, testmode=False):
        self.GCNlink = "https://gcn.gsfc.nasa.gov/gcn3/{}.gcn3".format(gcnNum)
        self.gcnNum = gcnNum
        
        # Check the availability of the given GCN
        gcn = simple_get(self.GCNlink)
        try:
            gcn = gcn.decode('utf8').strip()
            self.gcn = gcn.split("\n")
        except:
            bs = BeautifulSoup(gcn, 'html.parser')
            gcn = bs.prettify()
            self.gcn = gcn.split("\n")
        
        # For testmode
        if testmode:
            raw_html = simple_get(self.GCNlink)
            bs = BeautifulSoup(raw_html, 'html.parser')
            print(bs.prettify())
            return
        
    def __ReadTitle__(self):
        
        # Find the event name
        if 'BALROG' in self.gcn[2]:
            self._report=False
            self._obs=False
            return
        try:
            tempE = re.findall("(GRB|G|S|IceCube)[-|\s]?([0-9]+)([a-zA-Z]+)[:| -]?", self.gcn[2][9:])[0]
            self._report = True
        except:
            try:
                tempE = re.findall("(GRB|G|S|IceCube)[-|\s]?([0-9]+).([0-9]+)", self.gcn[2][9:])[0]
                tempE = (tempE[0], tempE[1], 'A')
                self._report = True
            except:
                self.__errorMessage__()
                self._report = False
                return
        
        # Define the event name
        if tempE[0] == "IceCube":
            self.Event = "IC {}{}".format(*tempE[1:])
        else:
            self.Event = "{} {}{}".format(*tempE)
            
        # Define the event date
        self.Date = "{}".format(tempE[1])
        
        contents = str("")
        for line in self.gcn[5:]:
            contents= contents+" "+str(line)
            
        self._contents = contents

    def __errorMessage__(self, errtype = 1):
        print("="*80)
        if errtype == 1:
            print("This GCN ({}) may not be related to GRB, GW, or neutrinos.".format(self.gcnNum))
            print("Please check the GCN.")
        elif errtype == 2:
            print("The observatory (GCN {}) is not likely to be in our list.".format(self.gcnNum))
            print("Please check the observatory list.")
        print("Title: {}".format(self.gcn[2][8:]))
        print("Link: {}".format(self.GCNlink))
        print("="*80)
        return
                
    def __scrapping__(self, verbose=False):
        
        self.__ReadTitle__()
        self.__Telescope__()

        self.Loc = [['',''],['','']]
        self._errC = ''
        self.Redshift = ['', '']
        self._dtr = {}
        if self._report:
            for tel in self._tel:
                self._dtr[tel] = self.__Detector__(tel)
        else:
            return
        
        if self._obs:
            self.__Localization__(verbose = verbose)
            self.__Trigger__(AG = (self._list_of_tel[tel][self.__Detector__(tel)[0]] in ['radio', 'optical']), verbose=verbose)
            self.__Redshift__(verbose = verbose)
        else:
            self.Trigger = "20{}-{}-{}".format(self.Date[:2], self.Date[2:4], self.Date[4:])
            return
        
    def __Telescope__(self):
        
        # Find telescopes
        self._tel = []
        
        for tel in list(self._list_of_tel.keys()):
            if self.gcn[2].find(tel)>=0:
                if tel == 'Zwicky Transient Facility':
                    self._tel.append('ZTF')
                elif tel == 'Dabancheng':
                    self._tel.append('HMT')
                elif tel == 'Discovery Channel Telescope':
                    self._tel.append('DCT')
                elif tel == 'Sardinia Radio Telescope':
                    self._tel.append('SRT')
                else:
                    self._tel.append(tel)
        
        # Remove follow-up telescopes
        if (np.size(self._tel) > 1) and "LIGO" in self._tel:
            self._tel.remove("LIGO") # for the LIGO follow-up observation

        if (np.size(self._tel) > 1) and "IceCube" in self._tel:
            self._tel.remove("IceCube") # for the IceCube follow-up observation
        
        if (np.size(self._tel) > 1) and "Swift" in self._tel:
            self._tel.remove("Swift") # for the Swift follow-up observation
        
        # Check whether a detector in the title. 
        if np.size(self._tel) == 0:
            for tel, dtrs in zip(list(self._list_of_tel.keys()), self._list_of_tel.values()):
                for dtr in list(dtrs.keys()):
                    if self.gcn[2].find(dtr)>=0:
                        self._tel.append(tel)
        
        # Find single reliable telescope?
        if np.size(self._tel) == 0:    
            self.__errorMessage__(2)
            self._report = False
            self._obs = False
            return
        
    def __Detector__(self, tel):

        # Find detectors
        dtrs = []
        flag=False
        if np.size(list(self._list_of_tel[tel].keys())) == 1:
            dtrs = list(self._list_of_tel[tel].keys())
        else:
            for dtr in list(self._list_of_tel[tel].keys()):
                if dtr in self.gcn[2]:
                    dtrs.append(dtr)
                    flag = True
            if not(flag):
                for dtr in list(self._list_of_tel[tel].keys()):
                    for line in self.gcn[2:]:
                        if (line.find(dtr)>=0) and (dtr not in dtrs):
                            dtrs.append(dtr)
                            flag=True
            if not(flag):
                dtrs = list(self._list_of_tel[tel].keys())[0]

        # Find one or more reliable detectors?

        if np.size(dtrs) == 1:
            dtrs = dtrs
        elif tel == '':
            self.__errorMessage__(2)
            self._report = False
            self._obs = False
            return
        
        # The event is detected by the observatories?
        self._multidtr = np.zeros([2,np.size(dtrs)], dtype = 'bool')
        self._multidtr[1] = True
        self._obs = True
        
        if (np.size(dtrs)>1): # For multiple detectors (usually Swift initial report)
            if len(re.findall(" (No|no|Upper|upper|limit|Retraction|inactive|SAA|Planned) ", self.gcn[2]))>0: # Check the title
                self._multidtr[1] = False
            else:
                for line in self.gcn[5:]: # Check the content
                    for i in range(np.size(dtrs)):
                        if line.find(dtrs[i])>=0:
                            self._multidtr[0][i] = True
                    if (sum(self._multidtr[0])>0) and (len(re.findall("(no.candidates|no.candidate|not.detect|Retraction|no.evidence|no.source|no.credible|No.credible|No.source|not.find.any|not detect.any)", line))>0):
                        self._multidtr[1][sum(self._multidtr[0])-1] = False 

                # for Swift
                if len(re.findall("(no)\s(XRT)?\s?(or)?\s?(UVOT)?\s?(data)", self._contents)) > 0:
                    exception = re.findall("(no)\s(XRT)?\s?(or)?\s?(UVOT)?\s?(data)", self._contents)[0]
                    if "XRT" in exception:
                        self._multidtr[1][1] = False
                    if "UVOT" in exception:
                        self._multidtr[1][2] = False
                if np.size(re.findall("\s(not.detect.any)\s?", self._contents))>0:
                    print("<< Exception >>")
                    self._multidtr[1] = False
                
                self._obs = sum(self._multidtr[1])

        else:
            if np.size(re.findall("\s(No|no|Upper|upper|limit|Retraction|inactive|SAA|Planned)\s?", self.gcn[2]))>0:
                self._obs = False 
            else:
                if np.size(re.findall("(\sno.candidates|\sno.candidate|no.credible|not.detect|not.found|not.detected|No.[a-zA-z\-\s]+.candidate|5\-sigma.upper.limit|not.find.any|no.evidence|no.optical|no.significant|do.not.detect|don't.find)\s?", self._contents))>0:
                    self._obs = False
                if np.size(re.findall("\s?(we.found|We.found|we..found|additional.events|we.obtained|well.detected|additional.[a-zA-Z.events])\s?", self._contents))>0:
                    if np.size(re.findall("\s(we.found.no)\s?", self._contents))>0:
                        print("<< Exception >>")
                        self._obs = False
                    else:
                        self._obs = True

            self._multidtr[1] = self._obs            
        return dtrs
        
    def __Trigger__(self, AG = False, verbose=False):
        
        T0Type = 0
        takeY = self.Date[:2]
        takeM = self.Date[2:4]
        takeD = self.Date[4:]
        
        if not(AG):
            monlist = ''
            for mon in self._month:
                monlist += mon+'|'
            for mon in self._month_f:
                monlist += mon+'|'

            T0Type = 0
            
            # Scrapping the trigger date
            for line, line2 in zip(self.gcn[5:-1], self.gcn[6:]):
                line = line+" "+line2
                if T0Type == 0:
                    tempT = re.findall("(20[0-9][0-9])?[-|\s]([0-9][0-9]|{})[-|,|\s]\s?([0-9][0-9]?)[\s,.]".format(monlist[:-1]), line)
                    if len(tempT)>0:
                        if verbose == 3:
                            print("Trigger date ({}) is adopted from \n'{}'".format(tempT, line))
                        if tempT[0][0] != '':
                            takeY = tempT[0][0][2:]
                        if tempT[0][1] in self._month:
                            takeM = "{:.2f}".format((self._month.index(tempT[0][1])+1)/100.)[-2:]
                        elif tempT[0][1] in self._month_f:
                            takeM = "{:.2f}".format((self._month_f.index(tempT[0][1])+1)/100.)[-2:]
                        else:
                            takeM = tempT[0][1]
                        if float(tempT[0][2]) < 10:
                            takeD = "{:.2f}".format(float(tempT[0][2])/100.)[-2:]
                        else:
                            takeD = tempT[0][2]
                        T0Type = 1
                        break
            
            # Scrapping the trigger time
            for line, line2 in zip(self.gcn[5:-1], self.gcn[6:]):
                line = line+" "+line2
                if T0Type <=1:
                    tempT = re.findall("[At|at]?\s?([0-9]+):([0-9]+):([0-9.]+)", line)
                    if len(tempT)>0:
                        if verbose == 3:
                            print("Trigger time ({}) is adopted from \n'{}'".format(tempT, line))
                        takeT = tempT[0]
                        T0Type = 2
                        break

            # Rearrange the trigger time
            if T0Type <= 1:
                self.Trigger = "20{}-{}-{}".format(takeY, takeM, takeD)
            if T0Type == 2:
                Time = "{}:{}:{}".format(*takeT)
                if Time != self.Loc[1][0]:
                    self.Trigger = "20{}-{}-{} {}".format(takeY, takeM, takeD, Time)
                else:
                    self.Trigger = "20{}-{}-{}".format(takeY, takeM, takeD)
        else:
            self.Trigger = "20{}-{}-{}".format(takeY, takeM, takeD)
            
        if self.Trigger[-1] == '.':
            self.Trigger = self.Trigger[:-1]
    
    def __Localization__(self, verbose=False):
        
        LocType = 0
        tempErr = [-1, 'degrees']

        # Scrapping the localization
        try:
            tempRa = re.findall("(RA|ra|Ra|RA[\sJ20\(\).]+|ra[\sJ20\(\).]+)\s?\s?\s?[:|=]\s?\s?\s?([0-9]+)[:|h|\s]\s?([0-9]+)[:|m|\s]\s?([0-9.]+)", self._contents[self._contents.find("RA"):])[-1][1:]
            LocType = 1
            if verbose==3:
                 print("RA is adopted from \n'{}'".format(self._contents[self._contents.find("RA"):][:20]))
        except:
            try:
                tempRa = re.findall("(RA|ra|Ra|RA[\sJ20\(\).]+|ra[\sJ20\(\).]+)\s?\s?\s?[:|=]\s?\s?([0-9.]+)", self._contents[self._contents.find("RA"):])[-1][1]                   
                LocType = 2
                if verbose==3:
                     print("RA is adopted from \n'{}'".format(self._contents[self._contents.find("RA"):][:20]))
            except:
                if self._tel=='MASTER':
                    try:
                        tempRa = re.findall("([0-9]+)[h|\s]\s?([0-9]+)[m|\s]\s?([0-9.]+)", self._contents)[0]
                        LocType = 1
                        if verbose==3:
                             print("RA is adopted from \n'{}'".format(self._contents[self._contents.find("RA"):][:20]))
                    except:
                        pass
                else:
                    pass
    
        try:
            tempDec = re.findall("(DEC|dec|Dec|DEC[\sJ20\(\).]+|Dec[\sJ20\(\).]+)\s?[:|=]\s?([\-|\+]?)([0-9]+)[:|d|\s|°]\s?([0-9]+)[:|m|'|\s|’]\s?([0-9.]+)", self._contents[self._contents.find("D"):])[-1][1:]
            if verbose==3:
                 print("DEC is adopted from \n'{}'".format(self._contents[self._contents.find("DEC"):][:20]))
        except:
            try:
                tempDec = re.findall("(DEC|dec|Dec|DEC[\sJ20\(\).]+|Dec[\sJ20\(\).]+)\s?[:|=]\s?([\-|\+]?)([0-9.]+)", self._contents[self._contents.find("Dec"):])[-1][1:]
                if verbose==3:
                    print("DEC is adopted from \n'{}'".format(self._contents[self._contents.find("DEC"):][:20]))
            except:
                if self._tel=='MASTER':
                    try:
                        tempDec = re.findall("([\-|\+]?)([0-9]+)[d|\s]\s?([0-9]+)[m|'|\s|’]\s?([0-9.]+)", self._contents)[0]
                        if verbose==3:
                             print("DEC is adopted from \n'{}'".format(self._contents[self._contents.find("DEC"):][:20]))
                    except:
                        pass
                else:
                    pass
    
        if LocType == 0:    
            try:
                tempLoc = re.findall("(RA, Dec|\(J2000\))\s+?([0-9]+)[:|h|\s]\s?([0-9]+)[:|m|\s]\s?([0-9.]+)\s+?([\-|\+]?)([0-9]+)[:|d|\s|°]\s?([0-9]+)[:|m|'|\s|’]\s?([0-9.]+)", self._contents)[-1][1:]
                tempRa = tempLoc[:3]
                tempDec = tempLoc[3:]
                LocType = 1
            except:
                try:
                    tempLoc = re.findall("(RA, Dec|[\sJ20\(\)]+)\s?[:|=]\s+?([0-9.]+)\,\s?([\-|\+]?)\s?([0-9.]+)", self._contents)[-1][1:]
                    tempRa = tempLoc[0]
                    tempDec = tempLoc[1:]
                    LocType = 2
                except:
                    pass

        # Rearrange the localization
        if LocType == 1:
            ra = "{}:{}:{}".format(*tempRa)
            dec = "{}{}:{}:{}".format(*tempDec)
            c = SkyCoord(ra, dec, unit=(u.hourangle,  u.deg))
            self.Loc = [["{:.3f}".format(c.ra.deg), "{:.3f}".format(c.dec.deg)], [ra, dec]]
        elif LocType == 2:
            ra = float("{}".format(tempRa))
            dec = float("{}{}".format(*tempDec))
            c = SkyCoord(ra*u.deg, dec*u.deg)
            self.Loc = [[ra, dec], ["{:.0f}:{:.0f}:{:.2f}".format(*c.ra.hms), "{:.0f}:{:.0f}:{:.2f}".format(c.dec.dms[0], abs(c.dec.dms[1]), abs(c.dec.dms[2]))]]

        # Scrapping the localization error
        for line in self.gcn[5:]:
            try:
                tempErr = re.findall('(uncertainty.of|error.radius.of|\+\/\-|pixel.error|uncertainty)\s?(about)?\s([0-9.]+)\s?(deg|degrees|arcmin|arcseconds|arcsec|arc sec|s|\")', line)[0][2:]
            except:
                try:
                    bitempErr = re.findall("\+([0-9.]+)\/\-([0-9.]+) (deg|degrees|arcmin|arcseconds|arc sec)", line)[0]
                    tempErr = [max(float(tempErr[0]), float(bitempErr[0]), float(bitempErr[1])), bitempErr[2]]
                except:
                    try:
                        quadtempErr = re.findall("\[\-([0-9.]+),\+([0-9.]+)\] \(?(deg|degrees|arcmin|arcseconds|arc sec)", line)[0]
                        tempErr = [max(float(tempErr[0]), float(quadtempErr[0]), float(quadtempErr[1])), quadtempErr[2]]
                    except:
                        pass
        
        # Rearrange the error
        if float(tempErr[0])> 0 and (LocType == 1 or LocType == 2):
            self._errC = float("{}".format(tempErr[0]))
            if tempErr[1] in ['degrees', 'deg']:
                self._errC*=60
            elif tempErr[1] in ['arcseconds', 'arcsec', 'arc sec', 's', '"']:
                self._errC/=60  
        elif self._tel == 'MASTER' and self._obs:
            self._errC = 0.01
        else:
            self._errC = ''
        
        try:
            self._errC = '{:.5f}'.format(self._errC)
        except:
            pass
        
    def __Redshift__(self, verbose=False):
        
        ZType = 0

        # Scrapping the redshift        
        for line in self.gcn[5:]:
            try:
                tempZ = re.findall("(z\s?=|redshift.of.about)\s?([0-9.]+)\s?\+?\/?\-?\s?([0-9.]+)?", line)[0][1:]
                if "mag" not in line:
                    ZType = 1
                if verbose==3:
                    print("Redshift is adopted from  \n'{}'".format(line))

            except:
                try:
                    tempZ = re.findall("estimate is ([0-9.]+)\s?\+?\/?\-?\s?([0-9.]+)?\s?(Mpc|kpc)", line)[0]
                    ZType = 2
                    if verbose==3:
                        print("Redshift is adopted from \n'{}'".format(line))
                except:
                    pass

        
        # Rearrange the redshift   
        if ZType == 1:
            if tempZ[0][-1] == ".":
                z = tempZ[0][:-1]
            else:
                z = tempZ[0]
                
            if tempZ[1] == "":
                zerr = 1/10**len(str(tempZ[0]).split(".")[1])
            elif tempZ[1][-1] == ".":
                zerr = tempZ[1][:-1]
            else:
                zerr = tempZ[1]
                
            self.Redshift = ['{:.4f}'.format(float(z)), float(zerr)]
        elif ZType == 2:
            if tempZ[2] == 'Mpc':
                unit = u.Mpc
            elif tempZ[2] == 'kpc':
                unit = u.kpc
            d0 = Distance(float(tempZ[0]), unit=unit)
            z0 = d0.compute_z()
            derr1 = Distance(float(tempZ[0])+float(tempZ[1]), unit=unit)
            zerr1 = derr1.compute_z()
            derr2 = Distance(float(tempZ[0])-float(tempZ[1]), unit=unit)
            zerr2 = derr2.compute_z()
            zerr = max(abs(zerr1-z0), abs(z0-zerr2))
            self.Redshift = ['{:.4f}'.format(z0), zerr]
            
            
    def __saveFile__(self, csvfile):
        fieldnames = ['GCN_number', 'Telescope', 'Detector', 'Event_name', 'Detection', 'Trigger_time',  'RA', 'RA(hms)', 'DEC', 'DEC(hms)', 'Error', 'Redshift', 'z_Error', 'MW']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        for tel in self._tel:
            for dtr, i in zip(self._dtr[tel], range(np.size(self._dtr[tel]))):
                writer.writerow({'GCN_number': self.gcnNum, 
                         'Telescope': tel,
                         'Detector': dtr,
                         'Event_name': self.Event,
                         'Detection': self._multidtr[1][i],
                         'Trigger_time': self.Trigger,
                         'RA': self.Loc[0][0] if (self._multidtr[1][i] != False) else '',
                         'RA(hms)': self.Loc[1][0] if (self._multidtr[1][i] != False) else '',
                         'DEC': self.Loc[0][1] if (self._multidtr[1][i] != False) else '',
                         'DEC(hms)': self.Loc[1][1] if (self._multidtr[1][i] != False) else '',
                         'Error': self._errC if (self._multidtr[1][i] != False) else '',
                         'Redshift': self.Redshift[0] if (self._multidtr[1][i] != False) else '',
                         'z_Error': self.Redshift[1] if (self._multidtr[1][i] != False) else '',
                         'MW': self._list_of_tel[tel][dtr] if (self._multidtr[1][i] != False) else ''})
            
    def Save_in_csv(self, verbose=False, save=True, online=False):

        self.__scrapping__(verbose=verbose)
        
        if self._report:
            if save:
                with open('GCN_list.csv', 'a') as csvfile:
                    self.__saveFile__(csvfile)
                            
            if online:
                with open('/Users/dtak/Documents/GitHub/NewGCN/GCN_list.csv', 'a') as csvfile:
                    self.__saveFile__(csvfile)
                
            if verbose>=1:
                print("="*80)
                if np.size(self._dtr[self._tel[0]]) > 1 and np.size(self._tel) == 1:
                    print("GCN {}: {} reported by {} with {} out of {} detectors".format(self.gcnNum, self.Event, self._tel[0], int(self._obs), np.size(self._dtr[self._tel[0]])))
                elif (np.size(self._tel) == 1) and (self._tel[0] == self._dtr[self._tel[0]][0]):
                    print("GCN {}: {} reported by {}".format(self.gcnNum, self.Event, self._tel[0]))
                elif np.size(self._tel) > 1:
                    print("GCN {}: {} reported by".format(self.gcnNum, self.Event), end=' ')
                    for tel in self._tel[:-1]:
                        print("{}".format(tel), end=', ')
                    print("and {}".format(self._tel[-1]))
                else:
                    print("GCN {}: {} reported by {}/{}".format(self.gcnNum, self.Event, self._tel[0], self._dtr[self._tel[0]][0]))
                if self._obs:
                    print("<<<< The observatory detected a transient. >>>>")
                    print("Trigger time: {}".format(self.Trigger))
                    if self.Loc[0][0] != '':
                        print("Localization: {}({}), {}({})".format(self.Loc[0][0], self.Loc[1][0], self.Loc[0][1], self.Loc[1][1]))
                        try:
                            print("Localization error: {:.4f} in arcmin".format(float(self._errC)))
                        except:
                            pass
                    if self.Redshift[0] != '':
                        print("Redshift: {} +/- {:.5f}".format(self.Redshift[0], self.Redshift[1]))
                else:
                    print("<<<< The observatory did NOT detect a transient. >>>>")
                if verbose>=2:
                    raw_html = simple_get(self.GCNlink)
                    bs = BeautifulSoup(raw_html, 'html.parser')
                    print("-"*80)
                    print(bs.prettify())
                print("="*80)    
                print("Completed")
            else:
                if verbose:
                    if verbose>=2:
                        raw_html = simple_get(self.GCNlink)
                        bs = BeautifulSoup(raw_html, 'html.parser')
                        print(bs.prettify())
                    print("="*80)
                    print("Error")
                
    def addObservatory(self, Tel, Dtr):
        self._list_of_tel[Tel] = Dtr
        data = json.dumps(self._list_of_tel)
        f = open("List_of_Tel.json","w")
        f.write(data)
        f.close()
        with open("List_of_Tel.json") as f:
            self._list_of_tel = json.load(f)
            
    def newFile(self):
         with open('GCN_list.csv', 'w') as csvfile:
            fieldnames = ['GCN_number', 'Telescope', 'Detector', 'Event_name', 'Detection', 'Trigger_time',  'RA', 'RA(hms)', 'DEC', 'DEC(hms)', 'Error', 'Redshift', 'z_Error', 'MW']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()


            