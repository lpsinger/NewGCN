import numpy as np
import matplotlib as plt
import re
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import Distance

from IPython.display import clear_output
import csv
import json

from xml.etree.ElementTree import Element, SubElement, Comment
from xml.etree import ElementTree
from xml.dom import minidom

def prettify(elem):
    """Return a pretty-printed XML string for the Element."""
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")



class MakeTable:

    def __init__(self, verbose=False, save=True, online=False):
        
        self.GCN_list = []
        
        with open('GCN_list.csv', newline='') as csvfile:
            fieldnames = ['GCN_number', 'Telescope', 'Detector', 'Event_name', 'Detection', 'Trigger_time',  'RA', 'RA(hms)', 'DEC', 'DEC(hms)', 'Error', 'Redshift', 'z_Error', 'MW']
            csv_list = csv.DictReader(csvfile, fieldnames=fieldnames)
            
            for row in csv_list:
                self.GCN_list.append(row)


        self.GCN_table = Element('table')

        self._evts = []
        for gcn in self.GCN_list:
            if gcn['Event_name'] not in self._evts:
                self._evts.append(gcn['Event_name'])

        self._evts = [[i[-7:], i] for i in self._evts]
        
        self._evts.sort(reverse=True)
        
        for evt in self._evts:
            self.__addEvts__(evt[1])
        
        if verbose:
            print(prettify(self.GCN_table))
            
        if save:
            with open("xml_table.xml", "w") as f:
                f.write(prettify(self.GCN_table))
            with open("./static/testapp/xml_table.xml", "w") as f:
                f.write(prettify(self.GCN_table))
                
        if online:
            with open("/Users/dtak/Documents/GitHub/NewGCN/xml_table.xml", "w") as f:
                f.write(prettify(self.GCN_table))    
            


    def __addEvts__(self, evt):
        comment = Comment('{}'.format(evt))
        self.GCN_table.append(comment)

        grouped_gcn = []
        for gcn in self.GCN_list:
            if gcn['Event_name'] == evt:
                grouped_gcn.append(list(gcn.values()))

        grouped_gcn = np.asarray(grouped_gcn)

        # Column 1: Event name
        temp = evt.split(' ')    
        if temp[0] == "GRB":
            evtlink = "/other/{}".format(temp[1])
            src = "1"
        elif temp[0] == "IC":
            evtlink = "/other/icecube_{}".format(temp[1])
            src = "3"
        else:
            evtlink = "/other/{}{}".format(temp[0],temp[1])
            src = "2"
            
        year = evt[-7:-5]

        row_init = SubElement(self.GCN_table, 'Transient')
        row_yr = SubElement(row_init, 'Year{}'.format(year))
        row_evt = SubElement(row_yr, 'Src{}'.format(src))

        col_evt = SubElement(row_evt, 'Event')
        col_evt.set('link', evtlink)
        col_evt.text = '{}'.format(evt)

        # Column 2: Earliest Trigger Time
        try:
            t0from = grouped_gcn[grouped_gcn[:,5] == np.sort([time for time in grouped_gcn[:,5] if len(time)>10])[0]][0]
        except:
            t0from = grouped_gcn[0]
            
        gcnlink = "/gcn3/{}".format(t0from[0])

        col_t0 = SubElement(row_evt, 'T0')
        col_t0.set('link', gcnlink)
        col_t0.set('reporter', self.__reportName__(t0from))
        col_t0.text = '{}'.format(t0from[5])

        # Column 3, 4, and 5: Localization
        try:
            locfrom = grouped_gcn[grouped_gcn[:,10] ==np.sort([err for err in grouped_gcn[:,10] if err != ''])[0]][0]
            gcnlink = "/gcn3/{}".format(locfrom[0])

            col_ra = SubElement(row_evt, 'RA')
            col_ra.set('link', gcnlink)
            col_ra.set('units', 'deg')
            col_ra.set('reporter', self.__reportName__(locfrom))
            col_ra.set('hms', locfrom[7])
            col_ra.text = '{}'.format(locfrom[6])

            col_dec = SubElement(row_evt, 'DEC')
            col_dec.set('link', gcnlink)
            col_dec.set('units', 'deg')
            col_dec.set('reporter', self.__reportName__(locfrom))
            col_dec.set('dms', locfrom[9])
            col_dec.text = '{}'.format(locfrom[8])

            err, err_u = self.__errorCircle__(locfrom)
            col_err = SubElement(row_evt, 'Error')
            col_err.set('link', gcnlink)
            col_err.set('units', err_u)
            col_err.set('reporter', self.__reportName__(locfrom))
            col_err.text = '{}'.format(err)

        except:
            col_ra = SubElement(row_evt, 'RA')
            col_ra.text = ' '

            col_dec = SubElement(row_evt, 'DEC')
            col_dec.text = ' '

            col_err = SubElement(row_evt, 'Error')
            col_err.text = ' '

        # Column 6: Redshift
        try:    
            zfrom = grouped_gcn[grouped_gcn[:,12] == np.sort([err for err in grouped_gcn[:,12] if err != ''])[0]][0]
            gcnlink = "/gcn3/{}".format(zfrom[0])

            col_z = SubElement(row_evt, 'Redshift')
            col_z.set('link', gcnlink)
            col_z.set('reporter', self.__reportName__(zfrom))
            col_z.set('error', zfrom[12])
            col_z.text = '{}'.format(zfrom[11])
        except:
            col_z = SubElement(row_evt, 'Redshift')
            col_z.text = ' '

        # Column 7: List of Observatories
        col_tel = SubElement(row_evt, 'Inst')

        reported = []
        for row in grouped_gcn:
            try:
                reported = np.vstack([reported, [self.__reportName__(row), row[4], row[0]]])
            except:
                reported = np.asarray([[self.__reportName__(row), row[4], row[0]]])

        reported_dtrs = list(set(reported[:,0]))
        reported_dtrs.sort()

        for dtr in reported_dtrs:
            if max(reported[reported[:,0] == dtr][:,1]) == "True":
                each_tel = SubElement(col_tel, 'Tel')
                each_tel.set('link', "/gcn3/{}".format(max(reported[reported[:,0] == dtr][:,2])))
                each_tel.set('obs', 'on')
                each_tel.text = dtr
                
        for dtr in reported_dtrs:
            if max(reported[reported[:,0] == dtr][:,1]) == "False":
                each_tel = SubElement(col_tel, 'Tel')
                each_tel.set('link', "/gcn3/{}".format(max(reported[reported[:,0] == dtr][:,2])))
                each_tel.set('obs', 'off')
                each_tel.text = dtr

        # Column 8: Multiwavelength
        col_mw = SubElement(row_evt, 'MW')
        mmw_obs = list(set(grouped_gcn[:,13]))
        mw_list = ['radio', 'optical', 'X-ray', '&gamma;-ray', 'HE', 'VHE']
        EM_flag = False
        for mw, i in zip(mw_list, range(6)):
            each_mw = SubElement(col_mw, 'Band')
            if mw in mmw_obs:
                each_mw.set('obs', '{}'.format(i))
                EM_flag = True
            else:
                each_mw.set('obs', 'off')
            each_mw.text = mw
        
        # Column 9: Multimessenger
        col_mm = SubElement(row_evt, 'MM')
        mm_list = ['EM', 'GW', '&nu;']
        for mm, i in zip(mm_list, range(3)):
            each_mm = SubElement(col_mm, 'Mgr')
            if mm == 'EM' and EM_flag:
                each_mm.set('obs', 'on')
            elif mm in mmw_obs:
                each_mm.set('obs', 'on')
            else:
                each_mm.set('obs', 'off')
            each_mm.text = mm
    
    def __reportName__(self, row):
        if np.size(row[2]) == 1 and row[1] != row[2]:
            reportName = "{}/{}".format(row[1], row[2])
        else:
            reportName = "{}".format(row[1])
        return reportName
    
    def __errorCircle__(self, row):
        if row[10] == '':
            err = ''
        elif float(row[10]) > 6:
            err = "{:.2f} &deg;".format(float(row[10])/60.)
            err_u = 'deg'
        elif float(row[10]) < 0.6:
            err = '{:.2f} "'.format(float(row[10])*60.)
            err_u = 'arcmin'
        else:
            err = "{:.2f} '".format(float(row[10]))
            err_u = 'arcsec'
        return err, err_u    
        
