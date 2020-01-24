from check_link import simple_get
from requests import get
from bs4 import BeautifulSoup, Tag
import csv
import numpy as np

def edit_csv(gcnNum, saved = True, verbose=False):

    fieldnames = ['GCN_number', 'Telescope', 'Detector', 'Event_name', 'Detection', 'Trigger_time',  'RA', 'RA(hms)', 'DEC', 'DEC(hms)', 'Error', 'Redshift', 'z_Error', 'MW']
    
    GCNlink = "https://gcn.gsfc.nasa.gov/gcn3/{}.gcn3".format(gcnNum)
    raw_html = simple_get(GCNlink)
    bs = BeautifulSoup(raw_html, 'html.parser')
    print(bs.prettify())

    if verbose: print("="*80)

    with open('GCN_list.csv', newline='') as csvfile:
        csv_list = csv.DictReader(csvfile, fieldnames=fieldnames)
        csv_new_list = []
        for row in csv_list:
            if float(row['GCN_number']) == gcnNum:
                for field in fieldnames:
                    entered = input("{}: {}, ok? ".format(field, row[field]))
                    if entered != '':
                        row[field] = entered
                        if verbose:
                            print("{} --> {}".format(row[field], entered))

                if saved:
                    np.save(gcnNum, row)
            
            csv_new_list.append(row)
            

    with open('GCN_list.csv', 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        for row in csv_new_list:
            writer.writerow(row)

    if verbose: 
        print("GCN {} is successfully updated.".format(gcnNum))
        print("="*80)

def edit_csv_auto(gcnNum, verbose=False):
    
    try:
        savedrow = np.load('{}.npy'.format(gcnNum))
        savedrow = savedrow.item()
    except:
        print("There is no modified GCN{}".format(gcnNum))
        print("Please do manually")
        return

    fieldnames = ['GCN_number', 'Telescope', 'Detector', 'Event_name', 'Detection', 'Trigger_time', 'RA', 'RA(hms)', 'DEC', 'DEC(hms)', 'Error', 'Redshift', 'z_Error', 'MW']    
    if verbose: print("="*80)

    with open('GCN_list.csv', newline='') as csvfile:
        csv_list = csv.DictReader(csvfile, fieldnames=fieldnames)
        csv_new_list = []
        for row in csv_list:
            if float(row['GCN_number']) == gcnNum:
                for field in fieldnames:
                    if (row[field] != savedrow[field]) and verbose:
                        print("{} --> {}".format(row[field], savedrow[field]))
                csv_new_list.append(savedrow)
            else:
                csv_new_list.append(row)

    with open('GCN_list.csv', 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        for row in csv_new_list:
            writer.writerow(row)

    if verbose: 
        print("GCN {} is successfully updated.".format(gcnNum))
        print("="*80)