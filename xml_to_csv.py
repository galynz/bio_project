# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 09:17:17 2016

@author: gal
"""

import xml.etree.ElementTree as ET
import glob
import csv
import sys

def xml_to_csv(xml_paths, output_path):
    keys = ['path']
    rows = []
    try:
        for path in glob.glob(xml_paths):
            iter_tree = ET.iterparse(path)
            row_dict = {'path' : path}
            for event, elem in iter_tree:
                key = elem.tag.split("}")[-1]
                #key = elem.tag
                value = elem.text
                n = 0
                while row_dict.has_key(key):
                    n+=1
                    key = key.split("#")[0] + "#%s" % n
                row_dict[key] = value
                if key not in keys:
                    keys.append(key)
            rows.append(row_dict)
            
        with open(output_path, 'wb') as f:
            #csv_writer = csv.DictWriter(f, sorted(list(keys)))
            csv_writer = csv.DictWriter(f, keys)
            csv_writer.writeheader()
            csv_writer.writerows(rows)
    except Exception, e:
        print e
        
def main():
    if len(sys.argv) != 3:
        print "./xml_to_csv.py <xml_paths> <output_path>"
        sys.exit(1)
    xml_to_csv(sys.argv[1], sys.argv[2])
    print "Saved", sys.argv[1]
    
if __name__ == "__main__":
    main()