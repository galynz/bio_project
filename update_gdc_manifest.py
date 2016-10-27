# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 17:45:26 2016

@author: gal
"""

import csv, os, hashlib, sys

def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), "b"):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def main():
    infile = open(sys.argv[1], 'rb')
    outfile = open(sys.argv[1] + '.final', 'wb')
    csv_infile = csv.DictReader(infile, dialect=csv.excel_tab)
    csv_outfile = csv.DictWriter(outfile, csv_infile.fieldnames, dialect=csv.excel_tab)
    for line in csv_infile:
        if os.path.exists(os.path.join([line['id'], line['filename']])):
            if md5(os.path.join([line['id'], line['filename']])) != line['md5']:
                csv_outfile.writerow(line)
        else:
            csv_outfile.writerow(line)
    infile.close()
    outfile.close()
    
    
if __name__ == "__main__":
    main()