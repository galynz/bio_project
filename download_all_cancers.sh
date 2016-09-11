#! /bin/bash
for cancer in $(cat cancer_list.txt); do
    echo $cancer;
    cd $cancer
    gdc-client download -m *.tsv
    cd ..
done
