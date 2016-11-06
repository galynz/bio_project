#! /bin/bash
for cancer in $(cat cancer_list.txt); do
    echo $cancer;
    cd $cancer
    cp ~/bio_project/update_gdc_manifest.py .
    python2.7 update_gdc_manifest.py gdc_manifest.2016-10-*.tsv
    gdc-client download -m *2016-10-*.tsv.final -t ~/bio_project/token.txt
    cd ..
done
