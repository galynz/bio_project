#! /bin/bash
for cancer in $(cat cancer_list.txt); do
    echo $cancer;
    python2.7 mutations_summary.py -c $cancer -p "/groups/nshomron/galyankovitz/storage_1_root/simple_nucleotide_variation/new_data/$cancer/*/*.maf.gz" --clinical_path "/groups/nshomron/galyankovitz/storage_1_root/simple_nucleotide_variation/new_data/$cancer/*/*.xml" -o "outputs_20160818/"$cancer"_20160818_without_silent_mut" -m "Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Missense_Mutation,Nonsense_Mutation,Splice_Site,Translation_Start_Site,Nonstop_Mutation,3'UTR,3'Flank,5'UTR,5'Flank,IGR,Intron,RNA,Targeted_Region"
done
