#! /bin/bash
today=$(date +"%Y%m%d")
mkdir "outputs_"$today
for cancer in $(cat cancer_list.txt); do
    echo $cancer;
    python2.7 mutations_summary.py -c $cancer -p "/groups/nshomron/galyankovitz/storage_1_root/simple_nucleotide_variation/new_data_by_disease/$cancer/*/*.maf.gz" --clinical_path "/groups/nshomron/galyankovitz/storage_1_root/simple_nucleotide_variation/new_data_by_disease/$cancer/*/*.xml" -o "outputs_"$today"/"$cancer"_"$today"_without_silent_mut" -m "Frame_Shift_Del,Frame_Shift_Ins,In_Frame_Del,In_Frame_Ins,Missense_Mutation,Nonsense_Mutation,Splice_Site,Translation_Start_Site,Nonstop_Mutation,3'UTR,3'Flank,5'UTR,5'Flank,IGR,Intron,RNA,Targeted_Region,De_novo_Start_InFrame,De_novo_Start_OutOfFrame"
done
