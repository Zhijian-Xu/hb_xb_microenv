#!/usr/bin/env bash

#select X...Y halogen bonds.

#awk '$13 == "Not_organic_symbol:O" || $13 == "Not_organic_symbol:N" || $13 == "Not_organic_symbol:S" ' halogen_PDB_XRay_NMR_2016-6-2.out >  halogen_PDB_XRay_NMR_2016-6-2_50671.txt

vdw_Cl=1.75
vdw_Br=1.85
vdw_I=1.98
vdw_O=1.52
vdw_N=1.55
vdw_S=1.80

#delta_d=0.5
#degree=120
file=$1
delta_d=$2
degree=$3
#for X...O
num=`awk -F ";" '$13 == "Not_organic_symbol:O" && $6 == "halogen_symbol:CL" && $15 < '"$(echo $vdw_Cl + $vdw_O + $delta_d | bc)"' && $15 > 1.5 && $16 > '$degree'' $file | wc -l`
awk -F ";" '$13 == "Not_organic_symbol:O" && $6 == "halogen_symbol:CL" && $15 < '"$(echo $vdw_Cl + $vdw_O + $delta_d | bc)"' && $15 > 1.5 && $16 > '$degree'' $file > halogen_PDB_XRay_NMR_O-Cl_${delta_d}_${degree}_${num}.txt
cat halogen_PDB_XRay_NMR_O-Cl_${delta_d}_${degree}_${num}.txt > halogen_PDB_XRay_NMR_${delta_d}_${degree}_tmp$$.txt
num=`awk -F ";" '$13 == "Not_organic_symbol:O" && $6 == "halogen_symbol:BR" && $15 < '"$(echo $vdw_Br + $vdw_O + $delta_d | bc)"' && $15 > 1.5 && $16 > '$degree'' $file | wc -l`
awk -F ";" '$13 == "Not_organic_symbol:O" && $6 == "halogen_symbol:BR" && $15 < '"$(echo $vdw_Br + $vdw_O + $delta_d | bc)"' && $15 > 1.5 && $16 > '$degree'' $file > halogen_PDB_XRay_NMR_O-Br_${delta_d}_${degree}_${num}.txt
cat halogen_PDB_XRay_NMR_O-Br_${delta_d}_${degree}_${num}.txt >> halogen_PDB_XRay_NMR_${delta_d}_${degree}_tmp$$.txt
num=`awk -F ";" '$13 == "Not_organic_symbol:O" && $6 == "halogen_symbol:I" && $15 < '"$(echo $vdw_I + $vdw_O + $delta_d | bc)"' && $15 > 1.5 && $16 > '$degree'' $file | wc -l`
awk -F ";" '$13 == "Not_organic_symbol:O" && $6 == "halogen_symbol:I" && $15 < '"$(echo $vdw_I + $vdw_O + $delta_d | bc)"' && $15 > 1.5 && $16 > '$degree'' $file > halogen_PDB_XRay_NMR_O-I_${delta_d}_${degree}_${num}.txt
cat halogen_PDB_XRay_NMR_O-I_${delta_d}_${degree}_${num}.txt >> halogen_PDB_XRay_NMR_${delta_d}_${degree}_tmp$$.txt


#for X...N
num=`awk -F ";" '$13 == "Not_organic_symbol:N" && $6 == "halogen_symbol:CL" && $15 < '"$(echo $vdw_Cl + $vdw_N + $delta_d | bc)"' && $15 > 1.5 && $16 > '$degree'' $file | wc -l`
awk -F ";" '$13 == "Not_organic_symbol:N" && $6 == "halogen_symbol:CL" && $15 < '"$(echo $vdw_Cl + $vdw_N + $delta_d | bc)"' && $15 > 1.5 && $16 > '$degree'' $file > halogen_PDB_XRay_NMR_N-Cl_${delta_d}_${degree}_${num}.txt
cat  halogen_PDB_XRay_NMR_N-Cl_${delta_d}_${degree}_${num}.txt >> halogen_PDB_XRay_NMR_${delta_d}_${degree}_tmp$$.txt
num=`awk -F ";" '$13 == "Not_organic_symbol:N" && $6 == "halogen_symbol:BR" && $15 < '"$(echo $vdw_Br + $vdw_N + $delta_d | bc)"' && $15 > 1.5 && $16 > '$degree'' $file | wc -l`
awk -F ";" '$13 == "Not_organic_symbol:N" && $6 == "halogen_symbol:BR" && $15 < '"$(echo $vdw_Br + $vdw_N + $delta_d | bc)"' && $15 > 1.5 && $16 > '$degree'' $file > halogen_PDB_XRay_NMR_N-Br_${delta_d}_${degree}_${num}.txt
cat halogen_PDB_XRay_NMR_N-Br_${delta_d}_${degree}_${num}.txt >> halogen_PDB_XRay_NMR_${delta_d}_${degree}_tmp$$.txt
num=`awk -F ";" '$13 == "Not_organic_symbol:N" && $6 == "halogen_symbol:I" && $15 < '"$(echo $vdw_I + $vdw_N + $delta_d | bc)"' && $15 > 1.5 && $16 > '$degree'' $file | wc -l`
awk -F ";" '$13 == "Not_organic_symbol:N" && $6 == "halogen_symbol:I" && $15 < '"$(echo $vdw_I + $vdw_N + $delta_d | bc)"' && $15 > 1.5 && $16 > '$degree'' $file > halogen_PDB_XRay_NMR_N-I_${delta_d}_${degree}_${num}.txt
cat halogen_PDB_XRay_NMR_N-I_${delta_d}_${degree}_${num}.txt >> halogen_PDB_XRay_NMR_${delta_d}_${degree}_tmp$$.txt

#for X...S
num=`awk -F ";" '$13 == "Not_organic_symbol:S" && $6 == "halogen_symbol:CL" && $15 < '"$(echo $vdw_Cl + $vdw_S + $delta_d | bc)"' && $15 > 1.5 && $16 > '$degree'' $file | wc -l`
awk -F ";" '$13 == "Not_organic_symbol:S" && $6 == "halogen_symbol:CL" && $15 < '"$(echo $vdw_Cl + $vdw_S + $delta_d | bc)"' && $15 > 1.5 && $16 > '$degree'' $file > halogen_PDB_XRay_NMR_S-Cl_${delta_d}_${degree}_${num}.txt
cat halogen_PDB_XRay_NMR_S-Cl_${delta_d}_${degree}_${num}.txt >> halogen_PDB_XRay_NMR_${delta_d}_${degree}_tmp$$.txt
num=`awk -F ";" '$13 == "Not_organic_symbol:S" && $6 == "halogen_symbol:BR" && $15 < '"$(echo $vdw_Br + $vdw_S + $delta_d | bc)"' && $15 > 1.5 && $16 > '$degree'' $file | wc -l`
awk -F ";" '$13 == "Not_organic_symbol:S" && $6 == "halogen_symbol:BR" && $15 < '"$(echo $vdw_Br + $vdw_S + $delta_d | bc)"' && $15 > 1.5 && $16 > '$degree'' $file > halogen_PDB_XRay_NMR_S-Br_${delta_d}_${degree}_${num}.txt
cat halogen_PDB_XRay_NMR_S-Br_${delta_d}_${degree}_${num}.txt >> halogen_PDB_XRay_NMR_${delta_d}_${degree}_tmp$$.txt
num=`awk -F ";" '$13 == "Not_organic_symbol:S" && $6 == "halogen_symbol:I" && $15 < '"$(echo $vdw_I + $vdw_S + $delta_d | bc)"' && $15 > 1.5 && $16 > '$degree'' $file | wc -l`
awk -F ";" '$13 == "Not_organic_symbol:S" && $6 == "halogen_symbol:I" && $15 < '"$(echo $vdw_I + $vdw_S + $delta_d | bc)"' && $15 > 1.5 && $16 > '$degree'' $file > halogen_PDB_XRay_NMR_S-I_${delta_d}_${degree}_${num}.txt
cat halogen_PDB_XRay_NMR_S-I_${delta_d}_${degree}_${num}.txt >> halogen_PDB_XRay_NMR_${delta_d}_${degree}_tmp$$.txt
XB_total_num=`awk -F ";" 'END{print NR}' halogen_PDB_XRay_NMR_${delta_d}_${degree}_tmp$$.txt`
echo "The total number of XBs are: $XB_total_num"

#XBs from X-Ray experiment
# awk -F ";" 'NR==FNR {hash[$1]=$0} NR != FNR && ($1 in hash) {print $0, hash[$1]}' pdb_X-Ray_105315_20160415.txt  XB_MainChain_${delta_d}_${degree}_pdb.txt > XB_MainChain_${delta_d}_${degree}_X-Ray_res.txt
# XB_XRay_MainChain_num=`awk 'END{print NR}' XB_MainChain_${delta_d}_${degree}_X-Ray_res.txt`
# echo "The number of XBs for MainChain from X-Ray are: $XB_XRay_MainChain_num"
# awk -F ";" 'NR==FNR {hash[$1]=$0} NR != FNR && ($1 in hash) {print $0, hash[$1]}' pdb_X-Ray_105315_20160415.txt  XB_SideChain_${delta_d}_${degree}_pdb.txt > XB_SideChain_${delta_d}_${degree}_X-Ray_res.txt
# XB_XRay_SideChain_num=`awk 'END{print NR}' XB_SideChain_${delta_d}_${degree}_X-Ray_res.txt`
# echo "The number of XBs for SideChain from X-Ray are: $XB_XRay_SideChain_num"



