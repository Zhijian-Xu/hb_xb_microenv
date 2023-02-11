#!/usr/bin/env bash

#select X hbonds
#pdb_id,pdb_file,residue name,residue number,residue atom name,residue atom elem,residue chain,residue atom index,residue atom alt,ligand name,ligand number,ligand atom name,ligand atom elem,ligand chain,ligand atom index,ligand atom alt,distance,angle
#1ctr,1ctr_TFP_20A,GLU,7,CA,C,A,48,,TFP,153,F1,F,A,978,,9.73,101.36
#要有氢键的话，N、O、S上的要能加H，所以主链上的O、PRO上的N\O，MET上的S，ASN上的所有O，GLN上的所有O，ASP、GLU、THR上所有的O（侧链羧基电离）都需要筛除

vdw_F=1.47
vdw_Cl=1.75
vdw_Br=1.85
vdw_I=1.98
vdw_O=1.52
vdw_N=1.55
vdw_S=1.80
exclude_O_res="PRO,ASN,GLN,ASP,GLU"

#delta_d=0.5
#degree=120
file=$1
delta_d=$2
degree=$3
#for X...O
num=`awk -F "," '$5 !="O" && "'$exclude_O_res'"!~$3 && $6 == "O" && $13 == "F" && $17 < '"$(echo $vdw_F + $vdw_O + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file | wc -l`
awk -F "," '$5 !="O" && "'$exclude_O_res'"!~$3 && $6 == "O" && $13 == "F" && $17 < '"$(echo $vdw_F + $vdw_O + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file > halogen_hbonds_drugbank_XRay_NMR_O-F_${delta_d}_${degree}_${num}.txt
cat halogen_hbonds_drugbank_XRay_NMR_O-F_${delta_d}_${degree}_${num}.txt > halogen_hbonds_drugbank_XRay_NMR_${delta_d}_${degree}_tmp$$.txt
num=`awk -F "," '$5 !="O" && "'$exclude_O_res'"!~$3 && $6 == "O" && $13 == "CL" && $17 < '"$(echo $vdw_Cl + $vdw_O + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file | wc -l`
awk -F "," '$5 !="O" && "'$exclude_O_res'"!~$3 && $6 == "O" && $13 == "CL" && $17 < '"$(echo $vdw_Cl + $vdw_O + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file > halogen_hbonds_drugbank_XRay_NMR_O-Cl_${delta_d}_${degree}_${num}.txt
cat halogen_hbonds_drugbank_XRay_NMR_O-Cl_${delta_d}_${degree}_${num}.txt >> halogen_hbonds_drugbank_XRay_NMR_${delta_d}_${degree}_tmp$$.txt
num=`awk -F "," '$5 !="O" && "'$exclude_O_res'"!~$3 && $6 == "O" && $13 == "BR" && $17 < '"$(echo $vdw_Br + $vdw_O + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file | wc -l`
awk -F "," '$5 !="O" && "'$exclude_O_res'"!~$3 && $6 == "O" && $13 == "BR" && $17 < '"$(echo $vdw_Br + $vdw_O + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file > halogen_hbonds_drugbank_XRay_NMR_O-Br_${delta_d}_${degree}_${num}.txt
cat halogen_hbonds_drugbank_XRay_NMR_O-Br_${delta_d}_${degree}_${num}.txt >> halogen_hbonds_drugbank_XRay_NMR_${delta_d}_${degree}_tmp$$.txt
num=`awk -F "," '$5 !="O" && "'$exclude_O_res'"!~$3 && $6 == "O" && $13 == "I" && $17 < '"$(echo $vdw_I + $vdw_O + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file | wc -l`
awk -F "," '$6 == "O" && $13 == "I" && $17 < '"$(echo $vdw_I + $vdw_O + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file > halogen_hbonds_drugbank_XRay_NMR_O-I_${delta_d}_${degree}_${num}.txt
cat halogen_hbonds_drugbank_XRay_NMR_O-I_${delta_d}_${degree}_${num}.txt >> halogen_hbonds_drugbank_XRay_NMR_${delta_d}_${degree}_tmp$$.txt


#for X...N
num=`awk -F "," '$3 != "PRO" && $6 == "N" && $13 == "F" && $17 < '"$(echo $vdw_F + $vdw_N + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file | wc -l`
awk -F "," '$3 != "PRO" && $6 == "N" && $13 == "F" && $17 < '"$(echo $vdw_F + $vdw_N + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file > halogen_hbonds_drugbank_XRay_NMR_N-F_${delta_d}_${degree}_${num}.txt
cat  halogen_hbonds_drugbank_XRay_NMR_N-F_${delta_d}_${degree}_${num}.txt >> halogen_hbonds_drugbank_XRay_NMR_${delta_d}_${degree}_tmp$$.txt
num=`awk -F "," '$3 != "PRO" && $6 == "N" && $13 == "CL" && $17 < '"$(echo $vdw_Cl + $vdw_N + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file | wc -l`
awk -F "," '$3 != "PRO" && $6 == "N" && $13 == "CL" && $17 < '"$(echo $vdw_Cl + $vdw_N + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file > halogen_hbonds_drugbank_XRay_NMR_N-Cl_${delta_d}_${degree}_${num}.txt
cat  halogen_hbonds_drugbank_XRay_NMR_N-Cl_${delta_d}_${degree}_${num}.txt >> halogen_hbonds_drugbank_XRay_NMR_${delta_d}_${degree}_tmp$$.txt
num=`awk -F ";" '$3 != "PRO" && $6 == "N" && $13 == "BR" && $17 < '"$(echo $vdw_Br + $vdw_N + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file | wc -l`
awk -F "," '$3 != "PRO" && $6 == "N" && $13 == "BR" && $17 < '"$(echo $vdw_Br + $vdw_N + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file > halogen_hbonds_drugbank_XRay_NMR_N-Br_${delta_d}_${degree}_${num}.txt
cat halogen_hbonds_drugbank_XRay_NMR_N-Br_${delta_d}_${degree}_${num}.txt >> halogen_hbonds_drugbank_XRay_NMR_${delta_d}_${degree}_tmp$$.txt
num=`awk -F "," '$3 != "PRO" && $6 == "N" && $13 == "I" && $17 < '"$(echo $vdw_I + $vdw_N + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file | wc -l`
awk -F "," '$3 != "PRO" && $6 == "N" && $13 == "I" && $17 < '"$(echo $vdw_I + $vdw_N + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file > halogen_hbonds_drugbank_XRay_NMR_N-I_${delta_d}_${degree}_${num}.txt
cat halogen_hbonds_drugbank_XRay_NMR_N-I_${delta_d}_${degree}_${num}.txt >> halogen_hbonds_drugbank_XRay_NMR_${delta_d}_${degree}_tmp$$.txt

#for X...S
num=`awk -F "," '$3 != "MET" && $6 == "S" && $13 == "F" && $17 < '"$(echo $vdw_F + $vdw_S + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file | wc -l`
awk -F "," '$3 != "MET" && $6 == "S" && $13 == "F" && $17 < '"$(echo $vdw_F + $vdw_S + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file > halogen_hbonds_drugbank_XRay_NMR_S-F_${delta_d}_${degree}_${num}.txt
cat halogen_hbonds_drugbank_XRay_NMR_S-F_${delta_d}_${degree}_${num}.txt >> halogen_hbonds_drugbank_XRay_NMR_${delta_d}_${degree}_tmp$$.txt
num=`awk -F "," '$3 != "MET" && $6 == "S" && $13 == "CL" && $17 < '"$(echo $vdw_Cl + $vdw_S + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file | wc -l`
awk -F "," '$3 != "MET" && $6 == "S" && $13 == "CL" && $17 < '"$(echo $vdw_Cl + $vdw_S + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file > halogen_hbonds_drugbank_XRay_NMR_S-Cl_${delta_d}_${degree}_${num}.txt
cat halogen_hbonds_drugbank_XRay_NMR_S-Cl_${delta_d}_${degree}_${num}.txt >> halogen_hbonds_drugbank_XRay_NMR_${delta_d}_${degree}_tmp$$.txt
num=`awk -F ";" '$3 != "MET" && $6 == "S" && $13 == "BR" && $17 < '"$(echo $vdw_Br + $vdw_S + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file | wc -l`
awk -F "," '$3 != "MET" && $6 == "S" && $13 == "BR" && $17 < '"$(echo $vdw_Br + $vdw_S + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file > halogen_hbonds_drugbank_XRay_NMR_S-Br_${delta_d}_${degree}_${num}.txt
cat halogen_hbonds_drugbank_XRay_NMR_S-Br_${delta_d}_${degree}_${num}.txt >> halogen_hbonds_drugbank_XRay_NMR_${delta_d}_${degree}_tmp$$.txt
num=`awk -F "," '$3 != "MET" && $6 == "S" && $13 == "I" && $17 < '"$(echo $vdw_I + $vdw_S + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file | wc -l`
awk -F "," '$3 != "MET" && $6 == "S" && $13 == "I" && $17 < '"$(echo $vdw_I + $vdw_S + $delta_d | bc)"' && $17 > 1.5 && $18 < '$degree'' $file > halogen_hbonds_drugbank_XRay_NMR_S-I_${delta_d}_${degree}_${num}.txt
cat halogen_hbonds_drugbank_XRay_NMR_S-I_${delta_d}_${degree}_${num}.txt >> halogen_hbonds_drugbank_XRay_NMR_${delta_d}_${degree}_tmp$$.txt
total_num=`awk -F "," 'END{print NR}' halogen_hbonds_drugbank_XRay_NMR_${delta_d}_${degree}_tmp$$.txt`
echo "The total number of hbonds for halogen elements are: $total_num"

#XBs from X-Ray experiment
# awk -F ";" 'NR==FNR {hash[$1]=$0} NR != FNR && ($1 in hash) {print $0, hash[$1]}' pdb_X-Ray_105315_20160415.txt  XB_MainChain_${delta_d}_${degree}_pdb.txt > XB_MainChain_${delta_d}_${degree}_X-Ray_res.txt
# XB_XRay_MainChain_num=`awk 'END{print NR}' XB_MainChain_${delta_d}_${degree}_X-Ray_res.txt`
# echo "The number of XBs for MainChain from X-Ray are: $XB_XRay_MainChain_num"
# awk -F ";" 'NR==FNR {hash[$1]=$0} NR != FNR && ($1 in hash) {print $0, hash[$1]}' pdb_X-Ray_105315_20160415.txt  XB_SideChain_${delta_d}_${degree}_pdb.txt > XB_SideChain_${delta_d}_${degree}_X-Ray_res.txt
# XB_XRay_SideChain_num=`awk 'END{print NR}' XB_SideChain_${delta_d}_${degree}_X-Ray_res.txt`
# echo "The number of XBs for SideChain from X-Ray are: $XB_XRay_SideChain_num"



