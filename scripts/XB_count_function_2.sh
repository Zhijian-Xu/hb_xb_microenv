#the interaction between protein mainchain and X
delta_d=$1
degree=$2

XB_total_num=`awk 'END{print NR}' halogen_PDB_XRay_NMR_${delta_d}_${degree}_tmp*.txt`
echo "The total number of XBs are: $XB_total_num"

for name in `cat ../20_AA_triple.txt`
do
#echo ${name}
side_chain_num=0
grep -i "Not_organic_resn:${name}" halogen_PDB_XRay_NMR_${delta_d}_${degree}_tmp*.txt | awk -F ";" '$12 == "Not_organic_name:O" || $12 == "Not_organic_name:N" || $12 == "Not_organic_name:OXT"' >> tmp$$
cat tmp$$ | sort -u > XB_protein_MainChain_${delta_d}_${degree}.txt
done
XB_MainChain_num=`awk 'END{print NR}' XB_protein_MainChain_${delta_d}_${degree}.txt`
echo "The number of XBs for MainChain are: $XB_MainChain_num"
rm tmp$$

#the interaction between protein sidechain and X
for name in `cat ../20_AA_triple.txt`
do
#echo ${name}
side_chain_num=0
grep -i "Not_organic_resn:${name}" halogen_PDB_XRay_NMR_${delta_d}_${degree}_tmp*.txt | awk -F ";" '$12 != "Not_organic_name:O" && $12 != "Not_organic_name:N" && $12 != "Not_organic_name:OXT"' >> tmp$$
cat tmp$$ | sort -u > XB_protein_SideChain_${delta_d}_${degree}.txt
done
XB_SideChain_num=`awk 'END{print NR}' XB_protein_SideChain_${delta_d}_${degree}.txt`
echo "The number of XBs for SideChain are: $XB_SideChain_num"
rm tmp$$

cat XB_protein_MainChain_${delta_d}_${degree}.txt XB_protein_SideChain_${delta_d}_${degree}.txt > XB_protein_${delta_d}_${degree}.txt
#求两个文件的差集,非标准残基也归到了not protein里面
sort halogen_PDB_XRay_NMR_${delta_d}_${degree}_tmp*.txt XB_protein_${delta_d}_${degree}.txt XB_protein_${delta_d}_${degree}.txt |uniq -u > XB_not_protein_${delta_d}_${degree}.txt
XB_not_protein_num=`awk 'END{print NR}' XB_not_protein_${delta_d}_${degree}.txt`
echo "The number of XBs for not protein are: $XB_not_protein_num"

#add pdb ID in the first column
awk -F ";" '{print substr($1, 1, 4)";", $0}' XB_protein_MainChain_${delta_d}_${degree}.txt > XB_protein_MainChain_${delta_d}_${degree}_pdb.txt
awk -F ";" '{print substr($1, 1, 4)";", $0}' XB_protein_SideChain_${delta_d}_${degree}.txt > XB_protein_SideChain_${delta_d}_${degree}_pdb.txt
awk -F ";" '{print substr($1, 1, 4)";", $0}' XB_protein_${delta_d}_${degree}.txt > XB_protein_${delta_d}_${degree}_pdb.txt
awk -F ";" '{print substr($1, 1, 4)";", $0}' XB_not_protein_${delta_d}_${degree}.txt > XB_not_protein_${delta_d}_${degree}_pdb.txt

