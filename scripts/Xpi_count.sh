out=`basename $1 .lst`
cp $1 ${out}_process.lst
sed -i "s/dist://" ${out}_process.lst
sed -i "s/theta://" ${out}_process.lst
sed -i "s/alpha://" ${out}_process.lst 
#halogen_PDB_XRay_NMR_X_pi.out:
#11gs_0001 halogen_chain:A halogen_resi:211 halogen_resn:EAA halogen_name:CL1 halogen_symbol:CL halogen_alt: Not_organic_hetatm:0 Not_organic_chain:A Not_organic_resi:7 Not_organic_resn:TYR Not_organic_name:CG Not_organic_symbol:C Not_organic_alt: dist:6.40 theta:111.72 alpha:27.84

#awk -F "_" '{print $1, $0}' halogen_PDB_XRay_NMR_X_pi.out > tmp1.txt
#11gs 11gs_0001 halogen_chain:A halogen_resi:211 halogen_resn:EAA halogen_name:CL1 halogen_symbol:CL halogen_alt: Not_organic_hetatm:0 Not_organic_chain:A Not_organic_resi:7 Not_organic_resn:TYR Not_organic_name:CG Not_organic_symbol:C Not_organic_alt: dist:6.40 theta:111.72 alpha:27.84

#awk -F "dist:" '{print $1, $2}' tmp1.txt > tmp2.txt
#awk -F "theta:" '{print $1, $2}' tmp2.txt > tmp3.txt
#awk -F "alpha:" '{print $1, $2}' tmp3.txt > halogen_PDB_XRay_NMR_X_pi.txt


num_Cl=`awk -F ";" '$6 == "halogen_symbol:CL" && $15 < 4.2 && $16 > 120 && $17 < 60' ${out}_process.lst | wc -l`
awk -F ";" '$6 == "halogen_symbol:CL" && $15 < 4.2 && $16 > 120 && $17 < 60' ${out}_process.lst > Cl-pi_${num_Cl}.txt
num_Br=`awk -F ";" '$6 == "halogen_symbol:BR" && $15 < 4.3 && $16 > 120 && $17 < 60' ${out}_process.lst | wc -l`
awk -F ";" '$6 == "halogen_symbol:BR" && $15 < 4.3 && $16 > 120 && $17 < 60' ${out}_process.lst > Br-pi_${num_Br}.txt
num_I=`awk -F ";" '$6 == "halogen_symbol:I" && $15 < 4.5 && $16 > 120 && $17 < 60' ${out}_process.lst | wc -l`
awk -F ";" '$6 == "halogen_symbol:I" && $15 < 4.5 && $16 > 120 && $17 < 60' ${out}_process.lst > I-pi_${num_I}.txt

cat Cl-pi_${num_Cl}.txt Br-pi_${num_Br}.txt I-pi_${num_I}.txt | sort -u > Cl_Br_I-pi.txt
#add pdb id in the first column
awk -F ";" '{print substr($1, 1, 4)";", $0}' Cl_Br_I-pi.txt > Cl_Br_I-pi_pdb.txt
