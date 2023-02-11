#python 3.8
#coding=utf-8
#create:  zjxu@mail.shcnc.ac.cn 4/10/2012
#modify: 4/18/2012 for NMR structure
#function: find the halogen in the org, determinte the atoms of Not_organic within 10.0 angstrom of the halogen. Print the distance of O...X and the angle of O...X-bond(X). The "neighbor" algebra in pymol will select the atom which bond to X.  To determine a O...X-C angle, I add "and e. c" to# sele1_bond, if you want X bond to an element not registrited to C, you should delete "and e. c"  
#Note: sele1, sele2, and sele_3 should be only one atom.
#modify polymer to Not_organic to include the X...water interaction
#selection =  "((" + model + " and chain " + halogen_chain + " and i. " + halogen_resi + " and name " + halogen_name + ") around 10.0)" + " and " + model      Other organics are also included except the one contarins X which is very usefull to include X...cofactor interaction.
#Update: change the input file from X_3566.lst to X_4133.lst in 4/3/2013 for un updated PDB survey.
#Update: change the input file from X_4133.lst to organic_ClBrI_2893.lst in 4/6/2014 for un updated PDB survey.
#Although the script considered FClBrI, the input PDB did not consider F on 4/6/2014.
#Update: change "e. c" to "e. c or e. C" because PyMOL1.8.1.0 could not identify element C using e. c. in 6/1/2016
#    change halogen_all_select. The old one is: #halogen_all_select = model + " and (e. F or e. Cl or e. Br or e. I) and organic"
#    #in PyMOL1.8.1.0 element Br could not be identified by e. Br if the element type is BR in PDB file.
#    halogen_all_select = model + " and ((e. Cl or e. CL or e. cl) or (e. Br or e. BR or e. br) or (e. I or e. i)) and organic"
#Update: add cmd.set('ignore_case', 1) on 3/15/2017.
#Update: add sys.argv on 11/7/2019
#update: for python 3.8 14/9/2022
#update: for F 21/10/2022

# import sys
# if (len(sys.argv) != 2):
  # print('Usage: pymol -rkqc halogen_PDB_XRay_NMR.py -- file_name')
  # print('For example: pymol -rkqc halogen_PDB_XRay_NMR.py -- organic_Cl_Br_I_5651.lst')
  # sys.exit()



from pymol import cmd
from glob import glob
PATH="/home/databank/zlp/PDB_halogen_bond_20220913/pdb_with_halogen_element/processed_pdb/"

cmd.set('ignore_case', 1)

f=open("dis_angle_PDB_F_Cl_Br_I_to_protein.csv","w")
f.write("pdb_id,pdb_file,residue name,residue number,residue atom name,residue atom elem,residue chain,residue atom index,residue atom alt,ligand name,ligand number,ligand atom name,ligand atom elem,ligand chain,ligand atom index,ligand atom alt,distance,angle" + "\n")
all_pdb_files = sorted(glob(PATH+"*.pdb"))
for pdb_file in all_pdb_files:
  cmd.delete("all")
  cmd.load(pdb_file)
  model = pdb_file.split("/")[-1].split(".")[0]
  #selection = model + " and (i. " + halogen_resi + " around 10.0)"
  #I select all the halogen in the pdb
  

  #halogen_all_select = model + " and (e. F or e. Cl or e. Br or e. I) and organic"
  #in PyMOL1.8.1.0 element Br could not be identified by e. Br if the element type is BR in PDB file.
  halogen_all_select = model + " and ((e. Cl or e. CL or e. cl) or (e. Br or e. BR or e. br) or (e. I or e. i) or (e. F or e. f)) and organic"

  halogen_all = cmd.select("halogen_all", halogen_all_select)
  halo_num = cmd.count_atoms("halogen_all")
  if halo_num > 0:
    halo_atm = cmd.get_model("halogen_all")
    for halo in halo_atm.atom:
      #print "the halogen"
      halogen_chain = halo.chain
      halogen_resi = halo.resi
      halogen_resn = halo.resn
      halogen_name = halo.name
      halogen_index = halo.index
      #print "halo.alt: " + halo.alt
      if halo.alt == "":
        sele1 = "halogen_all" + " and chain " + halogen_chain + " and i. " +  halogen_resi + " and name " + halogen_name 
      else:
        sele1 = "halogen_all" + " and chain " + halogen_chain + " and i. " +  halogen_resi + " and name " + halogen_name + " and alt " + halo.alt 
      # print "sele1: ", cmd.count_atoms(sele1)
      #sele1_bond = cmd.select("bond_halogen", "(neighbor " + sele1 + " ) and e. c")
      #sele1_bond = cmd.select("bond_halogen", "(neighbor " + sele1 + " ) and e. c and " + model)
      sele1_bond = cmd.select("bond_halogen", "(neighbor " + sele1 + " ) and (e. c or e. C) and " + model)
      X_bond_atom = cmd.get_model("bond_halogen")
      # print "bond_halogen: ", cmd.count_atoms("bond_halogen")
      for X_nbr_atom in X_bond_atom.atom:
        #print sele1_bond, X_nbr_atom.resi, X_nbr_atom.name
        if X_nbr_atom.alt == "":
          sele_3 = "bond_halogen and chain %s and resi %s and name %s"% (X_nbr_atom.chain, X_nbr_atom.resi, X_nbr_atom.name)
        else:
          sele_3 = "bond_halogen and chain %s and resi %s and name %s and alt %s"% (X_nbr_atom.chain, X_nbr_atom.resi, X_nbr_atom.name, X_nbr_atom.alt)
        # print "sele_3: ", cmd.count_atoms(sele_3)
  
      #the Not_organic within 10.0 angstrom of halogen
        selection =  "((" + model + " and chain " + halogen_chain + " and i. " + halogen_resi + " and name " + halogen_name + ") around 10.0)" + " and " + model + " and not " + "((" + model + " and chain " + halogen_chain + " and i. " + halogen_resi + "))"
        # print"selection: ", selection
        cmd.select("halo_10.0A_Not_organic", selection)
  
        halo_Not_organic_atom_num = cmd.count_atoms("halo_10.0A_Not_organic")
        # print "halo_Not_organic_atom_num: ", halo_Not_organic_atom_num
        if halo_Not_organic_atom_num > 0:
          atoms = cmd.get_model("halo_10.0A_Not_organic")
          for at in atoms.atom:
            if at.alt == "":
              sele2="halo_10.0A_Not_organic and chain %s and resi %s and name %s"% (at.chain, at.resi, at.name)
            else:
              sele2="halo_10.0A_Not_organic and chain %s and resi %s and name %s and alt %s"% (at.chain, at.resi, at.name, at.alt)
            #print sele1, sele2
            # print "sele2: ", cmd.count_atoms(sele2)
            dst = cmd.distance("tmp", sele1, sele2)
            ang = cmd.angle("tmp_angle", sele2, sele1, sele_3)
          #print at.chain, at.resn, at.resi, at.name, at.index, "%8.3f"%dst
          #print at.hetatm
          #print atoms.atom[0].__dict__
            #print "sele1: " + sele1
            #print "sele2: " + sele2
            pdb_id = model.split("_")[0]
            #pdb_id,pdb_file,residue name,residue number,residue atom name,residue atom elem,residue chain,residue atom index,residue atom alt,ligand name,ligand number,ligand atom name,ligand atom elem,ligand chain,ligand atom index,ligand atom alt,distance,angle
            f.write(pdb_id+","+model+ "," +at.resn+ ","+str(at.resi)+","+at.name+","+at.symbol+","+at.chain+","+str(at.index)+","+at.alt+","+halogen_resn+","+str(halogen_resi)+","+halogen_name+","+halo.symbol+","+halo.chain+","+str(halo.index)+","+halo.alt+","+"%.2f"%dst+","+"%.2f"%ang+"\n")
        #cmd.delete("halo_10.0A_Not_organic")

f.close()
