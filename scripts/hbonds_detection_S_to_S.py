#coding=utf-8
#python3.8
#识别complex中的氢键，并输出到相应的文件中
#for ligand element S to protein S

from concurrent.futures import process
from pymol import cmd
import os
import pandas as pd

def list_hb(selection,selection2=None,cutoff=3.2,angle=55,mode=1,hb_list_name='hbonds'):
  """
  USAGE

  list_hb selection, [selection2 (default=None)], [cutoff (default=3.2)], [angle (default=55)], [hb_list_name]

  The script automatically adds a requirement that atoms in the selection (and selection2 if used)
  must be either of the elements N or O.

  If mode is set to 0, then no angle cutoff is used, otherwise the angle cutoff is used
  and defaults to 55 degrees.

  e.g.
  To get a list of all H-bonds within chain A of an object
    list_hb 1abc & c. a &! r. hoh, cutoff=3.2, hb_list_name=abc-hbonds

  To get a list of H-bonds between chain B and everything else:
    list_hb 1tl9 & c. b, 1tl9 &! c. b

  """
  cutoff=float(cutoff)
  angle=float(angle)
  mode=float(mode)
  #model = cmd.get_names("objects")
# ensure only N and O atoms are in the selection
  #selection = selection + " & e. n+o"
  if not selection2:
    hb = cmd.find_pairs(selection,selection,mode=mode,cutoff=cutoff,angle=angle)
  else:
    #selection2 = selection2 + " & e. n+o"
    hb = cmd.find_pairs(selection,selection2,mode=mode,cutoff=cutoff,angle=angle)
#find_paris这个函数对于加H和不加的结构，输出的结果都是一样的
# sort the list for easier reading
#hb.sort(lambda x,y:(cmp(x[0][1],y[0][1])))
  resi_atom_info = []
  lig_atom_info = []
  distance_all = []
  for pairs in hb:
    #print cmd.get_names("objects")

    res_space = {'residue_atom_info': []}
    cmd.iterate("%s and index %s" % (pairs[0][0],pairs[0][1]),"residue_atom_info.append([model,resn,resi,name,elem,chain,index,alt])",space=res_space) #提取蛋白残基原子的信息
    lig_space = {'lig_atom_info': []}
    cmd.iterate("%s and index %s" % (pairs[1][0],pairs[1][1]),"lig_atom_info.append([model,resn,resi,name,elem,chain,index,alt])",space=lig_space)#提取配体原子信息
    distance = round(float(cmd.distance(hb_list_name,"%s and index %s" % (pairs[0][0],pairs[0][1]),"%s and index %s" % (pairs[1][0],pairs[1][1]))),4)


    resi_atom_info.append(list(map(lambda x:str(x),res_space["residue_atom_info"][0])))
    lig_atom_info.append(list(map(lambda x:str(x),lig_space["lig_atom_info"][0])))
    distance_all.append(str(distance))

  return resi_atom_info, lig_atom_info, distance_all
  

  
  


pdb_file = os.popen("ls /home/databank/zlp/PDB_halogen_bond_20220913/pdb_with_halogen_element/processed_pdb/")
pdb_file_list = pdb_file.read().split("\n")
#hbond_df = pd.DataFrame(index=pdb_file_list,columns=["pdb_file", "residue_number", "residue atom index","residue atom name", "ligand_number", "ligand atom index","ligand atom name", "distance"])
#hbond_df = pd.DataFrame(pdb_file_list)
      
found_no_hbonds = []
with open("hbonds_statistic_PDB_S_to_S.csv","w") as f:
  f.write("pdb_id,pdb_file,residue name,residue number,residue atom name,residue atom elem,residue chain,residue atom index,residue atom alt,ligand name,ligand number,ligand atom name,ligand atom elem,ligand chain,ligand atom index,ligand atom alt,distance" + "\n")
  for pdb in pdb_file_list:
    if pdb: #如果不是空值的话
      path = "/home/databank/zlp/PDB_halogen_bond_20220913/pdb_with_halogen_element/processed_pdb/" + pdb
      selection = "polymer.protein and e. s" #蛋白重原子
      selection2 = "organic and e. s" #残基重原子
      cmd.delete("all")
      cmd.load(path)
      # cmd.select("sele1","organic and e. n+o")
      lig_heavy_atoms = atom = cmd.count_atoms("organic")
      if lig_heavy_atoms == 0:
        print(pdb + " has no ligands")
      resi_atom, lig_atom, distances = list_hb(selection, selection2, cutoff=3.8,angle=55,mode=1,hb_list_name='hbonds') #这里的角度就相当于125度
      #print(pdb, resi_atom, lig_atom)
      pairs = len(resi_atom)
      if pairs != len(lig_atom) or pairs != len(distances): #蛋白上的原子数和配体原子数应该是一一对应的，如果个数对不上，说明有问题啊
        print("error")
      
      if pairs == 0:
        found_no_hbonds.append(pdb)
      
      for pair in range(pairs):
        f.write(resi_atom[pair][0].split("_")[0]+",")
        f.write(",".join(resi_atom[pair]) + ",")
        f.write(",".join(lig_atom[pair][1:]) + ",")
        f.write(distances[pair] + "\n")
        
print("These systems have no hbonds.")
print(len(found_no_hbonds), found_no_hbonds) 
# os.system("rm -rf found_no_hbonds_2/*") 
# for pdb in found_no_hbonds:
    # path = "processed_selected_pdb_2/" + pdb
    # os.environ["path"] = path
    # os.system("cp $path found_no_hbonds_2/ ")








