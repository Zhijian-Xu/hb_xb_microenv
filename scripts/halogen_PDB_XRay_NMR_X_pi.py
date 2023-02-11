#create:  zjxu@simm.ac.cn 4/10/2012
#modify: 4/18/2012 for NMR structure
#function: find the halogen in the org, determinte the residues of Not_organic within 8 angstrom of the halogen. Print the distance of pi...X and the angle of pi...X-bond(X). The "neighbor" algebra in pymol will select the atom which bond to X.  To determine a O...X-C angle, I add "and e. c" to# sele1_bond, if you want X bond to an element not registrited to C, you should delete "and e. c"  
#Note: sele1, sele2, and sele_3 should be only one atom.
#modify polymer to Not_organic to include the X...water interaction
#selection =  "((" + model + " and chain " + halogen_chain + " and i. " + halogen_resi + " and name " + halogen_name + ") around 8.0)" + " and " + model      Other organics are also included except the one contarins X which is very usefull to include X...cofactor interaction.
#modify: 4/19/2012 for X...pi interaction
#modify: change the input file to X_4133.lst in 4/3/2013 for un updated PDB survey.
#Update: change the input file to organic_ClBrI_2893.lst in 4/6/2014 for un updated PDB survey. 
#Update: set ignore_case, on in 7/29/2016 (The ignore_case setting (default: on, except in PyMOL 1.8.0.0 - 1.8.0.4) controls whether PyMOL does case sensitive matching of atomic identifiers and selection operators in the selection language.).
#Update: add sys.argv on 10/22/2019

import sys
if (len(sys.argv) != 2):
  print('Usage: pymol -rkqc halogen_PDB_XRay_NMR_X_pi.py -- file_name')
  print('For example: pymol -rkqc halogen_PDB_XRay_NMR_X_pi.py -- organic_Cl_Br_I_5651.lst')
  sys.exit()

from pymol import cmd
cmd.set('ignore_case', 1)
from pymol import stored
import math
#sys.path.append("/home/zjxu/sh/pymol")
import com
DEG_TO_RAD=3.1415926/180.0
RAD_TO_DEG=180.0/3.1415926

file = sys.argv[1]
PATH="/home/databank/zlp/PDB_halogen_bond_20220913/pdb_with_halogen_element/processed_pdb/"

def calcDistance(lx,ly,lz,px,py,pz):
  return sqrt((lx-px)**2+(ly-py)**2+(lz-pz)**2);

def calcVectorLength(*vector):
  len = math.sqrt((vector[0])**2+(vector[1])**2+(vector[2])**2);
  return len

def multiply2vectors(*vector_two):
  vector1 = vector_two[0:3]
  vector2 = vector_two[3:6]
  return (vector1[0]*vector2[0]+vector1[1]*vector2[1]+vector1[2]*vector2[2])

def calcVectorsAngle(*vector_two):
  vector1 = vector_two[0:3]
  vector2 = vector_two[3:6]
  len1=calcVectorLength(*vector1)
  len2=calcVectorLength(*vector2)
  if(len1 < 0.005):
     len1 = 0.005
  if(len2 < 0.005):
     len2 = 0.005
  dotMetrix = multiply2vectors(*vector_two)
  cos_angle = dotMetrix/(len1*len2)
  if(cos_angle > 1.0):
    #maybe there is a rounding error, sometimes cos_angle is a little bitter than 1.0
    print("cos_angle > 1.0, this cos_angle is:" + str(cos_angle))
    cos_angle = 1.0
  angle = math.acos(cos_angle)
  return angle

def calcNormalVector(*vector_two):
  vector1 = vector_two[0:3]
  vector2 = vector_two[3:6]
  nv=[]
  nvx=vector1[1]*vector2[2]-vector1[2]*vector2[1]
  nvy=vector1[2]*vector2[0]-vector1[0]*vector2[2]
  nvz=vector1[0]*vector2[1]-vector1[1]*vector2[0]
  nv=[nvx,nvy,nvz]
  return nv

def normalize(*vector):
  vlen=calcVectorLength(*vector)
  if(vlen<0.005):
    vlen=0.005
  nvector=(vector[0]/vlen,vector[1]/vlen,vector[2]/vlen)
  return nvector

#get the coordinates of the selection
def get_coord_zj(selection_tmp):
  stored.xyz = []
  cmd.iterate_state(1,selection_tmp,"stored.xyz.append([x,y,z])")
  sele_coord = [stored.xyz[0][0], stored.xyz[0][1],stored.xyz[0][2]]
  return sele_coord


#fp = open('organic_Cl_Br_I_3710.lst')
fp = open(file)
outfile = file.split("/")[-1].split(".")[0]+"_halogen_PDB_XRay_NMR_X_pi.lst"
f = open(outfile,"w")
num=0
for eachline in fp:
  #print eachline,
  eachline = eachline.strip()
  #eachline.split(" ")
  #print eachline
  tmp_list = []
  #aa.append(eachline[0])
  tmp_list = eachline.split(",")
  model = tmp_list[0]
  #print "model", model
  #halogen_resi = tmp_list[1]
  
  #using pymol to print out the environment around 8.0 Angstrom of halogen
  #cmd.load(PATH + model + ".pdb")
  cmd.load(PATH + model + ".pdb")
  #selection = model + " and (i. " + halogen_resi + " around 8.0)"
  #I select all the halogen in the pdb

  #halogen_all_select = model + " and (e. F or e. Cl or e. Br or e. I) and organic"
  halogen_all_select = model + " and (e. Cl or e. Br or e. I) and organic"
  #halogen_all_select = model + " and ((e. Cl or e. CL or e. cl) or (e. Br or e. BR or e. br) or (e. I or e. i)) and organic"

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
  
      #the Not_organic within 8.0 angstrom of halogen
        #selection =  "((" + model + " and chain " + halogen_chain + " and i. " + halogen_resi + " and name " + halogen_name + ") around vectorcom8.0)" + " and Not_organic"
        #selection =  "((" + model + " and chain " + halogen_chain + " and i. " + halogen_resi + " and name " + halogen_name + ") around 8.0)" + " and Not_organic and " + model
        #selection =  "((" + model + " and chain " + halogen_chain + " and i. " + halogen_resi + " and name " + halogen_name + ") around 8.0)" + " and " + model + " and not " + "((" + model + " and chain " + halogen_chain + " and i. " + halogen_resi + "))"
        #for Phe, Tyr, His, and Try, add 'polymer' to the selection.
        # model within 8.0 of X and not organic (which X resides)
        # I add byres to inculde all the atoms in the ring.
        selection =  "(byres ((" + model + " and chain " + halogen_chain + " and i. " + halogen_resi + " and name " + halogen_name + ") around 8.0)" + " and polymer and " + model + " and not " + "((" + model + " and chain " + halogen_chain + " and i. " + halogen_resi + ")) )" + " and r. phe+tyr+his+trp and name CG"
        #print"selection: ", selection
        #cmd.select("halo_8.0A_Not_organic", selection)
        #CG exist in phe+tyr+his+tyr, so it was used to find the ring
        cmd.select("halo_8.0A_Not_organic", selection)
  
        halo_Not_organic_atom_num = cmd.count_atoms("halo_8.0A_Not_organic")
        #print "halo_Not_organic_atom_num: ", halo_Not_organic_atom_num
        if halo_Not_organic_atom_num > 0:
          atoms = cmd.get_model("halo_8.0A_Not_organic")
          for at in atoms.atom:
            #print "CG_atom: ", at.name
            if at.alt == "":
              sele_ring = model + " and ( (byres(halo_8.0A_Not_organic and chain %s and resi %s and name %s)) and not name c+n+o+ca+cb+oh and not e. H )"% (at.chain, at.resi, at.name)
            else:
              sele_ring = model + " and ( (byres(halo_8.0A_Not_organic and chain %s and resi %s and name %s and alt %s)) and not name c+n+o+ca+cb+oh and not e. H)"% (at.chain, at.resi, at.name, at.alt)
            #print "sele_ring: ", sele_ring
            #print "sele_ring: ", cmd.count_atoms(sele_ring)

            #cmd.pseudoatom("sele2", sele_ring)
            #pseodoatom could not used in NRM mode structure in the command mode. I donot know why. So com was used to determine the center of the ring.
            #sele2 is the center of ring
            com.com(sele_ring, object="sele2")
            #print "sele2: ", cmd.count_atoms("sele2")

            #print "sele1: ", sele1
            #for the plane vector
            try: #exclude some index error
              sele2_xyz = get_coord_zj("sele2")
              CG_xyz = get_coord_zj(sele_ring + " and name CG")
              CD2_xyz = get_coord_zj(sele_ring + " and name CD2")
              #print CG_xyz, CD2_xyz, sele2_xyz
              vec1 = [sele2_xyz[0]-CG_xyz[0], sele2_xyz[1]-CG_xyz[1], sele2_xyz[2]-CG_xyz[2]]
              vec2 = [sele2_xyz[0]-CD2_xyz[0], sele2_xyz[1]-CD2_xyz[1], sele2_xyz[2]-CD2_xyz[2]]
              vec1_2 = vec1 + vec2
              nv=calcNormalVector(*vec1_2)
            except Exception as e:
              print(model + " has problems", e.__class__.__name__,e)
              continue

            sele1_xyz=get_coord_zj(sele1)
            sele1_sele2 = [sele1_xyz[0]-sele2_xyz[0], sele1_xyz[1]-sele2_xyz[1], sele1_xyz[2]-sele2_xyz[2]] 
            alpha_vector = nv + sele1_sele2
            alpha = calcVectorsAngle(*alpha_vector)
            alpha = alpha*RAD_TO_DEG
            if(alpha > 90.0):
              alpha = 180.0 - alpha
            dst = cmd.distance("tmp", sele1, "sele2")
            theta_ang = cmd.angle("tmp_angle", "sele2", sele1, sele_3)
          #print at.chain, at.resn, at.resi, at.name, at.index, "%8.3f"%dst
          #print at.hetatm
          #print atoms.atom[0].__dict__
            #print "sele1: " + sele1
            #print "sele2: " + sele2
            f.write(model +";"+ "halogen_chain:" + halogen_chain +";"+ "halogen_resi:" + halogen_resi+";"+"halogen_resn:" + halogen_resn+";" + "halogen_name:" + halogen_name+";" + "halogen_symbol:" + halo.symbol+";" + "halogen_alt:" + halo.alt+";" + "Not_organic_hetatm:%i"%at.hetatm+";" + "Not_organic_chain:" + at.chain+";" + "Not_organic_resi:" + at.resi+";" + "Not_organic_resn:" + at.resn+";" + "Not_organic_name:" + at.name+";" + "Not_organic_symbol:" + at.symbol+";" + "Not_organic_alt:" + at.alt+";" + "dist:%.2f"%dst+";" + "theta:%.2f"%theta_ang+";" + "alpha:%.2f"%alpha + ";" + "halogen_index:" + str(halogen_index) + "\n")
            cmd.delete("sele2")
        #cmd.delete("halo_8.0A_Not_organic")
  num += 1
  cmd.delete("all")
print(num)
f.close()
fp.close()
