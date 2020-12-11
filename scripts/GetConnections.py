import os,sys

#AMBERHOME = "/home/huanglei/tools/amber11/"

class Mol:
      def __init__(self,bonds,angles,diheds,improps,ljs=None):
	self.bonds = bonds
	self.angles = angles
	self.diheds = diheds
	self.improps = improps
	self.ljs = ljs


def runAC(q):
    os.system("$AMBERHOME/AmberTools/bin/antechamber -i mol-opt.pdb -fi pdb -o mol-ac -fo charmm -j both -c bcc -at gaff -nc "+str(q))


def getAtomCounts():
    fin = open("mol-ac.rtf","r")
    atomCounts = {}	
    for line in fin:
	token = line.rstrip().split()
	try:
	    if token[0]=="ATOM":
		atm = filter(str.isalpha, token[1]).upper()
		if not atomCounts.has_key(atm): atomCounts[atm]=1
		else: atomCounts[atm] = atomCounts[atm]+1
	except IndexError: continue
    return atomCounts	


def getConnections(atomCounts):
    atoms=[];bonds=[];angles=[];diheds=[];improps=[]
    atomF=1; bondF=1; angleF=1; dihedF=1; impropF=1
    fin = open("mol-ac.rtf","r")
    i=0
    for line in fin:
	token = line.rstrip().split()
	try:
	   if token[0]=="ATOM":
		atomF=0
		i+=1
	   elif token[0]=="BOND": 
		atomF=1
		bondF=0
	   elif token[0]=="ANGL":
		bondF=1
		angleF=0
	   elif token[0]=="DIHE":
		angleF=1
		dihedF=0
	   elif token[0]=="IMPH":
		dihedF=1
		impropF=0
	   if atomF==0:
		atm = filter(str.isalpha, token[1]).upper() 
		try: 
		     n = int(filter(str.isdigit, token[1]))+1
		     atm = atm+str(n)	
		except ValueError: 
		     if atomCounts[atm]>1:
			atm = atm+str(1)	 
		atoms.append(atm)	
	   elif bondF==0: 
		atm1 = filter(str.isalpha, token[1]).upper()
		atm2 = filter(str.isalpha, token[2]).upper()
		try: 
		     n1 = int(filter(str.isdigit, token[1]))+1 
		     atm1 = atm1+str(n1)		
		except ValueError: 
		     if atomCounts[atm1]>1:
			atm1 = atm1+str(1)
		try: 
		     n2 = int(filter(str.isdigit, token[2]))+1
		     atm2 = atm2+str(n2)	
		except ValueError: 
		     if atomCounts[atm2]>1:
			atm2 = atm2+str(1)
		bonds.append(atm1+"-"+atm2)
	   elif angleF==0: 
		atm1 = filter(str.isalpha, token[1]).upper()
		atm2 = filter(str.isalpha, token[2]).upper()
		atm3 = filter(str.isalpha, token[3]).upper()
		try: 
		     n1 = int(filter(str.isdigit, token[1]))+1 
		     atm1 = atm1+str(n1)	
		except ValueError: 
		     if atomCounts[atm1]>1:	
			atm1 = atm1+str(1)
		try: 
		     n2 = int(filter(str.isdigit, token[2]))+1
		     atm2 = atm2+str(n2)	
		except ValueError: 
		     if atomCounts[atm2]>1:
			atm2 = atm2+str(1) 
		try: 
		     n3 = int(filter(str.isdigit, token[3]))+1
		     atm3 = atm3+str(n3)	
		except ValueError: 
		     if atomCounts[atm3]>1:	
			atm3 = atm3+str(1)
		angles.append(atm1+"-"+atm2+"-"+atm3)
	   elif dihedF==0: 
                atm1 = filter(str.isalpha, token[1]).upper()
                atm2 = filter(str.isalpha, token[2]).upper()
                atm3 = filter(str.isalpha, token[3]).upper()
                atm4 = filter(str.isalpha, token[4]).upper()
                try: 
		     n1 = int(filter(str.isdigit, token[1]))+1
		     atm1 = atm1+str(n1) 	
                except ValueError: 
		     if atomCounts[atm1]>1:     	
			atm1 = atm1+str(1)
                try: 
		     n2 = int(filter(str.isdigit, token[2]))+1
		     atm2 = atm2+str(n2)	
                except ValueError: 
		     if atomCounts[atm2]>1:
			atm2 = atm2+str(1)	
                try: 
		     n3 = int(filter(str.isdigit, token[3]))+1
		     atm3 = atm3+str(n3)	
                except ValueError: 
		     if atomCounts[atm3]>1:	
			atm3 = atm3+str(1)
                try: 
		     n4 = int(filter(str.isdigit, token[4]))+1
		     atm4 = atm4+str(n4)	
                except ValueError: 
		     if atomCounts[atm4]>1:	
			atm4 = atm4+str(1)
		diheds.append(atm1+"-"+atm2+"-"+atm3+"-"+atm4)
	   elif impropF==0:
                atm1 = filter(str.isalpha, token[1]).upper()
                atm2 = filter(str.isalpha, token[2]).upper()
                atm3 = filter(str.isalpha, token[3]).upper()
                atm4 = filter(str.isalpha, token[4]).upper()
                try: 
		     n1 = int(filter(str.isdigit, token[1]))+1
		     atm1 = atm1+str(n1)	   	
                except ValueError: 
		     if atomCounts[atm1]>1:	
			atm1 = atm1+str(1)
                try: 
		     n2 = int(filter(str.isdigit, token[2]))+1
		     atm2 = atm2+str(n2)	
                except ValueError: 
		     if atomCounts[atm2]>1:	
			atm2 = atm2+str(1)
                try: 
		     n3 = int(filter(str.isdigit, token[3]))+1
		     atm3 = atm3+str(n3)	
                except ValueError: 
		     if atomCounts[atm3]>1:	
			atm3 = atm3+str(1)
                try: 
		     n4 = int(filter(str.isdigit, token[4]))+1
		     atm4 = atm4+str(n4)	
                except ValueError: 
		     if atomCounts[atm4]>1:	
			atm4 = atm4+str(1)
                improps.append(atm1+"-"+atm2+"-"+atm3+"-"+atm4)
	except IndexError: continue
    fin.close()

    return atoms,bonds,angles,diheds,improps

def getTypes(cgenRtf):
    fin = open(cgenRtf,"r")
    flag=0
    aTypes = {}
    for line in fin:
	token = line.rstrip().split()
	try:
	   if token[0]=="ATOM" and flag==0:
	   	aTypes[token[1]] = token[2]
    	   if token[0]=="BOND": flag=1
	except IndexError: continue
    return aTypes

def getConnectionTypes(cons,aTypes):
    conTypes=[]
    for ats in cons:
	at = ats.split("-")
    	cTypeF = '-'.join([aTypes[a] for a in at])
	at.reverse()
        cTypeR = '-'.join([aTypes[a] for a in at])
	if cTypeF in conTypes or cTypeR in conTypes: continue
	else: conTypes.append(cTypeF)
    return conTypes

def getPrm(cGenPrm,mol):	
    atomF=1; bondF=1; angleF=1; dihedF=1; impropF=1; ljF=1	
    bonds={};angles={};diheds={};improps={};ljs={}
    fin = open(cGenPrm,"r")
    for line in fin:
	token = line.rstrip().split()
	try:
	    if token[0]=="BONDS":bondF=0
	    elif token[0]=="ANGLES":
		bondF=1
		angleF=0
	    elif token[0]=="DIHEDRALS":
		angleF=1
		dihedF=0
	    elif token[0]=="IMPROPERS":
		dihedF=1
		impropF=0
	    elif token[0]=="NONBONDED":
		impropF=1
		ljF=0
	    if bondF==0 and token[0]!="BONDS":
		bonds[token[0]+"-"+token[1]] = [token[2],token[3]]
	    elif angleF==0 and token[0]!="ANGLES":
		if len(token)>5 and token[5][0]!="!":
			angles[token[0]+"-"+token[1]+"-"+token[2]] = [token[3],token[4],token[5],token[6]]		
		else: angles[token[0]+"-"+token[1]+"-"+token[2]] = [token[3],token[4]]		
	    elif dihedF==0 and token[0]!="DIHEDRALS": 	
		if not diheds.has_key(token[0]+"-"+token[1]+"-"+token[2]+"-"+token[3]):
		   diheds[token[0]+"-"+token[1]+"-"+token[2]+"-"+token[3]] = [[token[4],token[5],token[6]]]
		else: diheds[token[0]+"-"+token[1]+"-"+token[2]+"-"+token[3]].append([token[4],token[5],token[6]])
	    elif impropF==0 and token[0]!="IMPROPERS":
		improps[token[0]+"-"+token[1]+"-"+token[2]+"-"+token[3]] = [token[4],token[5],token[6]]
	    elif ljF==0 and token[0]!="NONBONDED" and token[0]!="END":
		ljs[token[0]] = [token[2],token[3]]						
	except IndexError: continue
    fin.close()
    myPMol = Mol(bonds,angles,diheds,improps,ljs)
    return myPMol

def compareParas(pExpected,pAvail):
    pAvailR =  [p.split("-") for p in pAvail]
    for p in pAvailR: p.reverse()
    pAvailR = ["-".join(p) for p in pAvailR]
    for p in pAvailR: pAvail.append(p)
    misses = set(pExpected)-set(pAvail)
    return list(misses)

def getMissingTypes(pTypes,pAvail,isDihed=False):
    misses = compareParas(pTypes,pAvail.keys())
    for m in misses:
	if len(m.split("-"))==2:
	   pAvail[m] = ["0.0","0.0"]
	elif len(m.split("-"))==3:
	   pAvail[m] = ["0.0","0.0"]
	elif len(m.split("-"))==4:
	   if isDihed: pAvail[m] = [["0.0","0.0","0.0"]]
	   else: pAvail[m] = ["0.0","0.0","0.0"] 	
    return pAvail 

def writePRM(mol):
    fout = open("mol-tmp.prm","w")
    print >> fout, "* MINI FORCE FIELD PARAMETER FILE."
    print >> fout, "*\n"
    print >> fout, "BONDS"

    for p in mol.bonds:
        bs = p.split("-")
        print >> fout, bs[0],bs[1],mol.bonds[p][0],mol.bonds[p][1]
    print >> fout, "\n"

    print >> fout, "ANGLES"
    for p in mol.angles:
        ang = p.split("-")
        if len(mol.angles[p])>2:
                print >> fout, ang[0],ang[1],ang[2],mol.angles[p][0],mol.angles[p][1],mol.angles[p][2],mol.angles[p][3]
        else:
                print >> fout, ang[0],ang[1],ang[2],mol.angles[p][0],mol.angles[p][1]
    print >> fout, "\n"

    print >> fout, "DIHEDRALS"
    for p in mol.diheds:
        ds = p.split("-")
	if len(mol.diheds[p]) > 1:
           for v in mol.diheds[p]:
                print >> fout, ds[0],ds[1],ds[2],ds[3],v[0],v[1],v[2]
        else: print >> fout, ds[0],ds[1],ds[2],ds[3],mol.diheds[p][0][0],mol.diheds[p][0][1],mol.diheds[p][0][2]
    print >> fout, "\n"

    print >> fout, "IMPROPERS"
    for p in mol.improps:
        ims = p.split("-")
        print >> fout, ims[0],ims[1],ims[2],ims[3],mol.improps[p][0],mol.improps[p][1],mol.improps[p][2]
    print >> fout, "\n"

    print >> fout, "NONBONDED  E14FAC  1.000000"
    for p in mol.ljs:
        print >> fout, p, "0.00", mol.ljs[p][0], mol.ljs[p][1]

    print >> fout, "\n"

    print >> fout, "END"
    fout.close()



if __name__=="__main__":
   q = sys.argv[1]
   runAC(q)
   aTypes = getTypes("mol.rtf")
   aCounts = getAtomCounts()
   atoms,bonds,angles,diheds,improps = getConnections(aCounts)		

   ptBonds = getConnectionTypes(bonds,aTypes) 
   ptAngles = getConnectionTypes(angles,aTypes) 
   ptDiheds = getConnectionTypes(diheds,aTypes) 
   ptImprops = getConnectionTypes(improps,aTypes)

   myMol = Mol(ptBonds,ptAngles,ptDiheds,ptImprops)
   os.system("cp mol.prm mol-cgen.prm")
   pAvailMol = getPrm("mol.prm",myMol)

   allPBonds = getMissingTypes(myMol.bonds,pAvailMol.bonds)
   allPAngles = getMissingTypes(myMol.angles,pAvailMol.angles)
   allPDiheds = getMissingTypes(myMol.diheds,pAvailMol.diheds,isDihed=True)
   
   pAvailMol.bonds = allPBonds 
   pAvailMol.angles = allPAngles
   pAvailMol.diheds = allPDiheds
  
   writePRM(pAvailMol)
   

