import os,sys

"""
This script needs to be run always after updating the LJ params
"""
class Mol:
    def __init__(self,rtfFile,prmFile):
	self.rtfFile = 	rtfFile
	self.prmFile = prmFile
	self.atomTypes = []
	self.mass ={}
	self.newAtomTypes = {}
	self.bonds = []  	
	self.angles = []	
	self.diheds = []
	self.imps = []
	self.nbs =[]
		
    def getAtomTypes(self):
    	fin=open(self.rtfFile,"r")	
    	for line in fin:
	    token=line.rstrip().split()
	    if len(token)>=3:
	   	if token[0]=="MASS":
		   atom = token[2].upper()	
		   if atom not in self.atomTypes:
		   	self.atomTypes.append(atom)
			self.mass[atom]=token[3]
	fin.close()

    def reNameAtomTypes(self):
    	for a in self.atomTypes:
	    if "_" in a:
		at = a.split("_")	
	    	newAtomType = at[0].upper()+at[1].upper()+"GP"
	    else:newAtomType = a.upper()+"GP"	
	    self.newAtomTypes[a.upper()]=newAtomType
 
    def updateRTF(self):
	fin = open(self.rtfFile)
	fout = open("mol_final.rtf","w")

	print >>fout,"* Topology File."
	print >>fout,"*"
	print >>fout,"   99   1 \n"

	###Writing all atom types under MASS here
	for a in self.atomTypes:
	    atomType = self.newAtomTypes[a.upper()]
	    atMass = self.mass[a]	  		
	    idx = -1		    
	    fout.write('{:9} {:4} {:8} {:10}\n'.format("MASS",idx,atomType,atMass))		

	###Updating ATOM section here
	for line in fin:
	    token=line.rstrip().split()
	    if len(token)>=1:
		if token[0]=="ATOM":	
		    atomName = token[1]
		    if "_" in token[2]:
		    	at = token[2].split("_")
			atomType = at[0].upper()+at[1].upper()+"GP"	
	    	    else:atomType =self.newAtomTypes[token[2].upper()] 
		    charge = token[3]	
		    fout.write('{:5} {:6} {:7} {:10}\n'.format("ATOM",atomName,atomType,charge))	
	    	else:
		    if token[0]=="MASS" or token[0]=="*" or token[0]=="99":continue
		    else:print>>fout,line.rstrip()
	    else:print>>fout,line.rstrip()	

	fout.close() 
	fin.close()
  				
    def getParams(self):
	fin = open(self.prmFile,"r")
	bondFlag=0;angleFlag=0;dihedFlag=0;impFlag=0;nbFlag=0	
	for line in fin:
	    token = line.rstrip().split()
	    try:	
	    	if token[0]=="BONDS":
		    bondFlag=1
		elif token[0]=="ANGLES":
		    bondFlag=0
		    angleFlag=1
		elif token[0]=="DIHEDRALS":
		    bondFlag=0
		    angleFlag=0
		    dihedFlag=1
		elif token[0]=="IMPROPERS":
		    bondFlag=0
		    angleFlag=0
		    dihedFlag=0
		    impFlag=1
		elif token[0]=="NONBONDED":dihedFlag=0;impFlag=0
		elif token[0]=="!" and token[1]=="(KCAL/MOL)":
		   nbFlag=1				
		elif bondFlag==1:
		   self.bonds.append([token[0],token[1],token[2],token[3]])  	
		elif angleFlag==1:
		   self.angles.append([token[0],token[1],token[2],token[3],token[4]])  	
		elif dihedFlag==1:
		   self.diheds.append([token[0],token[1],token[2],token[3],token[4],token[5],token[6]])  	
		elif impFlag==1:
		   self.imps.append([token[0],token[1],token[2],token[3],token[4],token[5],token[6]])  	
		elif nbFlag==1:		
		   self.nbs.append([token[0],token[1],token[2],token[3],token[4],token[5],token[6]])		
	    except IndexError:continue
	fin.close()				


    def writeUpdatedPrm(self):
	fout = open("mol_final.prm","w")
	print>>fout,"* Force Field Parameter File."
	print>>fout,"*\n"

	print>>fout,"ATOMS"
	###Writing all atom types under MASS here
        for a in self.atomTypes:
            atomType = self.newAtomTypes[a.upper()]
            atMass = self.mass[a]
            idx = -1 ### charmm/namd determines its own index/prevents any index clashes
            fout.write('{:9} {:4} {:8} {:10}\n'.format("MASS",idx,atomType,atMass))

	print>>fout,"\nBONDS"
	for b in self.bonds:
	    try: b1 = self.newAtomTypes[b[0].upper()]
	    except KeyError: b1 =b[0]
	    try: b2 = self.newAtomTypes[b[1].upper()]
	    except KeyError: b2 =b[1]
	    k = b[2]
	    b0 = b[3]		
	    fout.write('{:6} {:6} {:6} {:6}\n'.format(b1,b2,k,b0))	
	fout.write("\n")

	print>>fout,"ANGLES"
	for a in self.angles:
	    try: a1 = self.newAtomTypes[a[0].upper()]
	    except KeyError: a1 =a[0]
	    try: a2 = self.newAtomTypes[a[1].upper()]
	    except KeyError: a2 =a[1]
	    try: a3 = self.newAtomTypes[a[2].upper()]
	    except KeyError: a3 =a[2]
	    k = a[3]
	    a0 = a[4]				
	    fout.write('{:6} {:6} {:6} {:8} {:8}\n'.format(a1,a2,a3,k,a0))
	fout.write("\n")

	print>>fout,"DIHEDRALS"
	for d in self.diheds:
	    try: d1 = self.newAtomTypes[d[0].upper()]
	    except KeyError: d1 = d[0]
	    try: d2 = self.newAtomTypes[d[1].upper()]
	    except KeyError: d2 = d[1]
	    try: d3 = self.newAtomTypes[d[2].upper()]
	    except KeyError: d3 = d[2]	
	    try: d4 = self.newAtomTypes[d[3].upper()]
	    except KeyError: d4 = d[3]	
	    k = d[4];n=d[5];pa=d[6]
	    fout.write('{:6} {:6} {:6} {:6} {:8} {:5} {:6}\n'.format(d1,d2,d3,d4,k,n,pa))
	fout.write("\n")

	print>>fout,"IMPROPERS"
	for i in self.imps:
	    try: i1 = self.newAtomTypes[i[0].upper()]
	    except KeyError: i1 = i[0]	
	    try: i2 = self.newAtomTypes[i[1].upper()]
	    except KeyError: i2 = i[1]	
	    try: i3 = self.newAtomTypes[i[2].upper()]
	    except KeyError: i3 = i[2]	
	    try: i4 = self.newAtomTypes[i[3].upper()]
	    except KeyError: i4 = i[3]	
	    k = i[4];n = i[5];pa = i[6]
	    fout.write('{:6} {:6} {:6} {:6} {:8} {:5} {:6}\n'.format(i1,i2,i3,i4,k,n,pa)) 	
	fout.write("\n")

	print>>fout,"NONBONDED  NBXMOD 5  GROUP SWITCH CDIEL -"
	print>>fout,"CUTNB 14.0  CTOFNB 12.0  CTONNB 10.0  EPS 1.0  E14FAC 0.83333333  WMIN 1.4"
	print>>fout,"!                Emin     Rmin/2              Emin/2     Rmin  (for 1-4's)"
	print>>fout,"!             (kcal/mol)    (A)"
	for n in self.nbs:	
	    try: at = self.newAtomTypes[n[0].upper()]   		
	    except KeyError: at = n[0]
	    fout.write('{:6} {:6} {:12} {:12} {:6} {:12} {:12}\n'.format(at,n[1],n[2],n[3],n[4],n[5],n[6]))	 

	print>>fout,"\nEND"
	fout.close()

    def runUpdate(self):
	self.getAtomTypes()
   	self.reNameAtomTypes()
   	self.updateRTF()
   	self.getParams()
   	self.writeUpdatedPrm()
   		    

if __name__=="__main__":	
   baseDir = sys.argv[1] #GAAMP working dir
   os.chdir(baseDir+"/results")
   
   prm = "";rtf ="" 
   if os.path.isfile("mol-tor.prm"):
	prm = "mol-tor.prm"
   elif os.path.isfile("mol.prm"):
	prm = "mol.prm"
   print "UPDATING_PRM:",prm
	
   if os.path.isfile("mol-tor.rtf"):
	rtf = "mol-tor.rtf"
   elif os.path.isfile("mol-esp-wat.rtf"):
	rtf = "mol-esp-wat.rtf"	
   elif os.path.isfile("mol-esp.rtf"):
	rtf = "mol-esp.rtf"
   elif os.path.isfile("mol.rtf"):		 	
	rtf = "mol.rtf"
   print "UPDATING_RTF:",rtf

   myRtf = Mol(rtf,prm)
   myRtf.runUpdate()	

   os.chdir(baseDir+"/../")    ### back to $GAAMPRUN	


