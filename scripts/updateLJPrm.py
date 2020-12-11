import os,sys
import itertools

BABELPATH=""
GAAMPDIR=""

def setBabelPath():
    fin = open(GAAMPDIR+"/opt/PATHS","r")
    for line in fin:
        token = line.rstrip().split("=")
        try:
            if token[0]=="BABELDIR":
                global BABELPATH
                BABELPATH = token[1].split('"')[1].split('"')[0]
        except IndexError: continue


def getConnections(ac_rtf):
    nameToIdx = {}; bonds={};angles=[];diheds=[];improps=[]
    atmIdx=0; group=0
    fin=open(ac_rtf,"r")
    for line in fin:
	token=line.rstrip().split()
	if len(token)>0:
	   if token[0]=="GROUP":group+=1
           if token[0]=="ATOM" and group==1:
              nameToIdx[token[1]] = atmIdx
              atmIdx+=1 
	   if token[0]=="BOND" and group==1:
	      if not bonds.has_key(caps(token[1])):bonds[caps(token[1])]=[caps(token[2])]
	      else:bonds[caps(token[1])].append(caps(token[2]))
	   elif token[0]=="ANGL" and group==1:
	      angles.append([caps(token[1]),caps(token[2]),caps(token[3])])	
	   elif token[0]=="DIHE" and group==1:
	      diheds.append([caps(token[1]),caps(token[2]),caps(token[3]),caps(token[4])])	
	   elif token[0]=="IMPH" and group==1:	
	      improps.append([caps(token[1]),caps(token[2]),caps(token[3]),caps(token[4])])	
    fin.close()
    return nameToIdx, bonds,angles,diheds,improps		

def detectAtomsInCycle(mol2):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol2", "mol2")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol,mol2)
    cycle = []
    for obatom in openbabel.OBMolAtomIter(mol):
        if obatom.IsInRing():
            cycle.append(obatom.GetIndex())
    return cycle


def isInCycle(atomIdx,cycle):
    if atomIdx in cycle: return True
    else: return False


def getAllAtomTypes(rtf):
    atoms=[]
    atomTypes=[]
    group=0
    fin=open(rtf,"r")	
    for line in fin:
	token=line.rstrip().split()
    	if len(token)>0:
	   if token[0]=="GROUP":group+=1
	   if token[0]=="ATOM" and group==1:
	       atoms.append(token[1])
	       atomTypes.append(token[2])
    return atoms,atomTypes

def createFinalAtomTypes(allAtomTypes,atom,newAtType):
    try:finalAtType = newAtType+allAtomTypes[atom].split("_")[1]
    except IndexError: finalAtType = newAtType
    return finalAtType

def expandAtomTypes(nameToIdx,cycle,bonds,atoms,allAtomTypes):
    #print atoms
    #print allAtomTypes
    allAtomTypes = {caps(atoms[i]):caps(allAtomTypes[i]) for i in range(len(atoms))}
    #print allAtomTypes
    #print bonds

    specialAtomTypes={}
    for atom in bonds:
        hCount=0
        if allAtomTypes[atom][:2]=="C3":
           ###Checks if atom is in cycle... assigns carbon and hydrogen types
           if isInCycle(nameToIdx[atom],cycle):
                newAtType = "C3R"
                finalAtType = createFinalAtomTypes(allAtomTypes,atom,newAtType)
                specialAtomTypes[atom] = finalAtType
                for at in bonds[atom]:
                    if allAtomTypes[at][0]=="H":
                        newAtType = "HC3R"
                        finalAtType = createFinalAtomTypes(allAtomTypes,at,newAtType)
                        specialAtomTypes[at] = finalAtType
           else:# assigns remaining atom types  
                for at in bonds[atom]:
                    if allAtomTypes[at][0]=="H":hCount+=1
                    if allAtomTypes[at][:2]=="C3" and at not in bonds: ### carbon bonded to a carbon and not having its bonds listed in the rtf 
                        newAtType = "C30" ### according to GAFF, must be a carbon without hydrogen
                        finalAtType = createFinalAtomTypes(allAtomTypes,at,newAtType)
                        specialAtomTypes[at] = finalAtType
                if hCount==0:
                    newAtType = "C30"
                    finalAtType = createFinalAtomTypes(allAtomTypes,atom,newAtType)
                    specialAtomTypes[atom]= finalAtType ### preserve the "_" from Lei!! 

                elif hCount==1:
                    newAtType = "C31"
                    finalAtType = createFinalAtomTypes(allAtomTypes,atom,newAtType)
                    specialAtomTypes[atom]= finalAtType
                    for at in bonds[atom]:
                        if allAtomTypes[at][:2]=="HC":
                           newAtType = "HC31"
                           finalAtType = createFinalAtomTypes(allAtomTypes,at,newAtType)
                           specialAtomTypes[at] = finalAtType
                elif hCount==2:
                     newAtType = "C32"
                     finalAtType = createFinalAtomTypes(allAtomTypes,atom,newAtType)
                     specialAtomTypes[atom] = finalAtType
                     for at in bonds[atom]:
                        if allAtomTypes[at][:2]=="HC":
                            newAtType = "HC32"
                            finalAtType = createFinalAtomTypes(allAtomTypes,at,newAtType)
                            specialAtomTypes[at] = finalAtType
                elif hCount==3:
                     newAtType = "C33"
                     finalAtType = createFinalAtomTypes(allAtomTypes,atom,newAtType)
                     specialAtomTypes[atom] = finalAtType
                     for at in bonds[atom]:
                         if allAtomTypes[at][:2]=="HC":
                             newAtType = "HC33"
                             finalAtType = createFinalAtomTypes(allAtomTypes,at,newAtType)
                             specialAtomTypes[at]= finalAtType
        elif allAtomTypes[atom][:2]=="CA":
            for at in bonds[atom]:
                if allAtomTypes[at][:2]=="CL":
                    newAtType = "CLA"
                    finalAtType = createFinalAtomTypes(allAtomTypes,at,newAtType)
                    specialAtomTypes[at]= finalAtType
                elif allAtomTypes[at][:2]=="OH":
                    newAtType = "OHP"
                    finalAtType = createFinalAtomTypes(allAtomTypes,at,newAtType)
                    specialAtomTypes[at]= finalAtType
        else: continue
    return specialAtomTypes



def combineAtomTypes(specialTypes,atoms,allTypes):
    allTypes = {atoms[i]:caps(allTypes[i]) for i in range(len(atoms))}
    atomTypes={}
    for a in allTypes:
	if specialTypes.has_key(a):atomTypes[a]=specialTypes[a]
	else:
	    try: 
		aType = allTypes[a].split('_')[0]+allTypes[a].split('_')[1]	 
	    	atomTypes[a] = aType
	    except IndexError: atomTypes[a] = allTypes[a]	
    return atomTypes


def caps(atName):
    return ''.join([i.capitalize() for i in atName])


def updateNonBonded(mapD1_D2,nbPrm,useOptLJ,isAdditive):
    newNBPrms={}	

    for atom in mapD1_D2:
	if useOptLJ:
	   oldType,newType = mapD1_D2[atom][0],mapD1_D2[atom][1]   
	   if newType[0]=="C" and len(newType)>=3: 
		newGPType = newType[:3]
	   elif newType[0]=="H" and len(newType)>=4: 
		newGPType = newType[:4]	
	   else: newGPType = newType	
	   if nbPrm.has_key(newGPType): 
		param = nbPrm[newGPType]
		newNBPrms[oldType] = param
           else: # if we are here, some expanded type with "_" scheme weren't found here and will be handeled
               if len(newGPType)==3: ### examples like CL0 and BR0
                    param = nbPrm[newGPType[:2]]
                    newNBPrms[oldType] = param
               elif len(newGPType)==2 and newGPType!="LP": ### examples like F0, O0, N0, and C0
                    param = nbPrm[newGPType[:1]]
                    newNBPrms[oldType] = param 
	else: 
	   oldType,newType = mapD1_D2[atom][0],mapD1_D2[atom][1]   
	   if nbPrm.has_key(oldType): 
		param = nbPrm[oldType]
		newNBPrms[oldType] = param
    if not isAdditive:
	#newNBPrms["D*"] = ["0.0","-0.0000","0.0000"]
	newNBPrms["DRUD"] = ["0.0","-0.0000","0.0000"]
	newNBPrms["LP"] = ["0.0","-0.0000","0.0000"]

    #else:
    #    newNBPrms["HT_W"] = ["0.0","-0.046000","0.224500"]    
    #    newNBPrms["OT_W"] = ["0.0","-0.152100","1.768200"]    
        
    return newNBPrms	


def readOptLJ(ljFile):
    lj={}
    fin=open(ljFile,"r")	
    for line in fin:
	token = line.rstrip().split()
	lj[token[0]] = ["0.00",token[1],token[2]]
    return lj	


def writeUpdatedPrm(prmFile,newNBPrm,isAdditive):
    fout = open("tmp.prm","w")
    fin = open(prmFile,"r")
    flag=0
    for line in fin:
        token = line.rstrip().split()
        try:
            if token[0]=="BONDS" or token[0]=="ANGLES" or token[0]=="DIHEDRALS" or token[0]=="IMPROPERS":
                print>>fout,"\n"
            if token[0]=="NONBONDED":flag=1
        except IndexError:continue
        if flag==0: print>>fout,line.rstrip();

    print>>fout,"\n","NONBONDED  E14FAC  1.000000"
    print>>fout,"!               EMIN   RMIN/2        EMIN/2   RMIN(FOR 1-4'S)"
    print>>fout,"!           (KCAL/MOL)  (A)"
    for nb in newNBPrm:
	z = float(newNBPrm[nb][0]); e= float(newNBPrm[nb][1]); r = float(newNBPrm[nb][2])
	fout.write('{:6} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}\n'.format(nb,z,e,r,z,float(e)/2,r))
   
    
    fout.write("\nEND")	
    fout.close()	


def mapDef1ToDef2(atomsDef1,allAtomTypesDef1,allAtomTypesDef2):
    mapAtom_Types_Def1 = {k:v for k,v in zip(atomsDef1,allAtomTypesDef1)}
    mapD1_D2={}
    for k in mapAtom_Types_Def1:
	mapD1_D2[k] = [mapAtom_Types_Def1[k],allAtomTypesDef2[k]]
    return mapD1_D2	


def run(mol2,ac_rtf,ac_prm,useOptLJ,isAdditive):

    cycle = detectAtomsInCycle(mol2) 

    nameToIdx, bondsDef1,angleDef1,dihedsDef1,impropsDef1 = getConnections(ac_rtf)
    atomsDef1,allAtomTypesDef1 = getAllAtomTypes(ac_rtf) ## different rtf b/c gaamp assigns "_" to expand some GAFF types 

    specialAtomTypesDef2 = expandAtomTypes(nameToIdx,cycle,bondsDef1,atomsDef1,allAtomTypesDef1)

    allAtomTypesDef2 = combineAtomTypes(specialAtomTypesDef2,atomsDef1,allAtomTypesDef1) 

    mapD1_D2 = mapDef1ToDef2(atomsDef1,allAtomTypesDef1,allAtomTypesDef2)    

    if useOptLJ: 
	if isAdditive:
    	   ljOptPrm = readOptLJ(GAAMPDIR+"/scripts/LJPARAMETERS.dat")
    	   newNBPrm = updateNonBonded(mapD1_D2,ljOptPrm,useOptLJ,isAdditive)	
    else:
	newNBPrm = updateNonBonded(mapD1_D2,gaffNBPrm,useOptLJ,isAdditive)

    writeUpdatedPrm(ac_prm,newNBPrm,isAdditive)

    os.system("mv tmp.prm mol.prm")


if __name__=="__main__":
   
    global GAAMPDIR
    GAAMPDIR = sys.argv[1]
  
    os.chdir("020-initial_parameters")
    os.system("cp mol.prm ac-mol.prm")

    setBabelPath()
    sys.path.append(BABELPATH+"/lib64/")
    import openbabel

    mol2 = sys.argv[2] 
    ac_rtf = sys.argv[3]
    ac_prm = sys.argv[4]
    useOptLJ = True
    isAdditive = True
    run(mol2,ac_rtf,ac_prm,useOptLJ,isAdditive) 

