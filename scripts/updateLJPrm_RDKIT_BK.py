import os,sys
import itertools
from rdkit import Chem

GAAMPSCRIPTS=""


atomicMass = {"LP":0.00000,"DRUD":"0.00000","H":1.008000,"C":12.01000,"O":16.00000,"N":14.00670,"F":18.998403,"S":32.06500,"CL":35.45300,"P":30.97376,"BR":79.90400}


#def setRDKitEnv():
#	os.system("source activate my-rdkit-env")
#	from rdkit import Chem

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
    """	
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol2", "mol2")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol,mol2)
    cycle = []
    for obatom in openbabel.OBMolAtomIter(mol):
        if obatom.IsInRing():
            cycle.append(obatom.GetIndex())
    """
    cycle = []	  
    m = Chem.MolFromMol2File(mol2)
    for i in range(m.GetNumAtoms()):
	if m.GetAtomWithIdx(i).IsInRing():
		cycle.append(i)			
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
    
    #print mapD1_D2

    for atom in mapD1_D2:
	if useOptLJ:
	   oldType,newType = mapD1_D2[atom][0],mapD1_D2[atom][1]   
	   if nbPrm.has_key(newType): 
		param = nbPrm[newType]
		newNBPrms[newType] = param
	   else: 
		param = nbPrm[newType[:3]]
		newNBPrms[newType] = param 
        
    return newNBPrms	


def readOptLJ(ljFile):
    lj={}
    fin=open(ljFile,"r")	
    for line in fin:
	token = line.rstrip().split()
	lj[token[0]] = ["0.00",token[1],token[2]]
    return lj	


def updateRTF(ac_rtf,atoms,allAtomTypes,isAdditive):
    atomTypes = {v:k for k,v in allAtomTypes.iteritems()}
    fout = open("tmp.rtf","w")

    totalCharge = [c.rstrip().split()[2] for c in open(ac_rtf,"r") if len(c.rstrip().split())>0 \
                and c.rstrip().split()[0]=="RESI" and c.rstrip().split()[1]=="MOL"][0]
    pCharge = {}
    alphaThole = {}
    fin = open(ac_rtf,"r")
    flag=0
    for line in fin:
        token = line.rstrip().split()
        if len(token)>0:
           if token[0]=="GROUP":flag+=1
           if flag==1 and token[0]=="ATOM":
                atom = token[1]
                pC = token[3]
                pCharge[atom]=pC
                try:
                    alpha = token[5]
                    thole = token[7]
                    alphaThole[atom] = [alpha,thole]
                except IndexError:continue
    fin.close()

    print>>fout,"* Topology File."
    print>>fout,"*"
    print>>fout,"   99   1"
    for a in atomTypes:
        if caps(a[:2])=="CL":mass = atomicMass["CL"]
        elif caps(a[:2])=="BR":mass = atomicMass["BR"]
        elif caps(a[:2])=="LP": mass = atomicMass["LP"]
        elif caps(a[0])=="C": mass = atomicMass["C"]
        elif caps(a[0])=="H": mass = atomicMass["H"]
        elif caps(a[0])=="O": mass = atomicMass["O"]
        else: mass = atomicMass[caps(a[0])]
        fout.write('{:9} {:4} {:8} {:10}\n'.format("MASS",-1,a,mass))
    if not isAdditive:fout.write('{:9} {:4} {:8} {:10}\n'.format("MASS",-1,"DRUD",0.0000))

    print >>fout,"\n"
    if not isAdditive: print>>fout, "AUTOGENERATE ANGLES DIHEDRALS DRUDE\n"
    else: print>>fout, "AUTOGENERATE ANGLES DIHEDRALS\n"

    print>>fout,"\n","RESI MOL ",totalCharge
    print>>fout,"GROUP"
    for a in atoms:
        #fout.write('{:5} {:6} {:7} {:10}\n'.format("ATOM",a,allAtomTypes[a],pCharge[a]))
        try:
            alpha,thole = alphaThole[a][0],alphaThole[a][1]
            fout.write('{:5} {:6} {:7} {:10} {:8} {:8} {:8} {:8}\n'.format("ATOM",a,allAtomTypes[a],pCharge[a],"ALPHA",alpha,"THOLE",thole))
        except KeyError:
            fout.write('{:5} {:6} {:7} {:10}\n'.format("ATOM",a,allAtomTypes[a],pCharge[a]))

    print>>fout,"\n"

    fin=open(ac_rtf,"r")
    for line in fin:
        token=line.rstrip().split()
        if len(token)>1:
           if token[0]=="RESI" and (token[1]=="SWM4" or token[1]=="TIP3"):break
           else:
               if token[0]=="BOND":print>>fout,line.rstrip() #and token[1]!="OH2":print>>fout,line.rstrip()
               elif token[0]=="LONEPAIR" or token[0]=="ANISOTROPY": print>>fout,line.rstrip()
    print>>fout,"\n"
    print>>fout,"END"
    fout.close()
    fin.close()
    os.system("mv tmp.rtf mol.rtf")


def writeUpdatedPrm(newBondPrm,newAnglePrm,newDihedPrm,newImpPrm,newNBPrm,isAdditive):
    fout = open("tmp.prm","w")
    print>>fout,"* FORCE FIELD PARAMETER FILE."
    print>>fout,"*","\n"
    print>>fout,"BONDS"
    for b in newBondPrm:
        b1 = b.split("-")[0];b2 = b.split("-")[1]
        k = float(newBondPrm[b][0]);b0= float(newBondPrm[b][1])
        fout.write('{:6} {:6} {:.4f} {:.4f}\n'.format(b1,b2,k,b0))

    if not isAdditive:fout.write('{:6} {:6} {:.4f} {:.4f}\n'.format("X","DRUD",500.00,0.000)) ### bonds params for Drude 

    print>>fout,"\n","ANGLES"
    for a in newAnglePrm:
        a1 = a.split("-")[0];a2 = a.split("-")[1]; a3 = a.split("-")[2]
        k = float(newAnglePrm[a][0]);a0 = float(newAnglePrm[a][1])
        fout.write('{:6} {:6} {:6} {:.4f} {:.4f}\n'.format(a1,a2,a3,k,a0))

    print>>fout,"\n","DIHEDRALS"
    for d in newDihedPrm:
        d1 = d.split("-")[0];d2 = d.split("-")[1];d3 = d.split("-")[2];d4 = d.split("-")[3]
        for val in newDihedPrm[d]:
            k = float(val[0]); n = float(val[1]); pa = float(val[2])
            fout.write('{:6} {:6} {:6} {:6} {:.4f} {:.3f} {:.4f}\n'.format(d1,d2,d3,d4,k,n,pa))

    print>>fout,"\n","IMPROPERS"
    for im in newImpPrm:
        im1 = im.split("-")[0];im2 = im.split("-")[1];im3 = im.split("-")[2];im4 = im.split("-")[3]
        for val in newImpPrm[im]:
            k = float(val[0]); n = float(val[1]); pa = float(val[2])
            fout.write('{:6} {:6} {:6} {:6} {:.4f} {:.3f} {:.4f}\n'.format(im1,im2,im3,im4,k,n,pa))

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
	mapD1_D2[caps(k)] = [mapAtom_Types_Def1[k],allAtomTypesDef2[k]]
    return mapD1_D2	


def readPRM(prm):
    fin=open(prm,"r")
    bondFlag=0;angleFlag=0;dihedFlag=0;impropFlag=0;nbFlag=0;
    bondRecord={};angleRecord={};dihedRecord={};impropRecord={};nbRecord={};
    for line in fin:
        token=line.rstrip().split()
        if len(token)>0:
           if token[0][:4]=="BOND":
              bondFlag=1
           elif token[0][:5]=="ANGLE":
              bondFlag=0;angleFlag=1
           elif token[0][:5]=="DIHED":
              angleFlag=0;dihedFlag=1
           elif token[0][:6]=="IMPROP" or token[0][:5]=="IMPHI":
              dihedFlag=0;impropFlag=1
           elif token[0]=="NONBONDED":
              dihedFlag=0;impropFlag=0;nbFlag=1
           if bondFlag==1 and token[0][:4]!="BOND" :
              bondRecord[token[0]+"-"+token[1]] = [token[2],token[3]]
           if angleFlag==1 and token[0][:5]!="ANGLE":
              angleRecord[token[0]+"-"+token[1]+"-"+token[2]] = [token[3],token[4]]
           if dihedFlag==1 and token[0][:5]!="DIHED":
              a1 = token[0]; a2 = token[1]; a3 = token[2]; a4 = token[3]
              if not dihedRecord.has_key(a1+"-"+a2+"-"+a3+"-"+a4):
                 dihedRecord[a1+"-"+a2+"-"+a3+"-"+a4] = [[token[4],token[5],token[6]]]
              else:dihedRecord[a1+"-"+a2+"-"+a3+"-"+a4].append([token[4],token[5],token[6]])
           if impropFlag==1 and (token[0][:6]!="IMPROP" or token[0][:5]=="IMPHI"):
              try:
                  a1 = token[0]; a2 = token[1]; a3 = token[2]; a4 = token[3]
                  impropRecord[a1+"-"+a2+"-"+a3+"-"+a4] = [[token[4],token[5],token[6]]]
              except IndexError:continue

           try:
                if nbFlag==1 and token[1]!="(KCAL/MOL)":
                   if token[0]!="NONBONDED" and token[0]!="!" and token[0]!="CUTNB":
                        if token[0]=="HT_W" or token[1]=="OT_W":nbRecord[token[0]] = [token[1],token[2],token[3]]
                        else:
                                a = token[0]
                                nbRecord[a] = [token[1],token[2],token[3]]
           except IndexError:continue
    fin.close()
    return bondRecord,angleRecord,dihedRecord,impropRecord,nbRecord

def findNewType(oldAtomType,mapD1_D2):
    newTypes=[]
    for atom in mapD1_D2:
        if mapD1_D2[atom][0]==oldAtomType:
           newTypes.append(mapD1_D2[atom][1])
    return newTypes


def updateBondPrm(orgBond,mapD1_D2,bondPrm,isAdditive):
    newBondParams={}
    for b in orgBond:
        for a in orgBond[b]:
            obt,nbt = mapD1_D2[b][0],mapD1_D2[b][1]
            oat,nat = mapD1_D2[a][0],mapD1_D2[a][1]
            if bondPrm.has_key(obt+"-"+oat):
                k = bondPrm[obt+"-"+oat][0]; b0 = bondPrm[obt+"-"+oat][1]
                if not newBondParams.has_key(nbt+"-"+nat) and not newBondParams.has_key(nat+"-"+nbt):
                   newBondParams[nbt+"-"+nat] = [k,b0]
            elif bondPrm.has_key(oat+"-"+obt):
                k = bondPrm[oat+"-"+obt][0]; b0 = bondPrm[oat+"-"+obt][1]
                if not newBondParams.has_key(nbt+"-"+nat) and not newBondParams.has_key(nat+"-"+nbt):
                   newBondParams[nbt+"-"+nat] = [k,b0]
            else:print "FATAL ERROR: Missing Bond paramters for "+obt+"-"+oat+" !!";sys.exit()
    """	
    if not isAdditive: ### adding LP parameters if any
        for atom in lpBondedAtoms:
            newAtomType = mapD1_D2[atom][1]
            k,b0 = "0.00","0.000"
            newBondParams[newAtomType+"-"+"LP"] = [k,b0]
    """	
    return newBondParams


def updateAnglePrm(orgAngle,mapD1_D2,anglePrm):

    #print orgAngle
    #print mapD1_D2

    newAnglePrms={}
    for angle in orgAngle:
        oldAtomType1,newAtomType1 = mapD1_D2[angle[0]][0],mapD1_D2[angle[0]][1]
        oldAtomType2,newAtomType2 = mapD1_D2[angle[1]][0],mapD1_D2[angle[1]][1]
        oldAtomType3,newAtomType3 = mapD1_D2[angle[2]][0],mapD1_D2[angle[2]][1]
        oldAngleType =  oldAtomType1+"-"+oldAtomType2+"-"+oldAtomType3
        oldAngleTypeRev =  oldAtomType3+"-"+oldAtomType2+"-"+oldAtomType1
        newAngleType = newAtomType1+"-"+newAtomType2+"-"+newAtomType3
        if anglePrm.has_key(oldAngleType): newAnglePrm = anglePrm[oldAngleType]
        elif anglePrm.has_key(oldAngleTypeRev): newAnglePrm = anglePrm[oldAngleTypeRev]
        else:print "FATAL ERROR: Missing Angle paramters for "+oldAngleType+" !!";sys.exit()
        newAnglePrms[newAngleType] = newAnglePrm
    return newAnglePrms


def updateDihedPrm(orgDihed,mapD1_D2,dihedPrm):
    newDihedPrms={}
    for dihed in orgDihed:
        oldAtomType1,newAtomType1 = mapD1_D2[dihed[0]][0],mapD1_D2[dihed[0]][1]
        oldAtomType2,newAtomType2 = mapD1_D2[dihed[1]][0],mapD1_D2[dihed[1]][1]
        oldAtomType3,newAtomType3 = mapD1_D2[dihed[2]][0],mapD1_D2[dihed[2]][1]
        oldAtomType4,newAtomType4 = mapD1_D2[dihed[3]][0],mapD1_D2[dihed[3]][1]
        oldDihedType =  oldAtomType1+"-"+oldAtomType2+"-"+oldAtomType3+"-"+oldAtomType4
        oldDihedTypeRev =  oldAtomType4+"-"+oldAtomType3+"-"+oldAtomType2+"-"+oldAtomType1
        newDihedType = newAtomType1+"-"+newAtomType2+"-"+newAtomType3+"-"+newAtomType4
        newDihedTypeRev = newAtomType4+"-"+newAtomType3+"-"+newAtomType2+"-"+newAtomType1
        if dihedPrm.has_key(oldDihedType): newDihedPrm = dihedPrm[oldDihedType]
        elif dihedPrm.has_key(oldDihedTypeRev):newDihedPrm = dihedPrm[oldDihedTypeRev]
        elif dihedPrm.has_key("X-"+oldAtomType2+"-"+oldAtomType3+"-X"):
                newDihedPrm = dihedPrm["X-"+oldAtomType2+"-"+oldAtomType3+"-X"]
        elif dihedPrm.has_key("X-"+oldAtomType3+"-"+oldAtomType2+"-X"):
                newDihedPrm = dihedPrm["X-"+oldAtomType3+"-"+oldAtomType2+"-X"]
        else:
            print "FATAL ERROR: Missing Dihedral paramters for "+newDihedType
            print "FATAL ERROR: Missing Dihedral paramters for "+oldDihedType+ "!!";sys.exit()

        if newDihedPrms.has_key(newDihedType) or newDihedPrms.has_key(newDihedTypeRev):
            if newDihedPrms.has_key(newDihedType):
                prms = newDihedPrms[newDihedType]
                flag = 0
                for p in prms:
                    if p[0]==newDihedPrm[0][0] and p[1]==newDihedPrm[0][1] and p[2]==newDihedPrm[0][2]:
                        flag = 1
                if flag==0: newDihedPrms[newDihedType].append(newDihedPrm[0])
            else:
                prms = newDihedPrms[newDihedTypeRev]
                flag = 0
                for p in prms:
                    if p[0]==newDihedPrm[0][0] and p[1]==newDihedPrm[0][1] and p[2]==newDihedPrm[0][2]:
                        flag = 1
                if flag==0: newDihedPrms[newDihedTypeRev].append(newDihedPrm[0])

        else: newDihedPrms[newDihedType] = newDihedPrm

    ###Lets add all the "X" records here. If an atom type is expanded, include all types
    for p in dihedPrm:
        atomTypes = p.split("-")
        if atomTypes[0]=="X" or atomTypes[1]=="X" or atomTypes[2]=="X" or atomTypes[3]=="X":
           if atomTypes[0]=="X": newType1 = "X"
           else: newType1 = findNewType(atomTypes[0],mapD1_D2)
           if atomTypes[1]=="X": newType2 = "X"
           else: newType2 = findNewType(atomTypes[1],mapD1_D2)
           if atomTypes[2]=="X": newType3 = "X"
           else: newType3 = findNewType(atomTypes[2],mapD1_D2)
           if atomTypes[3]=="X": newType4 = "X"
           else: newType4 = findNewType(atomTypes[3],mapD1_D2)
           for combination in itertools.product(newType1,newType2,newType3,newType4):
                newDihedType = combination[0]+"-"+combination[1]+"-"+combination[2]+"-"+combination[3]
                newDihedTypeRev = combination[3]+"-"+combination[2]+"-"+combination[1]+"-"+combination[0]
                oldDihedType = atomTypes[0]+"-"+atomTypes[1]+"-"+atomTypes[2]+"-"+atomTypes[3]
                if dihedPrm.has_key(oldDihedType): newDihedPrm = dihedPrm[oldDihedType]
                else:
                    print "FATAL ERROR: Missing Dihedral paramters for "+newDihedType
                    print "FATAL ERROR: Missing Dihedral paramters for "+oldDihedType+ "!!";sys.exit()

                #if newDihedPrms.has_key(newDihedType) or newDihedPrms.has_key(newDihedTypeRev):continue
                #else: newDihedPrms[newDihedType] = newDihedPrm

                if newDihedPrms.has_key(newDihedType) or newDihedPrms.has_key(newDihedTypeRev):
                    if newDihedPrms.has_key(newDihedType):
                        prms = newDihedPrms[newDihedType]
                        flag = 0
                        for p in prms:
                            if p[0]==newDihedPrm[0][0] and p[1]==newDihedPrm[0][1] and p[2]==newDihedPrm[0][2]:
                                flag = 1
                        if flag==0: newDihedPrms[newDihedType].append(newDihedPrm[0])
                    else:
                        prms = newDihedPrms[newDihedTypeRev]
                        flag = 0
                        for p in prms:
                            if p[0]==newDihedPrm[0][0] and p[1]==newDihedPrm[0][1] and p[2]==newDihedPrm[0][2]:
                                flag = 1
                        if flag==0:
                            newDihedPrms[newDihedTypeRev].append(newDihedPrm[0])

                else: newDihedPrms[newDihedType] = newDihedPrm

                #newDihedPrms[newDihedType] = newDihedPrm

    return newDihedPrms


def updateImpropPrm(orgImprop,mapD1_D2,impropPrm):
    newImpropPrms={}
    for improp in orgImprop:
        oldAtomType1,newAtomType1 = mapD1_D2[improp[0]][0],mapD1_D2[improp[0]][1]
        oldAtomType2,newAtomType2 = mapD1_D2[improp[1]][0],mapD1_D2[improp[1]][1]
        oldAtomType3,newAtomType3 = mapD1_D2[improp[2]][0],mapD1_D2[improp[2]][1]
        oldAtomType4,newAtomType4 = mapD1_D2[improp[3]][0],mapD1_D2[improp[3]][1]
        oldImpropType =  oldAtomType1+"-"+oldAtomType2+"-"+oldAtomType3+"-"+oldAtomType4
        oldImpropTypeRev =  oldAtomType4+"-"+oldAtomType3+"-"+oldAtomType2+"-"+oldAtomType1
        newImpropType = newAtomType1+"-"+newAtomType2+"-"+newAtomType3+"-"+newAtomType4
        newImpropTypeRev = newAtomType4+"-"+newAtomType3+"-"+newAtomType2+"-"+newAtomType1
        if impropPrm.has_key(oldImpropType): newImpropPrm = impropPrm[oldImpropType]
        elif impropPrm.has_key(oldImpropTypeRev):newImpropPrm = impropPrm[oldImpropTypeRev]
        elif impropPrm.has_key("X-"+oldAtomType2+"-"+oldAtomType3+"-X"):
                newImpropPrm = impropPrm["X-"+oldAtomType2+"-"+oldAtomType3+"-X"]
        elif impropPrm.has_key("X-"+oldAtomType3+"-"+oldAtomType2+"-X"):
                newImpropPrm = impropPrm["X-"+oldAtomType3+"-"+oldAtomType2+"-X"]
        elif impropPrm.has_key("X-X-"+oldAtomType3+"-"+oldAtomType4):
                newImpropPrm = impropPrm["X-X-"+oldAtomType3+"-"+oldAtomType4]
        elif impropPrm.has_key(oldAtomType1+"-"+oldAtomType2+"-X-X"):
                newImpropPrm = impropPrm[oldAtomType1+"-"+oldAtomType2+"-X-X"]
        else:print "FATAL ERROR: Missing Improper paramters for "+oldImpropType+" !!";sys.exit()

        #if not newImpropPrms.has_key(newImpropType) or not newImpropPrms.has_key(newImpropTypeRev):
        #    newImpropPrms[newImpropType] = newImpropPrm

        if newImpropPrms.has_key(newImpropType) or newImpropPrms.has_key(newImpropTypeRev):
            if newImpropPrms.has_key(newImpropType):
                prms = newImpropPrms[newImpropType]
                flag = 0
                for p in prms:
                    if p[0]==newImpropPrm[0][0] and p[1]==newImpropPrm[0][1] and p[2]==newImpropPrm[0][2]:
                        flag = 1
                if flag==0: newImpropPrms[newImpropType].append(newImpropPrm[0])
            else:
                prms = newImpropPrms[newImpropTypeRev]
                flag = 0
                for p in prms:
                    if p[0]==newImpropPrm[0][0] and p[1]==newImpropPrm[0][1] and p[2]==newImpropPrm[0][2]:
                        flag = 1
                if flag==0: newImpropPrms[newImpropTypeRev].append(newImpropPrm[0])
        else: newImpropPrms[newImpropType] = newImpropPrm

    ###Lets add all the "X" records here. If an atom type is expanded, include all types
    for p in impropPrm:
        atomTypes = p.split("-")
        if atomTypes[0]=="X" or atomTypes[1]=="X" or atomTypes[2]=="X" or atomTypes[3]=="X":
           if atomTypes[0]=="X": newType1 = "X"
           else: newType1 = findNewType(atomTypes[0],mapD1_D2)
           if atomTypes[1]=="X": newType2 = "X"
           else: newType2 = findNewType(atomTypes[1],mapD1_D2)
           if atomTypes[2]=="X": newType3 = "X"
           else: newType3 = findNewType(atomTypes[2],mapD1_D2)
           if atomTypes[3]=="X": newType4 = "X"
           else: newType4 = findNewType(atomTypes[3],mapD1_D2)
           for combination in itertools.product(newType1,newType2,newType3,newType4):
                newImpropType = combination[0]+"-"+combination[1]+"-"+combination[2]+"-"+combination[3]
                newImpropTypeRev = combination[3]+"-"+combination[2]+"-"+combination[1]+"-"+combination[0]
                oldImpropType = atomTypes[0]+"-"+atomTypes[1]+"-"+atomTypes[2]+"-"+atomTypes[3]
                if impropPrm.has_key(oldImpropType): newImpropPrm = impropPrm[oldImpropType]
                else:
                    print "FATAL ERROR: Missing Dihedral paramters for "+newImpropType
                    print "FATAL ERROR: Missing Dihedral paramters for "+oldImpropType+ "!!";sys.exit()

                #if not newImpropPrms.has_key(newImpropType) or not newImpropPrms.has_key(newImpropTypeRev):
                #    newImpropPrms[newImpropType] = newImpropPrm

                if newImpropPrms.has_key(newImpropType) or newImpropPrms.has_key(newImpropTypeRev):
                    if newImpropPrms.has_key(newImpropType):
                        prms = newImpropPrms[newImpropType]
                        flag = 0
                        for p in prms:
                            if p[0]==newImpropPrm[0][0] and p[1]==newImpropPrm[0][1] and p[2]==newImpropPrm[0][2]:
                                flag = 1
                        if flag==0: newImpropPrms[newImpropType].append(newImpropPrm[0])
                    else:
                        prms = newImpropPrms[newImpropTypeRev]
                        flag = 0
                        for p in prms:
                            if p[0]==newImpropPrm[0][0] and p[1]==newImpropPrm[0][1] and p[2]==newImpropPrm[0][2]:
                                flag = 1
                        if flag==0: newImpropPrms[newImpropTypeRev].append(newImpropPrm[0])
                else: newImpropPrms[newImpropType] = newImpropPrm

                #newImpropPrms[newImpropType] = newImpropPrm

    return newImpropPrms


def writeUpdatedPrm(newBondPrm,newAnglePrm,newDihedPrm,newImpPrm,newNBPrm,isAdditive):
    fout = open("tmp.prm","w")
    print>>fout,"* FORCE FIELD PARAMETER FILE."
    print>>fout,"*","\n"
    print>>fout,"BONDS"
    for b in newBondPrm:
        b1 = b.split("-")[0];b2 = b.split("-")[1]
        k = float(newBondPrm[b][0]);b0= float(newBondPrm[b][1])
        fout.write('{:6} {:6} {:.4f} {:.4f}\n'.format(b1,b2,k,b0))

    if not isAdditive:fout.write('{:6} {:6} {:.4f} {:.4f}\n'.format("X","DRUD",500.00,0.000)) ### bonds params for Drude 

    print>>fout,"\n","ANGLES"
    for a in newAnglePrm:
        a1 = a.split("-")[0];a2 = a.split("-")[1]; a3 = a.split("-")[2]
        k = float(newAnglePrm[a][0]);a0 = float(newAnglePrm[a][1])
        fout.write('{:6} {:6} {:6} {:.4f} {:.4f}\n'.format(a1,a2,a3,k,a0))

    print>>fout,"\n","DIHEDRALS"
    for d in newDihedPrm:
        d1 = d.split("-")[0];d2 = d.split("-")[1];d3 = d.split("-")[2];d4 = d.split("-")[3]
        for val in newDihedPrm[d]:
            k = float(val[0]); n = float(val[1]); pa = float(val[2])
            fout.write('{:6} {:6} {:6} {:6} {:.4f} {:.3f} {:.4f}\n'.format(d1,d2,d3,d4,k,n,pa))

    print>>fout,"\n","IMPROPERS"
    for im in newImpPrm:
        im1 = im.split("-")[0];im2 = im.split("-")[1];im3 = im.split("-")[2];im4 = im.split("-")[3]
        for val in newImpPrm[im]:
            k = float(val[0]); n = float(val[1]); pa = float(val[2])
            fout.write('{:6} {:6} {:6} {:6} {:.4f} {:.3f} {:.4f}\n'.format(im1,im2,im3,im4,k,n,pa))

    print>>fout,"\n","NONBONDED  E14FAC  1.000000"
    print>>fout,"!               EMIN   RMIN/2        EMIN/2   RMIN(FOR 1-4'S)"
    print>>fout,"!           (KCAL/MOL)  (A)"
    for nb in newNBPrm:
        z = float(newNBPrm[nb][0]); e= float(newNBPrm[nb][1]); r = float(newNBPrm[nb][2])
        fout.write('{:6} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}\n'.format(nb,z,e,r,z,float(e)/2,r))

    fout.write("\nEND")
    fout.close()



def run(mol2,ac_rtf,ac_prm,useOptLJ,isAdditive):

    cycle = detectAtomsInCycle(mol2) 

    nameToIdx, bondsDef1,angleDef1,dihedsDef1,impropsDef1 = getConnections(ac_rtf)
    atomsDef1,allAtomTypesDef1 = getAllAtomTypes(ac_rtf) ## different rtf b/c gaamp assigns "_" to expand some GAFF types 

    specialAtomTypesDef2 = expandAtomTypes(nameToIdx,cycle,bondsDef1,atomsDef1,allAtomTypesDef1)

    allAtomTypesDef2 = combineAtomTypes(specialAtomTypesDef2,atomsDef1,allAtomTypesDef1) 

    mapD1_D2 = mapDef1ToDef2(atomsDef1,allAtomTypesDef1,allAtomTypesDef2)    

    updateRTF(ac_rtf,atomsDef1,allAtomTypesDef2,isAdditive)	

    bondPrm,anglePrm,dihedPrm,impropPrm,gaffNBPrm = readPRM(ac_prm)
    newBondPrm = updateBondPrm(bondsDef1,mapD1_D2,bondPrm,isAdditive)
    newAnglePrm = updateAnglePrm(angleDef1,mapD1_D2,anglePrm)
	
    newDihedPrm = updateDihedPrm(dihedsDef1,mapD1_D2,dihedPrm)
    newImpPrm = updateImpropPrm(impropsDef1,mapD1_D2,impropPrm)

    if useOptLJ: 
	if isAdditive:
    	   ljOptPrm = readOptLJ(GAAMPSCRIPTS+"/LJPARAMETERS.dat")
    	   newNBPrm = updateNonBonded(mapD1_D2,ljOptPrm,useOptLJ,isAdditive)	
	   #print "LJs..."
	   #print newNBPrm
    else:
	newNBPrm = updateNonBonded(mapD1_D2,gaffNBPrm,useOptLJ,isAdditive)

    writeUpdatedPrm(newBondPrm,newAnglePrm,newDihedPrm,newImpPrm,newNBPrm,isAdditive)	
    #writeUpdatedPrm(ac_prm,newNBPrm,isAdditive)

    os.system("mv tmp.prm mol.prm")


if __name__=="__main__":
   
    global GAAMPSCRIPTS
    GAAMPSCRIPTS = sys.argv[1] 

    os.chdir("020-initial_parameters")
 	
    os.system("cp mol.prm ac-mol.prm")
    os.system("cp mol.rtf ac-mol.rtf")

    #setRDKitEnv()	
	
    mol2 = sys.argv[2] 
    ac_rtf = sys.argv[3]
    ac_prm = sys.argv[4]
    useOptLJ = True
    isAdditive = True
    run(mol2,ac_rtf,ac_prm,useOptLJ,isAdditive) 

