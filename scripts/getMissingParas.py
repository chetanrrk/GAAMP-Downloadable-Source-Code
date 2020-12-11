import os,sys

class Master:
    def __init__(self,bonds,angles,diheds,improps,ljs):	
    	self.bonds = bonds
    	self.angles = angles
    	self.diheds = diheds
    	self.improps = improps
	self.ljs = ljs

def readPara(prmF):
    fin = open(prmF,"r")
    bonds = {}	
    angles = {}
    diheds = {}	
    improps = {}
    ljs = {}
    bondF = 1; angleF=1; dihedF = 1; impropsF = 1; nonBond = 1	
    for line in fin:
	token = line.rstrip().split()
	try:
	    if token[0]=="BONDS":
		bondF = 0
	    elif token[0]=="ANGLES":
		angleF = 0
		bondF = 1
	    elif token[0]=="DIHEDRALS":
		dihedF = 0
		angleF = 1
	    elif token[0]=="IMPROPERS":
		impropsF = 0
		dihedF = 1		
	    elif (token[0]=="HGA1" and token[1]=="0.0") or (token[0]=="NONBONDED"):
		nonBond = 0
		impropsF = 1
	    elif token[0]=="NBFIX" or token[0]=="END":
		nbFix = 0
		nonBond = 1

	    if bondF ==0 and token[:4]!="BOND":				
		bonds[token[0]+"-"+token[1]] = [token[2],token[3]]
	    elif angleF==0 and token[0]!="ANGLES":
		if len(token)==5:
			angles[token[0]+"-"+token[1]+"-"+token[2]] = [token[3],token[4]]
		elif len(token)>5 and token[5][0]!="!": 
			angles[token[0]+"-"+token[1]+"-"+token[2]] = [token[3],token[4],token[5],token[6]]
		else: angles[token[0]+"-"+token[1]+"-"+token[2]] = [token[3],token[4]] 
	    elif dihedF==0 and token[0]!="DIHEDRALS":
		if not diheds.has_key(token[0]+"-"+token[1]+"-"+token[2]+"-"+token[3]):
			diheds[token[0]+"-"+token[1]+"-"+token[2]+"-"+token[3]] = [[token[4],token[5],token[6]]]		
		else: diheds[token[0]+"-"+token[1]+"-"+token[2]+"-"+token[3]].append([token[4],token[5],token[6]])	
	    elif impropsF==0 and token[0]!="IMPROPERS":
		improps[token[0]+"-"+token[1]+"-"+token[2]+"-"+token[3]] = [token[4],token[5],token[6]] 		
	    elif nonBond==0 and token[0]!="NONBONDED": 
		ljs[token[0]] =  [token[2],token[3]]		
		
	except IndexError: continue

    fin.close()		
    
    return bonds, angles, diheds, improps, ljs

def readAtomTypeMap(aMap):
    atomMap = {}
    fin = open(aMap,"r")
    for line in fin:
	try:
	    token = line.rstrip().split()
	    if token[0]=="CGenFF": continue	  
	    atomMap[token[0]] = token[1]	
	except IndexError: continue		
    fin.close()
    return atomMap


def getParasFromMaster(params,master,isDihed=True):

    vals = {}	
    if len(params)<1: return vals
    if len(params.keys()[0].split("-"))==2: ## check bonds; both ways
	for ps in params:
	    p = ps.split("-")	
	    ab = p[0]+"-"+p[1]	
	    ba = p[1]+"-"+p[0]
	    if master.bonds.has_key(ab): vals[ab] = master.bonds[ab]
	    elif master.bonds.has_key(ba): vals[ba] = master.bonds[ba]	   	
	    else: vals[ps] = params[ps] ### cgennff guessed param	
    elif len(params.keys()[0].split("-")) == 3: ## check angles; both ways
	for ps in params:
	    p = ps.split("-")	
	    abc = p[0]+"-"+p[1]+"-"+p[2]
	    cba = p[2]+"-"+p[1]+"-"+p[0]
	    if master.angles.has_key(abc): 
		vals[abc] = master.angles[abc]
	    elif master.angles.has_key(cba): 
		vals[cba] = master.angles[cba]
	    else: vals[ps] = params[ps]	### cgenff guessed param
    elif len(params.keys()[0].split("-")) == 4: ## check diheds/improps; both ways
	for ps in params:
	    p = ps.split("-")	
	    abcd = p[0]+"-"+p[1]+"-"+p[2]+"-"+p[3]	
	    dcba = p[3]+"-"+p[2]+"-"+p[1]+"-"+p[0]	
	    if isDihed:	
	    	if master.diheds.has_key(abcd): vals[abcd] = master.diheds[abcd]
	    	elif master.diheds.has_key(dcba): vals[dcba] = master.diheds[dcba]
		else: vals[ps] = params[ps]   ### cgenff guessed param
	    else:	
	    	if master.improps.has_key(abcd): vals[abcd] = master.improps[abcd]
	    	elif master.improps.has_key(dcba): vals[dcba] = master.improps[dcba]
		else: vals[ps] = params[ps]   ### cgenff guessed param
    elif len(params.keys()[0].split("-"))==1:
	for ps in params:
	    vals[ps] = master.ljs[ps]
    return vals	

def readAlexLJ():
    fin = open("/lcrc/project/Drude/chetan/OLD_GAAMP/cgenff_gaamp/scripts/drude_lj_alex_2018.txt","r")
    ljs = {}
    for line in fin:
	token=line.rstrip().split()
	if token[0]!="Type":
		ljs[token[0]] = [token[1],token[2]]
    fin.close()
    return ljs

def getAtomTypes():
    fin = open("mol.rtf","r")	
    atomTypes = []	
    for line in fin:
	token = line.rstrip().split()
	if len(token)<1 : continue	
	if token[0]=="MASS": atomTypes.append(token[2])
    return atomTypes

def getMolDrudeLJ(atomTypes,atomMap,drudeMasterLJ,cGenMaster,molLJ):
    ljs = {}
    wat = {"HT_W":["-0.046000","0.224500"],"OT_W":["-0.152100","1.768200"]}
    for a in atomTypes:
	if a == "HT_W" or a=="OT_W":
		ljs[a] = wat[a]
	else:
	   if drudeMasterLJ.has_key(atomMap[a]):
		ljs[a] = drudeMasterLJ[atomMap[a]]
	   elif cGenMaster.ljs.has_key(atomMap[a]):
		ljs[a] = cGenMaster.ljs[atomMap[a]]
	   elif molLJ.has_key(a):
		ljs[a] = molLJ[a]
	   else:
		print "problem finding lj for ", a
		sys.exit()
    return ljs


if __name__=="__main__":
    #mol = os.getcwd().split("cgenff_gaamp/")[1].split("/")[0]
    bonds, angles, diheds, improps, ljs = readPara("/lcrc/project/Drude/chetan/software/cgenff/par_all36_cgenff.prm")	
    atomTypes = getAtomTypes()
    cGenMaster = Master(bonds, angles, diheds, improps, ljs)	
    #atomMap = readAtomTypeMap("/lcrc/project/Drude/chetan/OLD_GAAMP/reparameterize/"+mol+"/cgenff/cgenff_drude.txt") ### need this type to get drude specific LJ params
    
    molParams = readPara("mol-tmp.prm")
    
    """	 
    print molParams[0]
    print molParams[1]
    print molParams[2]
    print molParams[3]
    print molParams[4]
    """	
	

    bonds = getParasFromMaster(molParams[0],cGenMaster)  
    angles = getParasFromMaster(molParams[1],cGenMaster)  
    diheds = getParasFromMaster(molParams[2],cGenMaster)
    improps = getParasFromMaster(molParams[3],cGenMaster,False)
    molLJ = getParasFromMaster(molParams[4],cGenMaster)

    #drudeMasterLJ = readAlexLJ()
    #molLJ = getMolDrudeLJ(atomTypes,atomMap,drudeMasterLJ,cGenMaster,molParams[-1])

    """ 
    print bonds
    print angles
    print diheds
    print improps
    print molLJ
    """

    fout = open("mol.prm","w")
    print >> fout, "* MINI FORCE FIELD PARAMETER FILE."
    print >> fout, "*\n"
    print >> fout, "BONDS"
    
    for p in bonds:
	bs = p.split("-")
	print >> fout, bs[0],bs[1],bonds[p][0],bonds[p][1]
    os.system("grep 'HT_W   HT_W' mol-cgen-org.prm >> mol.prm")
    os.system("grep 'OT_W   HT_W' mol-cgen-org.prm >> mol.prm")
    print >> fout, "\n"

    print >> fout, "ANGLES"
    for p in angles:
	ang = p.split("-") 
	if len(angles[p])>2:
	    	print >> fout, ang[0],ang[1],ang[2],angles[p][0],angles[p][1],angles[p][2],angles[p][3] 
	else:
	    	print >> fout, ang[0],ang[1],ang[2],angles[p][0],angles[p][1] 
    os.system("grep 'HT_W   OT_W   HT_W' mol-cgen-org.prm >> mol.prm")
    print >> fout, "\n"
  
    print >> fout, "DIHEDRALS"
    for p in diheds:
	ds = p.split("-")
	if len(diheds[p]) > 1:
	   for v in diheds[p]:
		print >> fout, ds[0],ds[1],ds[2],ds[3],v[0],v[1],v[2]
	else: print >> fout, ds[0],ds[1],ds[2],ds[3],diheds[p][0][0],diheds[p][0][1],diheds[p][0][2]
    print >> fout, "\n"


    print >> fout, "IMPROPERS"
    for p in improps:
        ims = p.split("-")
        print >> fout, ims[0],ims[1],ims[2],ims[3],improps[p][0],improps[p][1],improps[p][2]
    print >> fout, "\n"

    print >> fout, "NONBONDED  E14FAC  1.000000"
    for p in molLJ:
	print >> fout, p, "0.00", molLJ[p][0], molLJ[p][1]

    print >> fout, "\n"

    print >> fout, "END"
    fout.close()


