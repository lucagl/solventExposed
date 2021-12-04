#            **MIT License**
#
#       Copyright (c) 2021 Luca Gagliardi 
#   Affiliation Istituto Italiano di tecnologia

#PROBLEM: hot to filter out atoms which are exposed towards internal cavities?
# the option fill cavity seems to un-affect exposed list..

import re
import numpy as np
import sys
import subprocess





##### NS UTILITIES


def setup_NSInput(confFile,radius,nThreads=1,grid_scale = 2,grid_selfInt = 2,maxProbes_selfInt=100, gridPerfil = 90, isSkin=False, accTriang=False, pInput = False):
    #Populate this function with other set_ups.. I could also inmagine an actual hardCoding of NS input file
    s=0
    # print(confFile)
    NS_input= open(confFile,'r')
    content = NS_input.readlines()
    s=0
    for line in content:
        # print(line)
        match = re.match("(Grid_scale\s*=\s*)(\d+)",line)
        match1 = re.match("(Grid_perfil\s*=\s*)(\d+)",line)
        match2 = re.match("(Self_Intersections_Grid_Coefficient\s*=\s*)(\d+)",line)
        match3 = re.match("(Max_Probes_Self_Intersections\s*=\s*)(\d+)",line)
        match4 = re.match("(Accurate_Triangulation\s*=\s*)([a-zA-Z]+)",line)
        match5 = re.match("(Vertex_Atom_Info\s*=\s*)([a-zA-Z]+)",line)
        match6 = re.match("(Compute_Vertex_Normals\s*=\s*)([a-zA-Z]+)",line)
        match7 = re.match("(Save_Mesh_MSMS_Format\s*=\s*)([a-zA-Z]+)",line)
        match8 = re.match("(XYZR_FileName\s*=\s*)(.+)",line)
        match9=re.match("(Surface\s*=\s*)(.+)",line)
        match10=re.match("(Cavity_Detection_Filling\s*=\s*)([a-zA-Z]+)",line)
        matchRadius = re.match("(Probe_Radius\s*=\s*)(\d*\.?\d+)",line)
        matchThread = re.match("(Number_thread\s*=\s*)(\d+)",line)
        if matchRadius:
            # print(line)
            newline = matchRadius.group(1)+re.sub('\d*\.?\d+',str(radius),matchRadius.group(2))+'\n'
            content[s] = newline
        elif matchThread:
            # print(line)
            newline = matchThread.group(1)+re.sub('\d+',str(nThreads),matchThread.group(2))+'\n'
            content[s] = newline
        elif(match):
            newline = match.group(1)+re.sub('\d+',str(grid_scale),match.group(2))+'\n'
            content[s] = newline
        elif(match1):
            newline = match1.group(1)+re.sub('\d+',str(gridPerfil),match1.group(2))+'\n'
            content[s] = newline
        elif(match2):
            newline = match2.group(1)+re.sub('\d+',str(grid_selfInt),match2.group(2))+'\n'
            content[s] = newline
        elif(match3):
            newline = match3.group(1)+re.sub('\d+',str(maxProbes_selfInt),match3.group(2))+'\n'
            content[s] = newline
        elif(match8):
            if(pInput):
                newline = match8.group(1)+re.sub('.+',"NanoShaper_Pocket_input.xyzr",match8.group(2))+'\n'
            else:
                newline = match8.group(1)+re.sub('.+',"NanoShaper_input.xyzr",match8.group(2))+'\n' 
            content[s] = newline
        elif(match4):
            if(accTriang):
                newline = match4.group(1)+re.sub('[a-zA-Z]+',"true",match4.group(2))+'\n'
            else:
                newline = match4.group(1)+re.sub('[a-zA-Z]+',"false",match4.group(2))+'\n'
            content[s] = newline
        elif(match5):
            newline = match5.group(1)+re.sub('[a-zA-Z]+',"false",match5.group(2))+'\n'
            content[s] = newline
        elif(match6):
            newline = match6.group(1)+re.sub('[a-zA-Z]+',"false",match6.group(2))+'\n'
            content[s] = newline
        elif(match7):
            newline = match7.group(1)+re.sub('[a-zA-Z]+',"false",match7.group(2))+'\n'
            content[s] = newline
        elif(match9):
            if(isSkin):
                newline = match9.group(1)+re.sub('.+',"skin",match9.group(2))+'\n'
            else:
                newline = match9.group(1)+re.sub('.+',"ses",match9.group(2))+'\n'
            content[s] = newline
        elif(match10):
            newline = match10.group(1)+re.sub('[a-zA-Z]+',"true",match10.group(2))+'\n'
            content[s] = newline
        s+=1
    NS_input.close()
    NS_input= open(confFile,'w')
    NS_input.writelines(content)
    NS_input.close()
    return









########
####

####### READING UTILITIES ###########




def fetchRes(filename):
    """
    CAREFUL: if chain not available, dummy chain A assumed
    """

    # global resMap # List of dictionaryìies mapping: atom label (line in file) <--> all infos 
    resMap = []
    content ={}
    # old_number = None


    comment =['#', 'CRYST[0-9]?']
    remark = ['REMARK']
    termination = ['TER', 'END']
    skip = comment+remark
    skip = '(?:% s)' % '|'.join(skip)


    # print("-- Loading PQR file --")
    try:
        inFile = open(filename,'r')
    except Exception:
            raise NameError("Cannot find PQR file")
    for line in inFile: 
        if(re.match(skip,line)): 
            pass 
        else:
            linegNOChain=re.match("(ATOM)\s*(\d+)\s*(\S+)\s+([A-Z]+)\s+(\-?\d+[A-Z]?)\s+(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s+(\-?\d*\.?\d+)\s+(\d*\.?\d+)",line)
            linegChain = re.match("(ATOM)\s*(\d+)\s*(\S+)\s+([A-Z]+)\s+([\w0-9]+)\s*(\-?\d+[A-Z]?)\s+(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s+(\-?\d*\.?\d+)\s+(\d*\.?\d+)",line)
            break

    if(linegChain):
        isChainID=1                                                        #resID
        matchPattern = "(ATOM)\s*(\d+)\s*(\S+)\s+([A-Z]+)\s+([\w0-9]+)\s*(\-?\d+[A-Z]?)\s+(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s+(\-?\d*\.?\d+)\s+(\d*\.?\d+)"
    elif(linegNOChain):
        isChainID =0                                        # resID
        matchPattern = "(ATOM)\s*(\d+)\s*(\S+)\s+([A-Z]+)\s+(\-?\d+[A-Z]?)\s+(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s+(\-?\d*\.?\d+)\s+(\d*\.?\d+)"
    else:
        raise NameError("Incorrect pqr file formatting")

    if(isChainID):
        resInd = 5
        chargeInd = 9
        rInd = 10
    else:
        resInd = 4
        chargeInd = 8
        rInd = 9
    nameInd = 3
    atomInd = 2
    coordInd = resInd +1
    chainInd = 4
    
    inFile.seek(0)
    comments=[]
    for line in inFile: 
        if(re.match(skip,line)): 
            comments.append(line)
            pass 
        else:   
            match=re.match(matchPattern,line)
            if match:
                lineg = match.groups()
                if(isChainID):
                    content = {'resName':lineg[nameInd],'resNum':lineg[resInd],'atomNumber':int(lineg[1]),'resAtom': lineg[atomInd],'resChain':lineg[chainInd],
                    'charge':float(lineg[chargeInd]),'coord':list(map(float, lineg[coordInd:coordInd+3])),'radius':float(lineg[rInd])}
                else:
                    content = {'resName':lineg[nameInd],'resNum':lineg[resInd],'atomNumber':int(lineg[1]),'resAtom': lineg[atomInd],'resChain':'A',
                    'charge':float(lineg[chargeInd]),'coord':list(map(float, lineg[coordInd:coordInd+3])),'radius':float(lineg[rInd])}
                resMap.append(content)
            # print(content)
    
    
    return resMap,comments

class exposedFetcher(object):
    def __init__(self) -> None:
        self.conf = 'surfaceConfiguration.prm'
        self.workingDir = 'temp/'
        self.indices=set()
        self.resMap =[]
        self.name =''
        self.comments = ''
        self.structureInit=False
        self.probeRadius =0 
    def setExposedAtoms(self,structure_name,rp=1.4,nThreads=1):
        self.name = structure_name
        self.probeRadius = rp
        print('settung up configuration file NS')
        setup_NSInput(self.workingDir+self.conf,radius=rp,nThreads=nThreads,accTriang=False)

        print('setting up structure file %s for NS'%structure_name)
        resMap,self.comments = fetchRes(structure_name+'.pqr') #also produces NS input file by default
        protein_atoms = np.empty((0,4))
        for i in resMap:
            c = np.append(np.asarray(i['coord']),i['radius'])
            protein_atoms = np.vstack([protein_atoms,c])
        np.savetxt(self.workingDir+"NanoShaper_input.xyzr",protein_atoms,delimiter="\t",fmt='%.4f')
        try:
            #call NS
            out = subprocess.check_output('./NanoShaper',cwd=self.workingDir)
        except subprocess.CalledProcessError as grepexc:                                                                                                   
            print ("error code", grepexc.returncode, grepexc.output)

        #CHECK EXPOSED
        exposedFile = open(self.workingDir+'exposedIndices.txt','r')
        for line in exposedFile.readlines():
            self.indices.add(int(line))
        print('ALL atoms:', len(resMap))
        print('ALL residues:',len(set([(d['resNum'],d['resName'],d['resChain']) for d in resMap])))
        self.map=set([(d['resNum'],d['resName'],d['resChain']) for d in list(filter(lambda x: x['atomNumber'] in self.indices, resMap))])
        print('EXPOSED residues:', len(self.map))

        self.resMap = resMap
        self.structureInit = True
    def checkExposed(self,atoms):
        '''
        Return list of exposed residues from a input residue map. 
        Useful to determine exposed subsets.
        '''
        if(not self.structureInit):
            raise FileNotFoundError('The structure was not passed for analysis')
        residuesFiltered = set(filter(lambda x : (x['resNum'],x['resName'],x['resChain']) in self.map, atoms))#[d for d in atoms if d['atomNumber'] in self.indices]
      
        return residuesFiltered


    def getExposed(self):
        '''
        Return list of exposed residues of the whole protein
        '''
        if(not self.structureInit):
            raise FileNotFoundError('The structure was not passed for analysis')
        return self.map

    def printExposedPQR(self):
        if(not self.structureInit):
            raise FileNotFoundError('The structure was not passed for analysis')
        outFile = open(self.name+'_exposed.pqr','w')
        #COPY ORIGINAL COMMENTS IN PQR..
        outFile.writelines(self.comments)
        outFile.write('REMARK Subset of exposed atoms at probe radius %.2f\n'%self.probeRadius)
        for ind,r in enumerate(self.resMap):
            if ((r['resNum'],r['resName'],r['resChain']) in self.map):
                outFile.write("{:<6s}{:>5d} {:<5s}{:>3s} {:1s}{:>5s}   {:>8.3f} {:>8.3f} {:>8.3f} {:>8.4f} {:>8.4f}\n".format('ATOM',r['atomNumber'],
                r['resAtom'],r['resName'],r['resChain'],r['resNum'],r['coord'][0],r['coord'][1],r['coord'][2],r['charge'],r['radius']))
            else:
                pass
        outFile.write('TER\n')
        outFile.write('END')
        outFile.close()
def main():
    import getopt
    '''
    Pass inline or prompted, a PQR file.
    Returns the list of exposed residues.
    OPTIONS: - Water radius
             - save to file only resName(maybe with some sorting in the number.. )
             USE SORTED function to return a copy and avoid in place changes..
             - save to file whole PQR of exposed atoms
    '''
    SES_RADIUS = 1.4
    saveFile = False
    savePQR = False
    name=''
    argv = sys.argv[1:]
    nThreads=1
    stdoutPrint = True
    try:
        opts, args = getopt.getopt(argv,"h",["radius=","name=","saveRES","savePQR","nThreads=","help"])
    except getopt.GetoptError:
        print ('<<ERROR>> Uncorrect formatting of options or unavaible option')
        sys.exit(2)

    # print('options:',opts)
    # print('arguments:',args)

    for opt, arg in opts:
        if opt in ["-h","--help"]:
            print('Standard usage:\npython3 expRes <options> <structure name>\n The structure must be in pqr format.\nOptions:')
            print('--name=<structure_name>: overwrites passed structure name' )
            print('--radius=<float>: Probe radius. Default 1.4 (water)')
            print('--saveRES: save to txt file exposed residue number, names and chain')
            print('--savePQR: save to pqr file exposed atoms of the original structure')
            print('or..\nIt can be used as an exteranl import with more functionalities available such as getting whatever subset of exposed residues or atoms via the function "checkExposed()"')
            input('\n')
            sys.exit()
        elif opt == '--radius':
            SES_RADIUS = float(arg)
            print("Overwritten default probe radius (1.4) --> %.2f"%SES_RADIUS)
        elif opt == '--nThreads':
            nThreads = int(arg)
            print("Overwritten default number of threads (1) --> %d"%nThreads)
        elif opt == '--save':
            saveFile = True
            stdoutPrint = False
        elif opt == '--savePQR':
            savePQR = True
            stdoutPrint = False
        elif opt=='--name':
            name=str(arg)
            match = re.match('([\w]*)',name) #accepts both <pqrname> and <pqrname.pqr>
            name=match.group(1)
            # print(name)
    if((len(args)>0)and (name=='')):
        name=str(args[0])
        match = re.match('([\w]*)',name) #accepts both <pqrname> and <pqrname.pqr>
        name=match.group(1)
    print("Structure to process:",name+'.pqr')
    print("Probe radius= %.2f"%SES_RADIUS)
    print("Number of threads= %d"%nThreads)
    input('continue?')
    exp = exposedFetcher()
    exp.setExposedAtoms(name,rp=SES_RADIUS,nThreads=nThreads)
    res = exp.getExposed()
    resList = sorted(res,key=lambda t:int(t[0])) #sort according to res number
    if(saveFile):
        outFile = open(name+'_resList.txt','w')
        for r in resList:
            resid=r[0]
            rname = r[1]
            rChain = r[2]
            
            outFile.write(resid+"\t"+rname+' '+rChain+'\n')
        outFile.close()
        print('Exposed residues saved in',name+'_exposed.pqr')
    if(savePQR):
        exp.printExposedPQR()
        print('Exposed atoms saved in',name+'_resList.txt')
    if(stdoutPrint):
        for r in resList:
            print(r)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("\nUser exit")
        sys.exit()