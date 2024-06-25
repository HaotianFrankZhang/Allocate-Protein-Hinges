import numpy as np
import os
def removeRepeat(list):
    newList = []
    for item in list:
        if item not in newList:
            newList.append(item)
    return newList

def splitChar(str):
    curr = ''
    flag = -10
    newStr = []
    for i in range(len(str)):
        if str[i] == '.':
            flag = i
        if i == flag + 4:
            newStr.append(curr)
            curr = ''
        curr += str[i]
    newStr.append(curr)
    return newStr

def spliteLine(line):
    newLine = []
    for item in line:
        currItem = splitChar(item)
        if len(currItem) > 1:
            newLine += currItem
        else:
            newLine.append(item)
    return newLine


def chainPosSplit(line):
    for i in range(len(line[4])):
        if line[4][i].isnumeric():
            break
    curr = line[:4] + [line[4][:i]] + [line[4][i:]] + line[5:]
    newCurr = []
    for item in curr:
        if item != '':
            newCurr.append(item)
    return newCurr

def checkNumeric(str):
    isNumber = True
    try:
        float(str)
    except ValueError:
        isNumber = False
    return isNumber

def getSubUnitName(currLine):
    number = 0
    for i in range(len(currLine)):
        if checkNumeric(currLine[i]):
            number += 1
            if number == 2:
                break
    if len(currLine[i-1]) == 1:
        return currLine[i-1]
    else:
        return currLine[i-2]


def isflaot(str):
    if str == 'NAN':
        return False
    else:
        try:
            float(str)
            try:
                int(str)
                return False
            except ValueError:
                return True
        except ValueError:
            return False

def AtomProteinSplit(word):
    split = -1
    isLast = False
    for i in range(len(word)):
        if word[i].isnumeric():
            split = i
        if split > 0 and word[i].isalpha():
            break

    return word[: split+1], word[split + 1:]

def getCoords(line):
    start = 0
    for i in range(2, len(line)):
        if isflaot(line[i]) and 'E' not in line[i]:
            start = i
            break
    coords = line[start:start+3]

    Atom = line[2]
    Residue = line[start - 3]
    if Atom != Residue:
        out = [Atom] + [line[3] + line[start-1]] + coords
    else:
        Atom, residue = AtomProteinSplit(Atom)
        out = [Atom] + [residue + line[start-1]] + coords
    return out

def checkIfWater(line):
    isWater = False
    for item in line:
        if 'HOH' in item:
            isWater = True
    return isWater

def splitTitle(word):
    split = -1
    for i in range(len(word)):
        if word[i].isnumeric():
            split = i
            break
    return word[:split], word[split :]

def findSmallMole(words, ref):
    isFind = False
    for eachRef in ref:
        if eachRef in words:
            isFind = True
            return eachRef
    if not isFind:
        return words

def file_name(file_dir):
    scPDB = {}
    PDB_id = []
    for root, dirs, files in os.walk(file_dir):
        if root[8:] not in scPDB.keys():
            scPDB[root[8:]] = files
            PDB_id.append(root[8:])
        else:
            print ("check dataset!")
    return scPDB, PDB_id

def getChainID(str):
    # print ('currstr is', str)
    str = str.strip()
    info = str.split(';')
    newInfo = []
    for each in info:
        if each != '':
            each = each.strip()
            newInfo.append(each)

    chain = []
    for item in newInfo:
        if 'CHAIN' in item:
            curr = item.split(':')
            chain.append(curr[-1].strip())
    # print (chain)
    return chain

def filterName(str):
    allInfo = str.split(';')
    # print(allInfo)
    newStr = ''
    for item in allInfo:
        if 'MOLECULE' in item or 'SYNONYM' in item:
            start = False
            end = False
            for char in item:
                if char == ';':
                    start = False
                    end = True
                if start and not end:
                    newStr += char
                if char == ':':
                    start = True
            newStr += ' '
            newStr += '/'
    # print (newStr[:-1].strip())
    return newStr[:-1].strip()

def removeDetail(str):
    newStr = ''
    isStart = True
    for char in str:
        if char == '(':
            isStart = False
        if isStart:
            newStr += char
        if char == ')':
            isStart = True
    return newStr.strip()

def listToStr(list):
    currStr = ''
    for item in list:
        currStr += item + ' '
    return currStr

def filterSpace(list, keyword):
    newList = []
    for item in list:
        if item != '' and item not in keyword:
            newList.append(item)
    return newList

# calculate the interactions
def deBugCoord(coord):
    newCoord = []
    for item in coord:
        if isflaot(item):
            newCoord.append(item)
        else:
            split = 0
            isSplit = False
            for i in range(len(item)):
                if item[i] == '.' and not isSplit:
                    split = i + 3
                    isSplit = True
            newCoord.append(item[:split + 1])
            newCoord.append(item[split + 1:])
    return newCoord[:3]

import numpy as np
def calDistance(coord1, coord2):
    coord1 = deBugCoord(coord1)
    coord2 = deBugCoord(coord2)
    coor1 = np.array(list(map(float, coord1)))
    coor2 = np.array(list(map(float, coord2)))
    dist = np.sqrt(np.sum((coor1 - coor2) ** 2))
    return dist

def countInteractions(interaction):
    stat = {}
    for item in interaction:
        if item[1] not in stat.keys():
            stat[item[1]] = [item[2]]
        else:
            stat[item[1]].append(item[2])
    count = {}
    for key in stat.keys():
        count[key] = [len(stat[key]), np.mean(stat[key])]
    return count

def Interactions(currMolecule, currProtein, molecule, protein, Threshold = 5):
    closeAtom = {}
    # print (molecule[currMolecule])
    for atom in molecule[currMolecule]:
        dist = []
        # print ('currMolecule', currMolecule)
        # print(atom)
        for proAtom in protein[currProtein]:
            currDist = calDistance(atom[2:], proAtom[2:])
            if currDist <= Threshold:
                dist.append([proAtom[0], proAtom[1], currDist])
        # print (atom)
        # print (atom[0])
        if atom[0] not in closeAtom.keys():
            closeAtom[atom[0]] = dist
        else:
            closeAtom[atom[0]] = closeAtom[atom[0]] + dist
    closedList = []
    # print (closeAtom)
    for atom in closeAtom:
        closedList += closeAtom[atom]
    stat = countInteractions(closedList)
    return stat

def getPosition(str):
    startPos = 0
    endPos = 0
    for i in range(len(str)):
        if str[i].isnumeric() and startPos == 0:
            startPos = i
        if str[i].isalpha() and startPos != 0:
            endPos = i - 1
        else:
            endPos = i
    return str[startPos: endPos + 1]

def maxOverlap(str1, str2):
    c = []
    for i in range(len(str1)):
        for j in range(len(str1)):
            if j + i <= len(str1):
                if str1[j:j + i + 1] in str2:
                    c.append(str1[j:j + i + 1])
    maxLength = 0
    overlap = ''
    # print ('c is', c)
    for item in c:
        if len(item) > maxLength:
            overlap = item
            maxLength = len(item)
    return overlap

def getResidueNumber(str):
    number = ''
    for char in str:
        if char.isnumeric():
            number += char
    return number

def merge_dictionaries(*dicts):
  res = dict()
  for d in dicts:
    res.update(d)
  return res

import re
def isfloat(num):
    pattern = re.compile(r'(.*)\.(.*)\.(.*)')
    if pattern.match(num):
        return False
    num = num.replace("-", "")
    return num.replace(".", "").isdigit()

def isCoordinate(str):
    if isfloat(str) and len(str) >= 5:
        if str[-4] == '.':
            return True
    else:
        curr = str.split('-')
        for each in curr:
            if not isfloat(each) or len(each) < 5:
                return False
        return True

def findResidueNumber(currline):
    j = 0
    for i in range(len(currline)):
        if isCoordinate(currline[i]):
            break
        j += 1
    return currline[j-1]

def getNumber(str):
    str = str[::-1]
    newStr = ''
    isStart = False
    for i in range(len(str)):
        if str[i].isnumeric():
            newStr += str[i]
            isStart = True
        if isStart and not str[i].isnumeric():
            break
    return newStr[::-1]

def getChain(line):
    # print (line)
    if len(line[4]) == 1 and line[4].isalpha():
        return line[4]
    else:
        for i in range(len(line)):
            if isCoordinate(line[i]):
                break
        if len(line[i-2]) == 1 and line[i-2].isalpha():
            return str(line[i-2])
        else:
            # print (i)
            return str(line[i-1][0])


def updateMol(currStr, index):
    '''
    if '_' not in str:
        return str + '_1'
    else:
        Info = str.split('_')
        index = int(Info[1])
        Moleculee = Info[0]
        return Moleculee + '_' + str(index + 1)
    '''
    return currStr + '_' + str(index)

def str2float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def getAlpha(str):
    newStr = ''
    for char in str:
        if char.isalpha():
            newStr += char
    return newStr

def findDrugID(line):
    PossibleIndex = []
    for i in range(len(line)):
        if str2float(line[i]):
            PossibleIndex.append(max(0, i - 1))
            PossibleIndex.append(max(0, i - 2))

    # print (PossibleIndex)
    for index in PossibleIndex:
        if len(getAlpha(line[index])) == 1:
            break
    return getAlpha(line[index])

def findSmallMolName(previousMol, currMol):
    isNew = True
    if len(previousMol) == 0:
        return currMol, isNew
    else:
        flagMol = currMol
        if currMol[-3:] == previousMol[-3:]:
            flagMol = previousMol
            isNew = False
        return flagMol, isNew

# test = ['HETATM', '1642', 'C82B1UN', 'A1001-13.602', '14.069', '28.325', '0.50', '20.47', 'C']
def splitCoordPosition(line):
    newLine = []
    for item in line:
        if '.' in item and len(item) >= 5 and '-' in item:
            if isflaot(item):
                newLine.append(item)
            else:
                currNew = item.split('-')
                # print (currNew)
                newLine.append(currNew[0])
                newLine.append('-' + currNew[1])
        else:
            newLine.append(item)
    return newLine
# print ('test split', splitCoordPosition(test))

def findAtoms(line):
    if line[0] == 'HETATM' or line[0] == 'ATOM':
        position = 2
    else:
        position = 1
    currAtom = line[position]
    isNumber = False
    Atom = ''
    isAtom = True
    Rest = ''

    for char in currAtom:
        if char.isalpha() and isNumber:
            isAtom = False
        if isAtom:
            if char.isalpha() and not isNumber:
                Atom += char
            if char.isnumeric():
                isNumber = True
                Atom += char
        else:
            Rest += char
    return Atom, Rest

def findLigandCoord(line):
    Coordinate = []
    j = 100
    for i in range(len(line)):
        if isCoordinate(line[i]):
            Coordinate.append(line[i])
            if i < j:
                j = i
    return Coordinate, j

def reverseStr(posiblePosition):
    posiblePosition = list(posiblePosition)
    posiblePosition.reverse()
    posiblePosition = ''.join(posiblePosition)
    return posiblePosition

def findMolPosition(line, pos):
    posiblePosition = line[pos - 1]
    posiblePosition = reverseStr(posiblePosition)
    # print ('pos position', posiblePosition)

    number = ''
    rest = ''
    isPos = True
    for char in posiblePosition:
        if char.isnumeric() and isPos:
            number += char
        else:
            rest += char
            isPos = False
    return reverseStr(number), rest

def findMolecule(line, Rest, rest, pos):
    # print ('check findMolecule', line, 'Rest', Rest, 'rest', rest, 'pos', pos)
    if len(Rest) >= 2:
        return Rest
    else:
        if rest == '':
            potentialPosition = pos - 3
        else:
            potentialPosition = pos - 2
        # print ('Result', line[potentialPosition], ' ')
        curr = line[potentialPosition]
        if len(curr) >= 2:
            return curr
        else:
            return line[potentialPosition + 1]

def findClosest(molList, mol):
    newlist = [ele for ele in reversed(molList)]
    dist = 0
    if molList == [] or mol not in newlist:
        return 0
    else:
        for item in newlist:
            dist += 1
            if item in mol or mol in item:
                break
        return dist

class PDB_file():
    def __init__(self, pdb_id):
        self.PDB = pdb_id
        self.PDBInfo = self.parseFile()
        self.title = self.getTitle()
        self.seqlength, self.sequence = self.calSeqLength()
        self.totalLength = self.getTotalLength()
        self.ligand = self.getKeyWords()
        # {chain:[[atom, residue, coordinate], []]}
        self.proteinCoordinate, self.ligandCoordinate = self.getCoordinate()
        self.residuePosition = {}
        self.proteinName, self.proteinName_details, self.chain = self.getProteinName() # lists of protein name
        self.Residue_GNM, self.MissingPosition = self.matchResidues()

    # Get basic information
    def parseFile(self):
        PDBInfo = []
        with open('./PDB/' + self.PDB + '.pdb', 'r') as rf:
            for line in rf.readlines():
                line = line.strip('\n')
                line = line.strip()
                line = line.split()
                PDBInfo.append(line)
        return PDBInfo

    def getTitle(self):
        currTitle = ''
        for line in self.PDBInfo:
            if line[0] == 'TITLE':
                for item in line[1:]:
                    # if not item.isnumeric():
                    currTitle += item + ' '
        return currTitle.strip()

    # sequence and sequence length
    def calSeqLength(self):
        SeqLength = []
        Sequence = {}
        for line in self.PDBInfo:
            if line[0] == 'SEQRES':
                if line[2] not in Sequence.keys():
                    Sequence[line[2]] = []
                Sequence[line[2]] += line[4:]
                newLine = []
                for i in range(len(line)):
                    if line[i] != '':
                        newLine.append(line[i])
                SeqLength.append(newLine[2:4])
        SeqLength = removeRepeat(SeqLength)  # [[chainID, number], [chainID, number], ...]
        return SeqLength, Sequence

    def getTotalLength(self):
        totalLength = 0
        for item in self.seqlength:
            totalLength += int(item[1])
        '''
        indexLength = {}
        for line in self.PDBInfo:
            if 'ATOM' == line[0]:
                currChain = getChain(line)
                # print (line)
                # print (currChain)
                if currChain not in indexLength.keys():
                    indexLength[currChain] = 0
                # print (getNumber(findResidueNumber(line)))
                # print ('index length is', indexLength)
                # print (findResidueNumber(line))
                indexLength[currChain] = max(indexLength[currChain], int(getNumber(findResidueNumber(line))))
        '''
        return totalLength

    def indexMissing(self):
        totalLength = 0
        for key in self.indexLength.keys():
            totalLength += self.indexLength[key]
        return 1 - (len(self.Residue_GNM) / totalLength)

    # get coordinate of protein and ligand
    def getCoordinate(self):
        ParseProtein = {}
        ParseMolecule = {}
        i = 0
        isProtein = False
        Ref_Mol = []
        index = 0
        previousMol = ''
        Prev_Chain = ''
        Ref_ID = {}
        for line in self.PDBInfo:
            line = spliteLine(line)
            if 'ENDMDL' == line[0]:
                break
            if line[0] == 'ATOM' and isProtein == False:
                isProtein = True
                i += 1
                line_split = chainPosSplit(line)
                ParseProtein[getChain(line)] = [getCoords(line_split)]
            elif line[0] == 'ATOM' and isProtein:
                line_split = chainPosSplit(line)
                if getChain(line) not in ParseProtein.keys():
                    ParseProtein[getChain(line)] = []
                ParseProtein[getChain(line)].append(getCoords(line_split))
            elif 'HETATM' in line[0] and not checkIfWater(line):
                if len(line[0]) > len('HETATM'):
                    title, id = splitTitle(line[0])
                    line = [title] + [id] + line[1:]
                if len(line) == 12:
                    Atom = line[2]
                    if len(line[3]) == 1:
                        currMol = line[4]
                    else:
                        currMol = line[3]
                    mol_ID = line[5]
                    Coords = line[6:9]
                    # print ('Coords', Coords)
                    # print ('molecule', currMol)
                else:
                    line = splitCoordPosition(line)
                    # print('currLine', line)
                    # print('coord', getCoords(line))
                    Atom, Rest = findAtoms(line)
                    # print('atom', Atom, Rest)
                    Coords, pos = findLigandCoord(line)
                    # print ('coordinate', Coords, pos)
                    mol_ID, rest = findMolPosition(line, pos)
                    # print ('position', mol_ID, rest)
                    currMol = findMolecule(line, Rest, rest, pos)

                # print ('currMol is', currMol)
                currChain = ''
                chainPlace = 0
                for item in line[:-1]:
                    if len(item) == 1 and item.isalpha() and chainPlace > 2:
                        currChain = item
                    chainPlace += 1

                if currChain == '':
                    if len(line[4]) == 1:
                        currChain = line[4]
                    elif len(rest) == 1:
                        currChain = rest
                    else:
                        currChain = line[pos - 1][0]
                # print ('currChain is', currChain)

                # print ('RefMol is', Ref_Mol)
                closestPos = findClosest(Ref_Mol, currMol)
                # print ('distance', closestPos)

                # different ligands are only based on their ID
                if mol_ID not in Ref_ID.keys():
                    Ref_ID[mol_ID] = [currMol, currChain]
                    index += 1
                elif closestPos >= 10 or currChain not in Ref_ID[mol_ID]:
                    index += 1
                    Ref_ID[mol_ID].append(currChain)

                # print (Ref_ID)
                # print('currIndex', index)

                previousMol = currMol
                # Prev_Chain = currChain
                Ref_Mol.append(currMol)

                # print ('getCoords', getCoords(line))
                currInfo = [Atom, currMol + mol_ID] + Coords
                # print ('currInfo is', currInfo)

                if Ref_ID[mol_ID][0] + '_' + str(index) not in ParseMolecule.keys():
                    ParseMolecule[Ref_ID[mol_ID][0] + '_' + str(index)] = [currInfo]
                else:
                    ParseMolecule[Ref_ID[mol_ID][0] + '_' + str(index)].append(currInfo)

        '''if line[0] == 'ATOM' and isProtein == False:
            isProtein = True

            line = splitCoordPosition(line)
            Atom, Rest = findAtoms(line)
            # print('atom', Atom, Rest)
            Coords, pos = findLigandCoord(line)
            # print ('coordinate', Coords, pos)
            mol_ID, rest = findMolPosition(line, pos)
            # print ('position', mol_ID, rest)
            currMol = findMolecule(line, Rest, rest, pos)
            currInfo = [Atom, currMol + mol_ID] + Coords
            ParseProtein[getChain(line)] = [currInfo]
        elif line[0] == 'ATOM' and isProtein:
            # print ('currLine', line)
            # line_split = chainPosSplit(line)
            # if getChain(line) not in ParseProtein.keys():
            #     ParseProtein[getChain(line)] = []
            # ParseProtein[getChain(line)].append(getCoords(line_split))
            # print('coord', getCoords(line_split))
            # line = splitCoordPosition(line)
            if len(line) < 12:
                Atom, Rest = findAtoms(line)
                # print('atom', Atom, Rest)
                Coords, pos = findLigandCoord(line)
                # print ('coordinate', Coords, pos)
                mol_ID, rest = findMolPosition(line, pos)
                # print ('position', mol_ID, rest)
                currMol = findMolecule(line, Rest, rest, pos)
                currInfo = [Atom, currMol + mol_ID] + Coords
                # print(line)
                # print('coord', currInfo)
                # print (getChain(line))
                if getChain(line) not in ParseProtein.keys():
                    ParseProtein[getChain(line)] = []
                ParseProtein[getChain(line)].append(currInfo)
            else:
                if getChain(line) not in ParseProtein.keys():
                    ParseProtein[getChain(line)] = []
                ParseProtein[getChain(line)].append(getCoords(line))
                # print(line, getCoords(line)) 
                
                
                '''
        '''
        if line[0] == 'ATOM' and isProtein == False:
            isProtein = True
            i += 1
            line_split = chainPosSplit(line)
            ParseProtein[getChain(line)] = [getCoords(line_split)]
        elif line[0] == 'ATOM' and isProtein:
            line_split = chainPosSplit(line)
            if getChain(line) not in ParseProtein.keys():
                ParseProtein[getChain(line)] = []
            ParseProtein[getChain(line)].append(getCoords(line_split))
        '''

        # print ('protein', ParseMolecule)
        return ParseProtein, ParseMolecule

    def getProteinName(self):
        currProtein = []
        curr = ''
        for line in self.PDBInfo:
            line = listToStr(line)
            if line.startswith('COMPND'):
                if 'MOL_ID' in line:
                    currProtein.append(curr)
                    curr = ''
                else:
                    curr += line + ' '
        currProtein.append(curr)

        filter = []
        chain = []
        for item in currProtein:
            if item != '':
                filter.append(filterName(item))
                chain.append(getChainID(item))
        MajorInfo = []
        for item in filter:
            MajorInfo.append(removeDetail(item))
        MajorInfo = list(set(MajorInfo))
        return MajorInfo, filter, chain

    def getKeyWords(self):
        currMatch = {}
        keyWord = ['HETSYN', 'HETNAM']
        shouldMerge = False
        for line in self.PDBInfo:
            if 'REVDAT' not in line:
                if keyWord[1] in line:
                    curr = filterSpace(line, keyWord)
                    if curr[0].isnumeric() and len(curr[0]) < 3:
                        curr = curr[1:]
                        shouldMerge = True
                    if curr[0] not in currMatch.keys():
                        currMatch[curr[0]] = []
                    currMatch[curr[0]] += curr[1:]
                    if shouldMerge:
                        temp =  currMatch[curr[0]][-2] + currMatch[curr[0]][-1]
                        currMatch[curr[0]] = currMatch[curr[0]][:-3]
                        currMatch[curr[0]].append(temp)
                        shouldMerge = False
                elif keyWord[0] in line:
                    curr = filterSpace(line, keyWord)
                    if curr[0].isnumeric() and len(curr[0]) < 3:
                        curr = curr[1:]
                        shouldMerge = True
                    if curr[0] not in currMatch.keys():
                        currMatch[curr[0]] = []
                    currMatch[curr[0]] += ['/'] + curr[1:]
                    if shouldMerge:
                        temp =  currMatch[curr[0]][-2] + currMatch[curr[0]][-1]
                        currMatch[curr[0]] = currMatch[curr[0]][:-3]
                        currMatch[curr[0]].append(temp)
                        shouldMerge = False
        return currMatch

    def calBinding(self, threshold=5):
        proteinList = list(self.proteinCoordinate.keys())
        moleculeList = list(self.ligandCoordinate.keys())
        # print (self.proteinCoordinate)
        # print (self.ligandCoordinate['2FD'])

        parsedPDB = {}
        for currMolecule in moleculeList:
            if currMolecule != 'HOH':
                for currProtein in proteinList:
                    currStat = Interactions(currMolecule, currProtein, self.ligandCoordinate, self.proteinCoordinate, threshold)
                    # print('currMol', currMolecule)
                    # print ('currProtein', currProtein)
                    # print ('currStat is', currStat)
                    if currMolecule not in parsedPDB.keys():
                        parsedPDB[currMolecule] = {}

                    parsedPDB[currMolecule][currProtein] = currStat
        # print (parsedPDB['NAI'])
        return parsedPDB

    def matchResidues(self):
        Residue_GNM = {}
        MissingPosition = []
        GNM_id = 0

        for chain in self.proteinCoordinate.keys():
            for eachAtom in self.proteinCoordinate[chain]:
                if eachAtom[0] == 'CA':
                    Residue_GNM[GNM_id] = [chain, eachAtom[1]]
                    GNM_id += 1

        # print (Residue_GNM)
        previousPos = ''
        previousRes = ''
        previousChain = ''
        for id in Residue_GNM.keys():
            currChain = Residue_GNM[id][0]
            if previousPos == '' or previousChain != currChain:
                previousPos = getPosition(Residue_GNM[id][1])
                previousRes = Residue_GNM[id][1]
            currPos = getPosition(Residue_GNM[id][1])

            # print (currPos, previousPos)
            if int(currPos) - int(previousPos) > 1:
                MissingPosition.append([Residue_GNM[id][0], previousRes, Residue_GNM[id][1], int(currPos) - int(previousPos) - 1])
            previousPos = currPos
            previousRes = Residue_GNM[id][1]
            previousChain = currChain
        return Residue_GNM, MissingPosition

    def calMissSeq(self):
        missingPorportion = {}
        totalMissing = 0
        for item in self.MissingPosition:
            if item[0] not in missingPorportion.keys():
                missingPorportion[item[0]] = 0
            missingPorportion[item[0]] += item[-1]
        for chain in missingPorportion.keys():
            for item in self.seqlength:
                if item[0] == chain:
                    length = item[-1]
                    break
            for eachChain in missingPorportion.keys():
                totalMissing += float(missingPorportion[eachChain])
            missingPorportion[chain] = (missingPorportion[chain]) / (float(length))
        totalMissing = totalMissing / self.totalLength
        return merge_dictionaries(missingPorportion, {'TOTAL': totalMissing})

def readPDBList(fileName):
    list = []
    with open(fileName, 'r') as rf:
        for line in rf.readlines():
            line = line.replace('\n', '')
            line = line.strip()
            list.append(line)
    return list

def readList(fileName):
    list = []
    with open(fileName, 'r') as rf:
        for line in rf.readlines():
            line = line.replace('\n', '')
            line = line.strip()
            line = line.split('\t')
            list.append(line)
    return list

def write2DList(fileName, data):
    with open(fileName, 'w') as wf:
        for line in data:
            for item in line:
                wf.write(item + '\t')
            wf.write('\n')

def writeBindingSites(allBinding, fileName='./Final_Binding5angstroms/cleanedAllBindingSites.txt'):
    with open(fileName, 'w') as wf:
        for pdb in allBinding.keys():
            for smal in allBinding[pdb].keys():
                for chain in allBinding[pdb][smal].keys():
                    wf.write(pdb + '\t' + smal + '\t')
                    wf.write("'" + chain + "'" + '\t')
                    for resi in allBinding[pdb][smal][chain].keys():
                        wf.write(resi + '\t')
                    wf.write('\n')

def writeBindingStrength(allBinding, fileName = './Final_Binding5angstroms/cleanBindingStrength.txt'):
    with open(fileName, 'w') as wf:
        for pdb in allBinding.keys():
            for smal in allBinding[pdb].keys():
                for chain in allBinding[pdb][smal].keys():
                    wf.write(pdb + '\t' + smal + '\t')
                    wf.write("'" + chain + "'" + '\t')
                    for resi in allBinding[pdb][smal][chain].keys():
                        strength = allBinding[pdb][smal][chain][resi]
                        wf.write(resi + '\t' + str(strength[0]) + '\t' + str(strength[1]) + '\t')
                    wf.write('\n')

def writeMissingResidue(allMissingProportion, filename='MissingSeqProportion.txt'):
    with open(filename, 'w') as wf:
        for key in allMissingProportion:
            wf.write(key + '\t')
            for eachChain in allMissingProportion[key]:
                wf.write(eachChain + '\t')
                wf.write(str(allMissingProportion[key][eachChain]) + ';\t')
            wf.write('\n')

def writeResidue_GNM(allResidue_GNM, filename='./Final_Binding4.5angstroms/Residue45A_CA.txt'):
    with open(filename, 'w') as wf:
        for pdb in allResidue_GNM.keys():
            wf.write(pdb + '\t')
            for title in allResidue_GNM[pdb].keys():
                wf.write(str(title) + '\t')
                for each in allResidue_GNM[pdb][title]:
                    wf.write(str(each + '\t'))
                wf.write(';')
            wf.write('\n')

# not used in parsing data
def parseTitle(currPDB):
    currAllNames = []
    for item in currPDB.proteinName:
        currNames = item.split(' ')
        updateNames = ''
        print('currname', currNames)
        CompoundPosition = []
        for i in range(len(currNames)):
            if currNames[i] == 'COMPND':
                CompoundPosition.append(i - 1)
                CompoundPosition.append(i)
                CompoundPosition.append(i + 1)
        for i in range(len(currNames)):
            if currNames[i] == 'COMPND':
                CompoundPosition.append(i)
    print (currAllNames)

def writediction_twoDArray(data, fileName):
    with open('./Final_Binding5angstroms/' + fileName, 'w') as wf:
        for pdb in data.keys():
            wf.write(pdb + '\t')
            for item in data[pdb]:
                for each in item:
                    wf.write(each + '\t')
                wf.write(';')
            wf.write('\n')

def writedictionary(data, fileName):
    with open('./Final_Binding5angstroms/' + fileName, 'w') as wf:
        for pdb in data.keys():
            wf.write(pdb + '\t')
            wf.write(str(data[pdb]))
            wf.write('\n')

def BindingSubUnitLength(Binding, Seqlength):
    LigandBindingLength = {}
    Length = {}
    for item in Seqlength:
        Length[item[0]] = int(item[1])
    for ligand in Binding.keys():
        for chain in Binding[ligand].keys():
            if ligand not in LigandBindingLength.keys():
                LigandBindingLength[ligand] = 0
            # print (LigandBindingLength)
            # print (Length, chain)
            if Binding[ligand][chain] != {}:
                LigandBindingLength[ligand] += Length[chain]
    return LigandBindingLength

def main():
    # finalList = readList('finalList_withoutAA.txt')
    finalList = readList('ids.txt')

    # print(finalList)
    allList = []
    for item in finalList:
        if item not in allList:
            allList.append(item[0])
    # allList = list(set(allList))

    # print (allList)
    # print (len(allList))

    # currPDB = PDB_file('4pwd')
    # print (currPDB.Residue_GNM)
    # print (currPDB.chain)
    # print (currPDB.calBinding().keys())
    # print (currPDB.calBinding(threshold=4.5))

    '''
    totalMissing = {}
    for i in range(len(allList)):
        print (str(i+1) + ' ' + allList[i] + ' Start!')
        curr = allList[i]
        currPDB = PDB_file(curr)
        if currPDB.totalLength <= 1000:
            totalMissing[allList[i]] = currPDB.indexMissing()
    print (totalMissing)
    writedictionary(totalMissing, 'totalMissingOnlyBasedOnIndex.txt')
    '''


    allBinding = {}
    allMissingProportion = {}
    allResidue_GNM = {}
    allProteinNames = {}
    allHETNAM = {}
    allProteinLength = {}

    for i in range(len(allList)):
        print (str(i+1) + ' ' + allList[i] + ' Start!')
        curr = allList[i]
        # curr = '1m0w'
        currPDB = PDB_file(curr)
        # print ('Binding', currPDB.calBinding())
        # print ('length', currPDB.seqlength)
        # print ('Results', max(BindingSubUnitLength(currPDB.calBinding(), currPDB.seqlength).values()))
        if allList[i] not in ['1rvv']: # max(BindingSubUnitLength(currPDB.calBinding(), currPDB.seqlength).values()) <= 1000:
            allBinding[curr] = currPDB.calBinding(threshold = 5)
            allMissingProportion[curr] = currPDB.calMissSeq()
            allResidue_GNM[curr] = currPDB.Residue_GNM
            # allProteinNames[curr] = [currPDB.proteinName]
            # allProteinLength[curr] = currPDB.seqlength
            # allHETNAM[curr] = []

            
            # for key in currPDB.ligand.keys():
            #     currHETNAM = []
            #     currHETNAM.append(key)
            #     for each in currPDB.ligand[key]:
            #         currHETNAM.append(each)
            #     allHETNAM[curr].append(currHETNAM)
            
    import os
    isExist = os.path.exists('./Final_Binding5angstroms')
    if not isExist:
       # Create a new directory because it does not exist
       os.makedirs('./Final_Binding5angstroms')
       print("The new directory is created!")

    # writediction_twoDArray(allProteinLength, 'proteinLength5A_sup.txt')
    # writediction_twoDArray(allProteinNames, 'proteinNames5A_sup.txt')
    # writediction_twoDArray(allHETNAM, 'ligandNames5A_sup.txt')
    writeBindingSites(allBinding, fileName='./Final_Binding5angstroms/BindingSites5A_ID_basedLigand.txt')
    writeMissingResidue(allMissingProportion, filename='./Final_Binding5angstroms/MissingSeqProportion.txt')
    writeResidue_GNM(allResidue_GNM, './Final_Binding5angstroms/Residue5A_CA.txt')
    writeBindingStrength(allBinding, fileName='./Final_Binding5angstroms/BindingStrength5A_ID_basedLigand.txt')




    '''
    checkList = []
    noStructurelist = []
    PDBlist = readPDBList('pdb_list.txt')
    allData = []

    for item in PDBlist:
        try:
            print (item)
            curr = PDB_file(item)
            currData = [item]
            currData.append(curr.title)
            # print (curr.ligand)
            for key in curr.ligand:
                currItem = ''
                for mol in curr.ligand[key]:
                    if mol != 'ION':
                        currItem += mol
                currData.append(currItem)
            allData.append(currData)
        except FileNotFoundError:
            noStructurelist.append([item])
        except:
            checkList.append([item])

    write2DList('PDB_KeyWords.txt', allData)
    write2DList('unkownStructure.txt', checkList)
    write2DList('noStructureList.txt', noStructurelist)
    '''

if "__name__" == main():
    main()