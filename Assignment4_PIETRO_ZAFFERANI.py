''' Markov Chains as a model to recognize CpG islands - by Pietro Zafferani '''

import numpy
import random

'''This function opens a fasta file and converts the sequence contained in a Python string. '''


def convertString(file) -> str:
    File = open(file)
    String = ''
    for line in File:
        if line[0] != '>':
            String += line[:-1]
    File.close()
    return String


'''This function opens a tsv file that contains the coordinates of each CpGI. It returns a list of tuples,
    each tuple is formed by the offset, the end and the length of the CpGI.'''


def get_coordinates(file) -> list:
    L = []
    File = open(file)
    for line in File:
        row = line.split(sep='\t')
        L.append((int(row[1]), int(row[2]), int(row[3])))
    return L


'''This function takes as parameter the string of the genome considered and the list of coordinates. At each round it 
    slices the genome and it yields a region corresponding to a CpGI. '''


def get_CpGIs(genome: str, coordinates: list) -> str:
    for trio in coordinates:
        # sclicing by CPGI start and end
        yield genome[trio[0] - 1:trio[1] - 1].upper()


'''This function creates an empty matrix in which the labels are all 4 DNA nucleotide. The matrix is used to compute 
    the dimers occurrences' probability both for the inside-model and the outside model. '''


def DNA_matrix(dna='ACGT') -> list:
    M = [['/']]
    for letter in dna:
        M[0].append(letter)
        M.append([letter])
        for zero in range(1, len(dna) + 1):
            M[-1].append(0)
    return M


'''This function prints the probabilities tables in a more readable way.'''


def showMatrix(M: list) -> print:
    for row in M:
        print(row)


'''This function creates a dictionary in which the keys are the nucleotide characters and their associated values
    correspond to the indexes in the probabilities tables. In this way we can retrieve the position of a specific 
    dimer in the table only by looking at its nucleotides. '''


def accessMatrix(dna='ACGT') -> dict:
    d = {}
    for i in range(len(dna)):
        d[dna[i]] = i + 1
    return d


'''This function returns the full probability table for the inside-model.'''


def InsideModel(genome: str, coordinates: list, dna='ACGT') -> list:
    # create matrix that stores the occurrences of each specific dimer
    SpecificM = DNA_matrix()
    # create matrix that stores the occurrences of any (A,C,G,T | X) dimer
    AspecificM = DNA_matrix()
    # create the registry to access the table directly
    registry = accessMatrix()

    # iterate over every CpGI
    for CpG in get_CpGIs(genome, coordinates):
        # iterate over every nucleotide in the CpGI
        for i in range(1, len(CpG)):

            # retrieve the coordinates of the dimer in the table
            column_ind = registry[CpG[i]]
            row_ind = registry[CpG[i - 1]]

            # store occurrence in the specific table
            SpecificM[row_ind][column_ind] += 1

            # store the occurrence in the aspecific table
            for element in range(1, len(dna) + 1):
                AspecificM[row_ind][element] += 1

    # create the final table that is still empty
    INSIDEmodel = DNA_matrix()

    # divide every cell of the specific table by each cell of the aspecific table
    for r in range(1, len(SpecificM[0])):
        for c in range(1, len(SpecificM)):
            INSIDEmodel[r][c] = round(SpecificM[r][c] / AspecificM[r][c], 2)

    # we obtained the probability table containing the percentage of occurrence of every dimer
    return INSIDEmodel


'''This function takes as parameter the string of the genome considered and the list of coordinates. At each round it 
    slices in a random point the genome and it yields a region with the same length of a CpGI present in the coordinates
    file. '''


def getNo_CpGI(genome: str, coordinates: list) -> str:
    i = 0
    Max = len(genome)
    while i < len(coordinates):
        # compute random offset
        ran_index = random.randint(0, Max)
        region = genome[ran_index:ran_index + coordinates[i][2]].upper()
        # discard the regions containing unknown sequences
        if 'N' not in region:
            yield region
            i += 1


'''This function returns the full probability table for the outside-model.'''


def OutsideModel(genome: str, coordinates: list, dna='ACGT') -> list:
    # create matrix that stores the occurrences of each specific dimer
    SpecificM = DNA_matrix()
    # create matrix that stores the occurrences of any (A,C,G,T | X) dimer
    AspecificM = DNA_matrix()
    # create the registry to access the table directly
    registry = accessMatrix()

    # iterate over every region drawn randomly
    for region in getNo_CpGI(genome, coordinates):
        # iterate over each dimer of the region
        for i in range(1, len(region)):

            # retrieve the coordinates of the dimer in the table
            column_ind = registry[region[i]]
            row_ind = registry[region[i - 1]]

            # store the occurrence of a given dimer in the specific table
            SpecificM[row_ind][column_ind] += 1

            # store the occurrence of a set of dimers in the aspecific table
            for element in range(1, len(dna) + 1):
                AspecificM[row_ind][element] += 1

    # create the final table that is still empty
    OUTSIDEmodel = DNA_matrix()

    # divide every cell of the specific table by each cell of the aspecific table
    for r in range(1, len(SpecificM[0])):
        for c in range(1, len(SpecificM)):
            OUTSIDEmodel[r][c] = round(SpecificM[r][c] / AspecificM[r][c], 2)

    # we obtained the probability table containing the percentage of occurrence of every dimer
    return OUTSIDEmodel


'''This function computes the logarithmic score of a query sequence given the 2 probabilistic models.'''


def queryScore(inModel: list, outModel: list, query: str) -> float:
    # create the registry
    registry = accessMatrix()
    # initialize both in and out scores, the probability of the first character i always 0,25
    inscore = numpy.log(0.25)
    outscore = numpy.log(0.25)
    P = query.upper()
    # iterate over every character of the query string
    for i in range(1, len(P)):
        # retrieve the coordinates in the tables
        column_ind = registry[P[i]]
        row_ind = registry[P[i - 1]]

        # obtain the probabilities of a given dimer in both tables
        In_CpG_prob = inModel[row_ind][column_ind]
        Out_CpG_prob = outModel[row_ind][column_ind]

        # update both scores
        inscore += numpy.log(In_CpG_prob)
        outscore += numpy.log(Out_CpG_prob)

    # return the final score
    return round(inscore - outscore, 4)


'''Given the list of the CpGI coordinates, this function returns their total average length.'''


def get_average(L: list) -> int:
    tot = 0
    quantity = len(L)
    for trio in L:
        tot += trio[2]

    return tot // quantity


'''This function returns a DNA sequence randomly chosen and of the desired length.'''


def getRegion(length: int) -> str:
    string = ''.join(random.choice('ACTG') for i in range(length))
    return string


'''This function is able to scan a long region of the genome and by a sliding-window approach it checks whether
    the window considered belongs to a CpGI.'''


def Window_scanning(Inmodel: list, Outmodel: list, coordinates: list, region: str) -> list:
    offsets = []
    # get the average CpGI length
    mean = get_average(coordinates)
    # iterate over each window
    for i in range(len(region) - mean):
        # define the window
        window = region[i:i + mean]
        # find the logarithmic score of the window
        score = queryScore(Inmodel, Outmodel, window)
        if score >= 0:
            offsets.append((i, score))

    # both offset and logarithmic score are returned in the output
    return offsets


'''This function's purpose is to interact with the user. It asks to insert a query string,
    as an output it returns the logarithmic score of the query. The execution may last till 3
    minutes because the function opens both files to retrieve the information needed.'''


def QUERY() -> print:
    q = input('Please insert a query string: ')

    # files are opened and probabilistic models are computed
    chr22 = convertString('chr22.fa')
    CpGIs = get_coordinates('CpGI_chr22')

    INsideModel = InsideModel(chr22, CpGIs)
    OUTsideModel = OutsideModel(chr22, CpGIs)

    # the result is computed
    res = queryScore(INsideModel, OUTsideModel, str(q))

    print('logarithm score: ' + str(res), '\n')
    if res >= 0:
        print('Your query is more likely to be inside a CpG island.')
    else:
        print('Your query is likely to be outside of a CPG island.')


'''This function interacts with the user, it allows either to generate randomly a 2000 nucleotide-long
    sequence or the user is free to insert their own DNA sequence. Then the function performs the sliding-
    window approach to the sequence and it reports both the offsets and the logarithmic score of each window
    that is contained in a CpGI.
    The execution may last till 3 minutes because the function opens both files to retrieve the information needed.'''


def WINDOW() -> print:

    first = input('Do you want to generate a random DNA sequence of 2000 nt? [yes/no]: ')

    if first == 'yes':
        region = getRegion(2000)

    else:
        region = input('Insert a DNA sequence equal or longer than 1000 nucleotides: ')

    chr22 = convertString('chr22.fa')
    CpGIs = get_coordinates('CpGI_chr22')

    INsideModel = InsideModel(chr22, CpGIs)
    OUTsideModel = OutsideModel(chr22, CpGIs)

    # the lists with offsets and scores is computed
    scan = Window_scanning(INsideModel, OUTsideModel, CpGIs, str(region))

    if len(scan) == 0:
        print('NO CpG islands were found.')
    else:
        print(scan)


'''
Examples of the probability models created in the module:

INsideModel =  [['/', 'A', 'C', 'G', 'T'],
               ['A', 0.19, 0.28, 0.4, 0.14],
               ['C', 0.19, 0.36, 0.25, 0.2],
               ['G', 0.17, 0.33, 0.36, 0.14],
               ['T', 0.09, 0.34, 0.38, 0.19]]

OUTsideModel =  [['/', 'A', 'C', 'G', 'T'],
                ['A', 0.29, 0.2, 0.29, 0.22],
                ['C', 0.32, 0.3, 0.07, 0.31],
                ['G', 0.25, 0.24, 0.3, 0.21],
                ['T', 0.18, 0.23, 0.3, 0.29]]
'''

if __name__ == '__main__':

    # execute one of the 2 functions
    # NOTE: the execution may require few minutes

    # WINDOW()
    # QUERY()
