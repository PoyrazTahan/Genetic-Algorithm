import Node
import Model
import math
from random import randint
import time
import matplotlib.pyplot as plt



stop_sugesstion = 300
population_size = 100
turned_on_gene_number = 10
mutation_percentage = 10
number_of_breedingAlgo = 2

#  irrelevant variables to the user
breeding_algo = 0
is_mutation = 0  # 0 for no mutation 1 for static mutation 2 for dynamic
stagnation = 0
unique_gene_number = 0
gene_list = []
distancesMatrix = []
individuals_list = []


def readGeneProperties(path):
    '''
    :param path:
    :return:
    '''
    global unique_gene_number

    counter = 0
    with open(path, 'r') as file:
        for line in file:
            n = Node.Node()
            for word in line.split(";"):
                if counter % 4 is 0:
                    n.index = int(word)  # set index of the current Node
                elif counter % 4 is 1:
                    n.demand = int(word)  # set X coordinate of the current Node
                elif counter % 4 is 2:
                    n.x = int(word)  # set X coordinate of the current Node
                else:

                    n.y = int(word)

                counter += 1
            gene_list.append(n)

    unique_gene_number = counter/4
    file.close()




def distanceMatrix():
    """
    This function creates distance matrix of the global nodeList
    :return:
    """

    for i in range(len(gene_list)):
        tmp = []
        for j in range(len(gene_list)):
            distance = math.sqrt(
                ((int(gene_list[i].x) - int(gene_list[j].x)) ** 2) + ((int(gene_list[i].y) - int(gene_list[j].y)) ** 2))
            # euclidian distance: sqrt( (x1-x2)^2 + (y1-y2)^2 )
            tmp.append(distance)
        distancesMatrix.append(tmp)


def initialGenerationCreation():
    """
    creates the first random generation
    :return:
    """

    for i in range(population_size):
        locatedModel = Model.LocatedModel(unique_gene_number)  # create default locatedModel

        #   define 1 chromosomes in DNA of current locatedModel randomly untill ones is equal to the pmedian
        counter = 0
        while counter < turned_on_gene_number:
            randomIndex = randint(0, len(locatedModel.DNA) - 1)
            if (locatedModel.DNA[randomIndex] is not 1):
                locatedModel.DNA[randomIndex] = 1
                counter += 1

        calcFitness(locatedModel)  # calculate fitness value for the current locatedModel
        individuals_list.append(locatedModel)  # Create Candidate goes to the List of Generation


def calcFitness(locatedModel):
    """

    :param locatedModel:
    :return:
    """
    sum = 0
    candidates = locatedModel.getOnes()  # get opened Node indexes which is based on 1 chromosomes in DNA of the locatedModel

    for i in range(len(locatedModel.DNA)):
        if (locatedModel.DNA[i] is 0):
            closestCandidateDistance = closestDitanceFinder(i, locatedModel, candidates[
                0])  # candidates[0] is used as reference to find closest Candidate Node
            sum = sum + (float(gene_list[i].demand) * float(closestCandidateDistance))

    locatedModel.fitness = sum


def closestDitanceFinder(demandNodeIndex, locatedModel, referenceCandidate):
    """
    finds the closest visited city to the considered city

    :param demandNodeIndex:
    :param locatedModel:
    :param referenceCandidate:
    :return:
    """

    closestDistance = distancesMatrix[demandNodeIndex][referenceCandidate]  # gives an initial distance

    # search whether there is any closer candidate Node to the current demand Node
    for i in locatedModel.getOnes():
        if ((float(distancesMatrix[demandNodeIndex][i]) - float(closestDistance) < 0)):
            closestDistance = distancesMatrix[demandNodeIndex][i]

    return closestDistance


# Iterations Methods
def naturalSelection():
    """
    start of the genetic process, natural selection

    :return:
    """
    protectedPercentage = (population_size * 30) / 100  # first 30% of generation will be saved
    startRemovingIndex = population_size - protectedPercentage  # last 30% of generation will be removed

    #   remove last 30% of locatedModelList[] / least fit individuals
    for i in range(population_size - 1, (startRemovingIndex - 1), -1):
        # starts from LAST index and remove locatedModel backwards (solution of index out of boundry issue)
        individuals_list.remove(individuals_list[i])

    populationAfterElection = len(individuals_list)  # get current population after removing last 30%

    # randomly remove last 40% of elected generation
    for i in range(populationAfterElection - 1, protectedPercentage, -1):
        r = randint(0, 1)
        if (r is 0):
            individuals_list.remove(individuals_list[i])

    reproduce()  # reproduce the generation with new baby locatedModels


def reproduce():
    """
    replenish the generation, selection of parent individuals
    :return:
    """
    currentPopulation = len(individuals_list)

    #   start from current population and reproduce baby till current population is equal to the 'population' parameter
    for i in range(currentPopulation, population_size):
        dad = individuals_list[randint(0, currentPopulation - 1)]  # choose a dad randomly from current population
        mom = individuals_list[randint(0, currentPopulation - 1)]  # choose a mom randomly from current population
        while dad is mom:
            mom = individuals_list[randint(0, currentPopulation - 1)]

        if breeding_algo is 1:
            individuals_list.append(breeding1(dad, mom))  # after get the new baby add it into populationArray[]
        elif breeding_algo is 2:
            individuals_list.append(breeding2(dad, mom))
        else:
            pass


def breeding1(dad, mom):
    """
    breeding algorithm

    :param dad:
    :param mom:
    :return:
    """
    baby = Model.LocatedModel(unique_gene_number)  # create a default  baby

    #   UNIFORM CROSSOVER chromosomes are chosen radomly form dad and mom
    for i in range(len(baby.DNA)):
        randomChromosome = randint(0, 1)
        if randomChromosome is 0:
            baby.DNA[i] = dad.DNA[i]
        else:
            baby.DNA[i] = mom.DNA[i]

    #   In this part fix amount of ones in DNA if there are more or less than pmedian
    onesIndexes = baby.getOnes()

    difference = turned_on_gene_number - len(
        onesIndexes)  # difference<0 means number of ones more than pmedian so we need to remove some ones randomly
    if difference < 0:
        for i in range(abs(difference)):
            r = randint(0, len(onesIndexes) - 1)
            randomDnaIndex = onesIndexes[r]
            baby.DNA[randomDnaIndex] = 0
            onesIndexes.remove(onesIndexes[r])
    elif difference > 0:
        dnaLength = len(baby.DNA)
        i = 0
        while i is not difference:
            randomDnaIndex = randint(0, dnaLength - 1)
            if (baby.DNA[randomDnaIndex] is 0):
                baby.DNA[randomDnaIndex] = 1
                i += 1
    else:
        pass

    if is_mutation == 0:
        pass
    elif is_mutation == 1:
        mutationStatic(baby)  # check mutation factor for baby
    else:
        mutationDynamic(baby)  # check mutation factor for baby

    calcFitness(baby)  # get fitness value of baby

    return baby


def breeding2(dad, mom):
    """
    breeding algorithm

    :param dad:
    :param mom:
    :return:
    """
    baby = Model.LocatedModel(unique_gene_number)  # create a default  baby

    while True:
        for i in range(len(baby.DNA)):
            randomChromosome = randint(0, 1)
            if randomChromosome is 0:
                baby.DNA[i] = dad.DNA[i]
            else:
                baby.DNA[i] = mom.DNA[i]

        #   In this part fix amount of ones in DNA if there are more or less than pmedian
        onesIndexes = baby.getOnes()

        difference = turned_on_gene_number - len(
            onesIndexes)  # difference<0 means number of ones more than pmedian so we need to remove some ones randomly
        if difference == 0:
            break

    if is_mutation == 0:
        pass
    elif is_mutation == 1:
        mutationStatic(baby)  # check mutation factor for baby
    else:
        mutationDynamic(baby)  # check mutation factor for baby

    calcFitness(baby)  # get fitness value of baby

    return baby


def mutationStatic(baby):
    onesIndexes = baby.getOnes()

    mutationFactor = randint(0, 100)  # check mutation 10% probability

    if mutationFactor < mutation_percentage:  # if mutationFactor less than 10 then turn one of the 1s into 0 and vice versa
        muttationChromosome = onesIndexes[randint(0, len(onesIndexes) - 1)]
        baby.DNA[muttationChromosome] = 0

        randomChromosome = randint(0, len(baby.DNA) - 1)

        # check if random chromosome is 1 then choose another one till chosen chromosome is 0 to turn into 1
        while baby.DNA[randomChromosome] is 1:
            if baby.DNA[randomChromosome] is 0:
                break
            randomChromosome = randint(0, len(baby.DNA) - 1)

        baby.DNA[randomChromosome] = 1


def mutationDynamic(baby):
    onesIndexes = baby.getOnes()
    zeroIndexes = baby.getZeros()

    if stagnation < 40:
        probMuatate = mutation_percentage
    else:
        probMuatate = mutation_percentage + 0.15 * stagnation

    mutationFactor = randint(0, 100)  # check mutation 10% probability

    if mutationFactor < probMuatate:  # if mutationFactor less than 10 then turn one of the 1s into 0 and vice versa
        muttationChromosome = [onesIndexes[randint(0, len(onesIndexes) - 1)]]
        complementryList = [zeroIndexes[randint(0, len(onesIndexes) - 1)]]

        if mutationFactor < probMuatate / 3:
            muttationChromosome.append(onesIndexes[randint(0, len(onesIndexes) - 1)])
            while muttationChromosome[1] == muttationChromosome[0]:
                muttationChromosome[1] = onesIndexes[randint(0, len(onesIndexes) - 1)]

            complementryList.append(zeroIndexes[randint(0, len(zeroIndexes) - 1)])
            while complementryList[1] == complementryList[0]:
                complementryList[1] = zeroIndexes[randint(0, len(zeroIndexes) - 1)]

            if mutationFactor < probMuatate / 10:
                muttationChromosome.append(onesIndexes[randint(0, len(onesIndexes) - 1)])
                while muttationChromosome[0] == muttationChromosome[2] or muttationChromosome[1] == muttationChromosome[
                    2]:
                    muttationChromosome[2] = onesIndexes[randint(0, len(onesIndexes) - 1)]

                complementryList.append(zeroIndexes[randint(0, len(zeroIndexes) - 1)])
                while complementryList[1] == complementryList[2] or complementryList[0] == complementryList[2]:
                    complementryList[2] = zeroIndexes[randint(0, len(zeroIndexes) - 1)]

        for i in muttationChromosome:
            baby.DNA[i] = 0

        for i in complementryList:
            baby.DNA[i] = 1


# ---------HELPER FUNCTIONS TO SORT GENERATIONS by QuickSort Algorithm-----------
def sortList(alist):
    quickSortHelper(alist, 0, len(alist) - 1)


def quickSortHelper(alist, first, last):
    if first < last:
        splitpoint = partition(alist, first, last)

        quickSortHelper(alist, first, splitpoint - 1)
        quickSortHelper(alist, splitpoint + 1, last)


def partition(alist, first, last):
    pivotvalue = alist[first].fitness

    leftmark = first + 1
    rightmark = last

    done = False
    while not done:

        while leftmark <= rightmark and alist[leftmark].fitness <= pivotvalue:
            leftmark = leftmark + 1

        while alist[rightmark].fitness >= pivotvalue and rightmark >= leftmark:
            rightmark = rightmark - 1

        if rightmark < leftmark:
            done = True
        else:
            temp = alist[leftmark]
            alist[leftmark] = alist[rightmark]
            alist[rightmark] = temp

    temp = alist[first]
    alist[first] = alist[rightmark]
    alist[rightmark] = temp

    return rightmark


def printInitialGeneration():
    print "\n\n_*_  The Best Fitness Value AT THE BEGINNING: " + str(individuals_list[0].fitness) + "  _*_"
    print "--LocatedModel Info--"
    individuals_list[0].toString()


def printResults(time):
    print "\n\n_*_  The Best Fitness Value AT THE END: " + str(individuals_list[0].fitness) + "  _*_"
    print "--LocatedModel Info--"
    individuals_list[0].toString()
    print "time: " + str(time)


def main():
    global individuals_list
    global breeding_algo
    global is_mutation

    # ----Process Preparation----#
    readGeneProperties('./TurkeyP1.csv')

    distanceMatrix()

    breeding_algo = 0
    array_generations = []
    array_avg_fitness = []


    record_best_fitness = []
    record_generation = []
    record_time = []
    record_avg_fitness = []

    for i in range(number_of_breedingAlgo):
        breeding_algo = i + 1
        for j in range(3):
            is_mutation = j

            start_time = time.time()

            individuals_list = []

            initialGenerationCreation()

            sortList(
                individuals_list)  # ilk neslin en iyi fitness degerli locatedmodelini ekrana basmak icin sort ettik
            printInitialGeneration()

            stagnation = 0
            lastBest = individuals_list[0].fitness
            # ----Process----#
            counter = 1
            stagnation = 0
            mix = 0

            while mix + stagnation < stop_sugesstion:  # sort, remove, reproduce

                naturalSelection()
                sortList(individuals_list)  # sort locatedModels by their fitness values in locatedModelList[]

                if lastBest == individuals_list[0].fitness:
                    stagnation = stagnation + 1
                else:
                    stagnation = 0

                print  "Generation: " + str(counter) + " \tbest fitness: " + str(
                    individuals_list[0].fitness) + "\tstagnation: " + str(stagnation)
                lastBest = individuals_list[0].fitness
                record_generation.append(lastBest)
                
                sum_fitness = 0
                for i in individuals_list:
                    sum_fitness += i.fitness
                record_avg_fitness.append(sum_fitness/population_size)
                counter = counter + 1
                if mix > 100 and stagnation == 0:
                    mix = mix - 50

                mix = mix + 1

            sortList(individuals_list)

            elapsed_time = time.time() - start_time
            record_time.append(elapsed_time)
            record_best_fitness.append(individuals_list[0].fitness)



           # average_fitness = sum_fitness/population_size
            
            printResults(elapsed_time)

            array_generations.append(record_generation)
            array_avg_fitness.append(record_avg_fitness)
            

    print "time for each algo: " + str(record_time)
    print "solution for each algo: " + str(record_best_fitness)

    
    plt.figure(1)
    plt.title('Evolvement of the each algorithm')
    for i in array_generations:
        plt.plot(i)
    plt.show()
    plt.figure(2)
    plt.title('Average Evolvement of the each algorithm')
    for i in array_avg_fitness:
        plt.plot(i)
    plt.show()
    plt.figure(3)
    plt.title('time spent in each Algorithm')
    plt.bar(range(6), record_time)
    plt.show()
    plt.figure(4)
    plt.title('Solution Comparison')
    plt.bar(range(6), record_best_fitness)
    plt.show()
    


    


if __name__ == '__main__':
    main()
