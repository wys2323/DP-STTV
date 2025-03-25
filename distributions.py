import math
import random
import numpy as np
from preprocess import dataset, string, para

sensitivity = 1
lengthMax = 15
speedMax = 17
speedMin = 3
cubeSum = para.cellCount * para.numtimeInterval
path = './data/output/' + dataset + string + '_index.txt'


def sign(t):
    if t < 0:
        return -1.0
    else:
        return 1.0


def LapLaceNoise(epsilon, s):
    d = np.random.random()
    uniform = d - 0.5
    nc = -s / epsilon * sign(uniform) * math.log(1.0 - 2.0 * math.fabs(uniform))
    return nc


def norm(List):
    return [_ / sum(List) for _ in List]


def narrow(List, paraSum, Sum, step):
    countSum = 0
    while countSum < Sum:
        countSum = 0
        for value in List:
            countSum += int(value * paraSum)
        paraSum += step
    paraSum = paraSum - step - step
    return paraSum, countSum


def enforce_consistency(List, Sum, step):
    countSum, paraSum = 0, Sum
    while countSum != Sum and step > pow(0.1, 13):
        paraSum, countSum = narrow(List, paraSum, Sum, step)
        step = step / 10
    paraSum += step * 10
    return paraSum


def startDis(epsilon):
    with open(path) as input:
        content = input.readlines()
    startList = [0 for _ in range(cubeSum)]
    traSum, i = int(len(content) / 2), 1
    while i < len(content):
        line = content[i][3:].split(';')[0:-1]
        startCube = [int(_) for _ in line[0].split(',')]
        startIndex = startCube[0] + startCube[1] * para.cellCount
        startList[startIndex] += 1
        i += 2
    startNoiseList = [0 for _ in range(cubeSum)]
    for i in range(cubeSum):
        noise = LapLaceNoise(epsilon, sensitivity)
        startNoiseList[i] = startList[i] + noise if startList[i] + noise > 0 else 0
    startNoiseList = norm(startNoiseList)
    countSum = enforce_consistency(startNoiseList, traSum, 1000)
    for i in range(cubeSum):
        startNoiseList[i] = int(startNoiseList[i] * countSum)
    return startNoiseList


def markov(epsilon):
    with open(path) as input:
        content = input.readlines()
    matrix, i = [[0 for _ in range(para.cellCount)] for _ in range(para.cellCount)], 1
    while i < len(content):
        line = content[i][3:].split(';')[0:-1]
        line = line[0:-1] if len(line) > 1 else line
        if len(line) > 1:
            j, value = 0, 1 / (len(line) - 1)
            while j < len(line) - 1:
                start = int(line[j].split(',')[0])
                next = int(line[j + 1].split(',')[0])
                matrix[start][next] += value
                j += 1
        i += 2
    matrixNoise = [[0 for _ in range(para.cellCount)] for _ in range(para.cellCount)]
    for i in range(para.cellCount):
        for j in para.neighbor[i]:
            if j == -1:
                break
            if j != i:
                noise = LapLaceNoise(epsilon, sensitivity)
                matrixNoise[i][j] = matrix[i][j] + noise if matrix[i][j] + noise > 0 else 0
        if sum(matrixNoise[i]):
            matrixNoise[i] = norm(matrixNoise[i])
        else:
            neighbors = [x for x in para.neighbor[i] if x != -1]
            if len(neighbors):
                for j in neighbors:
                    if j != i:
                        matrixNoise[i][j] = 1 / (len(neighbors))
    return matrixNoise


def lengthDis():
    with open(path) as input:
        content = input.readlines()
    lengthList = [[] for _ in range(cubeSum)]
    i, maxLength = 1, 0
    while i < len(content):
        line = content[i][3:].split(';')[0:-1]
        line = line[0:-1] if len(line) > 1 else line
        startCube = [int(_) for _ in line[0].split(',')]
        startIndex = startCube[0] + startCube[1] * para.cellCount
        lengthList[startIndex].append(len(line))
        if maxLength < len(line):
            maxLength = len(line)
        i += 2
    para.maxLength = lengthMax
    return lengthList


def speedDis():
    with open(path) as input:
        content = input.readlines()
    speedList = [[[] for _ in range(para.cellCount)] for _ in range(cubeSum)]
    i, maxSpeed = 1, 0
    while i < len(content):
        line = content[i][3:].split(';')[0:-1]
        speed = float(line[-1])
        if speed:
            startCube = [int(_) for _ in line[0].split(',')]
            startIndex = startCube[0] + startCube[1] * para.cellCount
            endCell = int(line[-2].split(',')[0])
            speedList[startIndex][endCell].append(speed)
            if maxSpeed < speed:
                maxSpeed = speed
        i += 2
    return speedList


def lengthLap(traSum, List, epsilon):
    # set
    maxValue = lengthMax + 1
    countLists = [0 for _ in range(1, maxValue)]
    # load
    for value in List:
        if value >= maxValue:
            countLists[maxValue - 2] += 1
        else:
            countLists[value - 1] += 1
    # add noise
    noiseLists = [0 for _ in range(1, maxValue)]
    for i in range(maxValue - 1):
        noise = LapLaceNoise(epsilon, sensitivity)
        noiseLists[i] = countLists[i] + noise if countLists[i] + noise > 0 else 0
    if sum(noiseLists):
        noiseLists = norm(noiseLists)
    else:
        noiseLists = [1 / (maxValue - 1) for _ in range(1, maxValue)]
    countSum = enforce_consistency(noiseLists, traSum, 1000)
    for i in range(maxValue - 1):
        noiseLists[i] = int(noiseLists[i] * countSum)
    return noiseLists


def getIntervalsIndex(intervals, value):
    index = len(intervals)
    for i in range(len(intervals)):
        if value <= intervals[i]:
            index = i
            break
    return index


def speedLap(List, epsilon):
    # set
    increment, flag = 2, 1
    maxValue = speedMax + 1
    intervals, value = [speedMin], speedMin + increment  # minSpeed
    while value < maxValue:
        intervals.append(value)
        value += increment
    dis = [0 for _ in range(len(intervals))]
    # load
    for value in List:
        index = getIntervalsIndex(intervals, value)
        index = index - 1 if index == len(intervals) else index
        dis[index] += 1
    if flag == 1:
        noiseDis = [0 for _ in range(len(intervals))]
        for i in range(len(intervals)):
            noise = LapLaceNoise(epsilon, sensitivity)
            noiseDis[i] = dis[i] + noise if dis[i] + noise > 0 else 0
    else:
        noiseDis = dis
    # norm
    if sum(noiseDis):
        normList = norm(noiseDis)
    else:
        normList = [1 / len(noiseDis) for _ in range(len(intervals))]
    if flag == 1:
        valueSelected = intervals[normList.index(max(normList))]
    else:
        valueSelected = random.choices(intervals, weights=normList, k=1)[0]
    return valueSelected
