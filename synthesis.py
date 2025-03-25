import time
import random
import datetime
from preprocess import para, dataset, string, distance_dif
from distributions import startDis, lengthDis, speedDis, markov, cubeSum, lengthLap, speedLap


def getNeighbors(cell):
    neighbors = []
    for value in para.neighbor[cell]:
        if value == -1:
            break
        if value != cell:
            neighbors.append(value)
    return neighbors


def transLoc(index):
    rowIndex = int(int(index) / para.cellW)
    columnIndex = int(int(index) % para.cellW)
    height = float((para.top - para.bottom) / para.cellH)
    width = float((para.right - para.left) / para.cellW)
    lat = float(para.bottom + rowIndex * height + height / 2)
    lon = float(para.left + columnIndex * width + width / 2)
    return lon, lat


def transTime(index):
    return para.startTime + datetime.timedelta(seconds=para.timestep * float(index) * 60)


def trans(cellTra, timeIndex, speed):
    startLon, startLat = transLoc(cellTra[0])
    startTime = transTime(timeIndex)
    tra = [[startLon, startLat, str(startTime)]]
    i, preLoc = 1, [startLon, startLat]
    while i < len(cellTra):
        if cellTra[i - 1] != cellTra[i]:
            currentLoc = transLoc(cellTra[i])
            routeLength = distance_dif(preLoc, currentLoc)
            currentTime = startTime + datetime.timedelta(seconds=int(routeLength / speed))
            tra.append([currentLoc[0], currentLoc[1], str(currentTime)])
            preLoc = currentLoc
            startTime = currentTime
        i += 1
    return tra


rate_start, rate_markov, rate_length, rate_speed = 1, 1, 1, 1
rate_sum = rate_start + rate_markov + rate_length + rate_speed


# start(cube)、length(cube)、markov(cell)、speed(cube,cell)
def synthesisTra(epsilon, run):
    # set epsilon
    eps_start, eps_markov = rate_start / rate_sum * epsilon, rate_markov / rate_sum * epsilon
    eps_speed, eps_length = rate_speed / rate_sum * epsilon, rate_length / rate_sum * epsilon
    # generate trajectories
    numTra = 0
    output = open('./data/output/syn/' + dataset + string + '_' + str(epsilon)
                  + '_' + str(run) + '_tra.txt', 'w', encoding='utf-8')
    output_loc = open('./data/output/syn/' + dataset + string + '_' + str(epsilon)
                      + '_' + str(run) + '_loc.txt', 'w', encoding='utf-8')
    cStartDis, cMarkovDis = startDis(eps_start), markov(eps_markov)
    cLengthDis, cSpeedDis = lengthDis(), speedDis()
    lengthList = [_ for _ in range(1, para.maxLength + 1)]
    for i in range(cubeSum):
        if cStartDis[i]:
            countList = lengthLap(cStartDis[i], cLengthDis[i], eps_length)
            for j in range(len(countList)):
                if countList[j]:
                    startCell, startTime, traLen = i % para.cellCount, int(i / para.cellCount), lengthList[j]
                    for k in range(countList[j]):
                        cellTra = [startCell]
                        if traLen > 2:
                            for n in range(traLen - 1):
                                neighbors = getNeighbors(cellTra[n])
                                ngbProbs = [cMarkovDis[cellTra[n]][_] for _ in neighbors]
                                if len(neighbors):
                                    currentCell = random.choices(neighbors, weights=ngbProbs, k=1)[0]
                                    cellTra.append(currentCell)
                                else:
                                    break
                        speed = 0
                        while speed == 0:
                            speed = speedLap(cSpeedDis[i][cellTra[-1]], eps_speed)
                        tra = trans(cellTra, startTime, speed)
                        output.write('#' + str(numTra) + ':\n' + '>0:')
                        output_loc.write('#' + str(numTra) + ':\n' + '>0:')
                        for point in tra:
                            p_tra = ','.join([str(_) for _ in point]) + ';'
                            p_loc = ','.join([str(_) for _ in point[0:2]]) + ';'
                            output.write(p_tra)
                            output_loc.write(p_loc)
                        output.write('\n')
                        output_loc.write('\n')
                        numTra += 1
    output.close()
    output_loc.close()


times = 1
for eps in [1]:
    for j in range(times):
        synthesisTra(eps, j + 1)
