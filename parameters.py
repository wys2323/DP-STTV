import os
import sys
import math
import pandas as pd
from datetime import datetime

sys.setrecursionlimit(6000)


def boundary():
    l_boundary = open('./data/parameters/boundary.txt', "r")
    if l_boundary is None:
        print("Can not open boundary file\n")
        return
    lines, data = l_boundary.readlines(), []
    for line in lines:
        data.append(line.split(' '))
    left = float(data[0][0])
    right = float(data[0][1])
    bottom = float(data[0][3])
    top = float(data[0][2])
    l_boundary.close()
    return left, right, bottom, top


def time():
    t_boundary = open('./data/parameters/time.txt', "r")
    lines, data = t_boundary.readlines(), []
    for line in lines:
        data.append(line.replace('\n', ''))
    time_boundary = pd.to_datetime(data)
    start_time = time_boundary[0].hour
    end_time = time_boundary[1].hour - 1
    start_hour = str(time_boundary[0])
    end_hour = str(time_boundary[1])
    start = start_hour.split(' ')
    end = end_hour.split(' ')
    startTime = datetime.strptime(start[1], '%H:%M:%S')
    endTime = datetime.strptime(end[1], '%H:%M:%S')
    t_boundary.close()
    return startTime, endTime


def TimeInterval():
    f_timeInterval = open('./data/parameters/timeInterval.txt', "r")
    if f_timeInterval is None:
        print("Can not open timeInterval file\n")
        return
    lines, data = f_timeInterval.readlines(), []
    for line in lines:
        data.append(line)
    numTimeInterval = int(data[0])
    f_timeInterval.close()
    return numTimeInterval


def cellSize():
    f_cellSize = open('./data/parameters/cellSize.txt', "r")
    if f_cellSize is None:
        print("Can not open cellSize file\n")
        return
    lines, data = f_cellSize.readlines(), []
    for line in lines:
        data.append(line.split(' '))
    cellH = int(data[0][0])
    cellW = int(data[0][1])
    f_cellSize.close()
    return cellH, cellW


def neighborFile():
    n_boundary = open('data/parameters/neighborFile_' + str(cellSize()[0]) + '.txt', 'a+')
    path = 'data/parameters/neighborFile_' + str(cellSize()[0]) + '.txt'
    if not os.path.getsize(path):
        pass
    else:
        n_boundary.close()
        n_boundary = open(path, 'r')
        if n_boundary is None:
            print("Can not open file\n")
            return
        lines, neighbor = n_boundary.readlines(), []
        for line in lines:
            count, data = 0, []
            while count < 9:
                data.append(int(line.split(' ')[count]))
                count += 1
            neighbor.append(data)
        return neighbor


class Parameter(object):
    def __init__(self):
        self.left = boundary()[0]
        self.right = boundary()[1]
        self.top = boundary()[2]
        self.bottom = boundary()[3]
        self.neighbor = neighborFile()
        self.startTime = time()[0]
        self.endTime = time()[1]
        self.cellH = cellSize()[0]
        self.cellW = cellSize()[1]
        self.numtimeInterval = TimeInterval()
        self.cellCount = self.cellW * self.cellH
        self.timestep = math.ceil((self.endTime - self.startTime).seconds / self.numtimeInterval / 60)  # min
        self.maxLength = 15
        self.epsilon = 1
        self.sensitivity = 1

    def show(self):
        print('********** 参数')
        print('left : %f' % self.left, '\t right : %f' % self.right, '\t bottom : %f' % self.bottom,
              '\t top : %f' % self.top)
        print('cellH : %d' % self.cellH, '\t cellW : %d' % self.cellW, '\t cellCount : %d' % self.cellCount)
        print('startTime : %s' % self.startTime, '\t endTime : %s' % self.endTime)
        print('numTimeInterval : %d' % self.numtimeInterval, '\t timestep : %d' % self.timestep + ' min')
        print('maxTLen : %d' % self.maxLength, '\t epsilon : %s' % self.epsilon)
        print('********** end')
