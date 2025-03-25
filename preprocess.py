import sys
import math
from haversine import haversine
from parameters import Parameter
from datetime import datetime, timedelta

sys.setrecursionlimit(900)


class Point(object):
    class Struct(object):
        def __init__(self, lon, lat, time, cellIndex, timeIndex, next):
            self.lon = lon
            self.lat = lat
            self.time = time
            self.cellIndex = cellIndex
            self.timeIndex = timeIndex
            self.next = next


def row(index):
    return int(index / para.cellW)


def col(index):
    return int(index - row(index) * para.cellH)


def getNeighbor(cellCount):
    outpath = open('data/parameters/neighborFile_' + str(cellCount) + '.txt', 'w')
    for i in range(0, para.cellCount):
        count = 0
        flag = 0
        row_i, col_i = row(i), col(i)
        for j in range(0, para.cellCount):
            row_j, col_j = row(j), col(j)
            if abs(row_i - row_j) <= 1 and abs(col_i - col_j) <= 1:
                count += 1
                if count <= 8:
                    print("%d" % j, end='', file=outpath)
                    print(' ', end='', file=outpath)
                if count == 9:
                    print("%d" % j, end='', file=outpath)
                    print('', file=outpath)
                    flag = 1
        if flag == 0:
            while count < 9:
                print("-1", end='', file=outpath)
                print(' ', end='', file=outpath)
                count += 1
            print('', file=outpath)
    outpath.close()


def read_tra(tra):
    coordinates = list(map(lambda s: tuple(s.split(',')),
                           filter(lambda l: len(l) > 1, tra.split(';'))))
    return coordinates


def TimeIndex(time):
    start = datetime.combine(time.date(), para.startTime.time())
    delta = (time - start).seconds / 60  # min
    timeIndex = int(delta / para.timestep)
    if (time - (timedelta(minutes=(timeIndex * para.timestep)) + start)).seconds < pow(10, -5):
        timeIndex -= 1
    if timeIndex < 0:
        timeIndex = 0
    return timeIndex


def CellIndex(longitude, latitude):
    incre = 0.0005
    longitude, latitude = float(longitude), float(latitude)
    height = float((para.top - para.bottom) / para.cellH)
    width = float((para.right - para.left) / para.cellW)
    rowIndex = int(math.floor((latitude - para.bottom) / height))
    columnIndex = int(math.floor((longitude - para.left) / width))
    rowIndex1 = int(math.floor((latitude - para.bottom - incre) / height))
    columnIndex1 = int(math.floor((longitude - para.left - incre) / width))
    rowIndex2 = int(math.floor((latitude - para.bottom + incre) / height))
    columnIndex2 = int(math.floor((longitude - para.left + incre) / width))
    if math.fabs(rowIndex * height + para.bottom - latitude) < pow(10, -6):
        rowIndex -= 1
    if math.fabs(columnIndex * width + para.left - longitude) < pow(10, -6):
        columnIndex -= 1
    if rowIndex1 != rowIndex2:
        if rowIndex1 < 0:
            rowIndex = rowIndex2
            latitude += incre
        else:
            rowIndex = rowIndex1
            latitude -= incre
    if columnIndex1 != columnIndex2:
        if columnIndex1 < 0:
            columnIndex = columnIndex2
            longitude += incre
        else:
            columnIndex = columnIndex1
            longitude -= incre

    if rowIndex < 0:
        rowIndex = 0
    if columnIndex < 0:
        columnIndex = 0
    if rowIndex >= para.cellH:
        rowIndex = para.cellH - 1
    if columnIndex >= para.cellW:
        columnIndex = para.cellW - 1
    cellIndex = rowIndex * para.cellW + columnIndex
    if cellIndex >= para.cellH * para.cellW:
        print('something is wrong \n')
    return cellIndex


def CellIndex1(longitude, latitude):
    incre = 0.0005
    longitude = float(longitude)
    latitude = float(latitude)
    height = float((para.top - para.bottom) / para.cellH)
    width = float((para.right - para.left) / para.cellW)
    rowIndex = int(math.floor((latitude - para.bottom) / height))
    columnIndex = int(math.floor((longitude - para.left) / width))
    rowIndex1 = int(math.floor((latitude - para.bottom - incre) / height))
    columnIndex1 = int(math.floor((longitude - para.left - incre) / width))
    rowIndex2 = int(math.floor((latitude - para.bottom + incre) / height))
    columnIndex2 = int(math.floor((longitude - para.left + incre) / width))
    if math.fabs(rowIndex * height + para.bottom - latitude) < pow(10, -6):
        rowIndex -= 1
    if math.fabs(columnIndex * width + para.left - longitude) < pow(10, -6):
        columnIndex -= 1
    if rowIndex1 != rowIndex2:
        if rowIndex1 < 0:
            rowIndex = rowIndex2
            latitude += incre
        else:
            rowIndex = rowIndex1
            latitude -= incre
    if columnIndex1 != columnIndex2:
        if columnIndex1 < 0:
            columnIndex = columnIndex2
            longitude += incre
        else:
            columnIndex = columnIndex1
            longitude -= incre
    if rowIndex < 0:
        rowIndex = 0
    if columnIndex < 0:
        columnIndex = 0
    if rowIndex >= para.cellH:
        rowIndex = para.cellH - 1
    if columnIndex >= para.cellW:
        columnIndex = para.cellW - 1
    cellIndex = rowIndex * para.cellW + columnIndex
    if cellIndex >= para.cellH * para.cellW:
        print('something is wrong \n')
    return cellIndex, longitude, latitude


def PointInsertion(head, tail):
    if head.cellIndex == tail.cellIndex:
        return
    else:
        head_row = int(head.cellIndex / para.cellW)
        head_column = int(head.cellIndex - head_row * para.cellH)
        tail_row = int(tail.cellIndex / para.cellH)
        tail_column = int(tail.cellIndex - tail_row * para.cellH)
        delta_T = (tail.time - head.time).seconds
        if math.fabs(head_row - tail_row) <= 1 and math.fabs(head_column - tail_column) <= 1:
            return
        else:
            mid = Point()
            mid_lon = (head.lon + tail.lon) / 2
            mid_lat = (head.lat + tail.lat) / 2
            mid_time = head.time + timedelta(seconds=int(delta_T / 2))

            mid.cellIndex, mid_lon, mid_lat = CellIndex1(mid_lon, mid_lat)

            mid.lon = mid_lon
            mid.lat = mid_lat
            mid.time = mid_time

            mid.next = tail
            head.next = mid

            PointInsertion(head, mid)
            PointInsertion(mid, tail)
    return


def interpolation():
    ori_out = open('./data/output/' + dataset + '_' + str(para.cellW) + '_out.txt', 'w')
    input_path = './data/output/' + dataset + '.txt'
    with open(input_path) as input:
        content = input.readlines()

    number = 0
    i = 1

    while i < len(content):
        s = content[i][3:]
        tra = read_tra(s)

        flag = 0
        for j in range(len(tra)):
            lon2 = float(tra[j][0])
            lat2 = float(tra[j][1])
            time2 = datetime.strptime(tra[j][2].replace('\n', ''), "%Y-%m-%d %H:%M:%S")
            if j == 0:
                lon1 = lon2
                lat1 = lat2
                time1 = time2

            head = Point()
            head.lon = lon1
            head.lat = lat1
            head.time = time1

            if para.left <= lon1 <= para.right and para.bottom <= lat1 <= para.top:
                head.cellIndex, head.lon, head.lat = CellIndex1(lon1, lat1)
            else:
                head.cellIndex = -1

            tail = Point()
            tail.lon = lon2
            tail.lat = lat2
            tail.time = time2

            if para.left <= lon2 <= para.right and para.bottom <= lat2 <= para.top:
                tail.cellIndex, tail.lon, tail.lat = CellIndex1(lon2, lat2)
            else:
                tail.cellIndex = -1

            if head.cellIndex == -1 or tail.cellIndex == -1:
                continue

            tail.next = None
            head.next = tail
            PointInsertion(head, tail)
            p = head
            while p:
                if p.cellIndex != -1:
                    if flag == 0:
                        if number != 0:
                            print('', file=ori_out)
                        print('#%d:' % number, file=ori_out)
                        print('>0:', end='', file=ori_out)
                        flag = 1
                        number += 1
                    print("%.8f,%.8f,%s;" % (p.lon, p.lat, p.time), end='', file=ori_out)
                p = p.next

            lon1 = tail.lon
            lat1 = tail.lat
            time1 = tail.time

        i += 2
    ori_out.close()


def distance_dif(sp, ep):
    return haversine((sp[1], sp[0]), (ep[1], ep[0]), unit='m')


def time_dif(st, et):
    start_time, end_time = datetime.strptime(st, '%H:%M:%S'), datetime.strptime(et, '%H:%M:%S')
    if end_time >= start_time:
        difference = end_time - start_time
    else:
        end_time += timedelta(days=1)
        difference = end_time - start_time
    return difference.total_seconds()


def remove_duplicates(tra):
    num, pre = 1, tra[0]
    coordinates = [tra[0]]
    while num < len(tra):
        if tra[num][2] != pre[2]:
            coordinates.append(tra[num])
            pre = tra[num]
        num += 1
    return coordinates


def calTraSpeed(tra):
    speedList = []
    for i in range(len(tra)):
        tra[i] = [float(tra[i][0]), float(tra[i][1]), tra[i][2].split(' ')[1]]
    tra = remove_duplicates(tra)
    for j in range(len(tra) - 1):
        distance = distance_dif(tra[j], tra[j + 1])  # m
        duration = time_dif(tra[j][2], tra[j + 1][2])  # s
        speedList.append(round(distance / duration, 2))
    if len(speedList):
        return round(sum(speedList) / len(speedList), 2)
    else:
        return 0


def getIndexSequences():
    with open('./data/output/' + dataset + string + '_index.txt', 'w') as output:
        input_path = './data/output/' + dataset + '_' + str(para.cellW) + '_out.txt'
        with open(input_path) as input:
            content = input.readlines()

        i = 1
        traLen = []
        while i < len(content):
            s = content[i][3:]
            tra = read_tra(s)
            count = 0

            output.write('#' + str(int(i / 2)) + ':\n')
            output.write('>0:')

            start_cell = CellIndex(tra[0][0], tra[0][1])
            time = datetime.strptime(tra[0][2], "%Y-%m-%d %H:%M:%S")
            start_time = TimeIndex(time)
            output.write(str(start_cell) + ',' + str(start_time) + ';')
            count += 1

            if len(tra) > 1:
                j = 1
                last = start_cell
                lastIndex = 0
                cellSpeed = []
                while j < len(tra):
                    cellindex = CellIndex(tra[j][0], tra[j][1])
                    if cellindex != last:
                        output.write(str(cellindex) + ';')
                        speed = calTraSpeed(tra[lastIndex:j])
                        cellSpeed.append(speed)
                        lastIndex = j
                        count += 1
                    last = cellindex
                    if j == len(tra) - 1 and lastIndex != j:
                        speed = calTraSpeed(tra[lastIndex:j])
                        cellSpeed.append(speed)
                    j += 1
                cellSpeed = [x for x in cellSpeed if x != 0]
                avgSpeed = round(sum(cellSpeed) / len(cellSpeed), 2) if cellSpeed else 0
                output.write(str(avgSpeed) + ';')
            output.write('\n')
            traLen.append(count)
            i += 2
    with open('./data/output/' + dataset + '_' + str(para.cellW) + '_traLenFile.txt', 'w', encoding='utf-8') as f:
        traLen = map(lambda x: str(x), traLen)
        for length in traLen:
            f.write(length + '\n')


def preProcess():
    getNeighbor(para.cellW)
    interpolation()
    getIndexSequences()


dataset = 'sample'
para = Parameter()
string = '_' + str(para.cellW) + '_' + str(para.numtimeInterval)
preProcess()
