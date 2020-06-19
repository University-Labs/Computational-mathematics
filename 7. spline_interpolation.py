# Интерполяция сплайнами
import pylab
from math import *


# Построение сплайна дефекта 1
def spline_interpolation(interpolation_x, interpolation_y):
    # Создание сплайна для каждого [xi-1; xi]
    # структура элемента - { x, a, b, c, d }
    spline = []
    for i in range(len(interpolation_x)):
        spline.append([-1, 0, 0, 0, 0])
        spline[i][0] = interpolation_x[i]

    # Находим коэфф-ты ai ( Si(xi) = ai )
    for i in range(len(interpolation_x)):
        spline[i][1] = interpolation_y[i]

    # c0 = 0, cn = 0
    spline[0][3] = 0
    spline[len(interpolation_x) - 1][3] = 0

    # Нахождение ci методом прогонки
    alpha = [0 for _ in range(len(interpolation_x) - 1)]
    betta = [0 for _ in range(len(interpolation_x) - 1)]

    for i in range(1, len(interpolation_x) - 1):
        hi = interpolation_x[i] - interpolation_x[i - 1]        # hi
        hi1 = interpolation_x[i + 1] - interpolation_x[i]       # hi+1

        # приведение уравнения к виду, пригодному для решения методом прогонки
        a = hi
        b = -(2 * (hi + hi1))
        c = hi1
        d = 6 * ((interpolation_y[i + 1] - interpolation_y[i]) / hi1 - (interpolation_y[i] - interpolation_y[i - 1])
                 / hi)

        alpha[i] = c / (b - a * alpha[i - 1])
        betta[i] = (a * betta[i - 1] - d) / (b - a * alpha[i - 1])
    for i in range(len(interpolation_x) - 2, 0, -1):        # обратный ход метода прогонки
                                                            # нахождение c[i] с конца на начало
        spline[i][3] = alpha[i] * spline[i + 1][3] + betta[i]

    # Нахождение bi и di
    for i in range(len(interpolation_x) - 1, 0, -1):
        h = interpolation_x[i] - interpolation_x[i - 1]
        spline[i][4] = (spline[i][3] - spline[i - 1][3]) / h
        spline[i][2] = (h / 2 * spline[i][3]) - \
                       (h ** 2 / 6 * spline[i][4]) + (interpolation_y[i] - interpolation_y[i - 1]) / h
    return spline


# Вычисление сплайна и значения функции по сплайну в заданной точке
def spline_values(splines, value_x):
    size = len(splines)
    if value_x <= splines[0][0]:            # значение x слева от минимального сплайна
        curr_spline = splines[0]
    elif value_x >= splines[size - 1][0]:   # значение x справа от максимального сплайна
        curr_spline = splines[size - 1]
    else:                                   # значение x в одном из внутренних Si
        # бинарный поиск Si, в который входит value_x
        min_c = 0
        max_c = size - 1
        while min_c + 1 < max_c:
            middle = (min_c + max_c) // 2
            if value_x <= splines[middle][0]:
                max_c = middle
            else:
                min_c = middle
        curr_spline = splines[max_c]
    # Значение сплайна в этой точке
    return curr_spline[1] + (curr_spline[2] * (value_x - curr_spline[0])) +\
           (curr_spline[3] / 2 * (value_x - curr_spline[0]) ** 2) + \
           (curr_spline[4] / 6 * (value_x - curr_spline[0]) ** 3)


def main():
    # Имя файла по умолчанию

    print("Введите 1, чтобы программа определяла значение функции в заданной точке")
    print("Введите 2, чтобы программа выводила по наименованию функции графики функции и сплайна")
    print("Введите 3, чтобы программа выводила по таблице графики сплайна")
    print("Введите 4, чтобы программа выводила кривую, заданную параметрически ")
    des = ''  # решение пользователя по режиму работы программы
    while True:
        des = input()
        if des == '1' or des == '2' or des == '3' or des == '4':
            break

    if des == '1':  # определяется значение функции в заданной точке
        interpolation_x = []                # Узлы интерполяции
        interpolation_y = []                # Значение функции в узлах интерполяции
        value_x = 0                         # Точка, в которой нужно найти приближенное значение
        value_y = 0                         # Значение функции в этой точке
        with open('input1.txt', "r") as inp:          # считывание входных данных из файла
            for line in inp.readlines():
                list_m = line.split()
                if len(list_m) == 2:
                    interpolation_x.append(float(list_m[0]))
                    interpolation_y.append(float(list_m[1]))
                else:
                    value_x = float(list_m[0])

        # Вывод значений
        print('Табличные значения: ')
        print('\t\tX\t\t\tY')
        for i in range(len(interpolation_x)):
            print('{0})\t{1:^10}\t{2:^10}'.format(i + 1, interpolation_x[i], interpolation_y[i]))

        # построение сплайна
        splines = spline_interpolation(interpolation_x, interpolation_y)
        # значение функции в заданной точке
        value_y = spline_values(splines, value_x)

        # Печать результата
        print("\nf({0}) = {1}\n".format(value_x, value_y))

    elif des == '2':  # Строится таблица значений по названию функции, выводится график этой функции и инт. сплайна
        function_name = ""  # Введённая функция
        interpolation_x = []  # Узлы интерполирования
        interpolation_y = []  # Значения функции
        x_list = []  # Список значений аргументов функции для красивого вывода на графике

        # Считывание из файла названия функции и координат x, в которых необходимо её исследовать
        with open('input2.txt', 'r') as inp:
            function_name = inp.readline().replace('\n', '')
            interpolation_x = list(map(float, inp.readline().split()))

        # Нахождение значений функции программно
        for t_values in interpolation_x:
            x = t_values
            interpolation_y.append(eval(function_name))

        # Вычисление списка аргументов, лежащих внутри отрезка [x_min; x_max]
        x_min = interpolation_x[0]  # Минимальный x
        x_max = interpolation_x[len(interpolation_x) - 1]  # Максимальный x

        tmp_val = x_min
        while tmp_val <= x_max:
            x_list.append(tmp_val)
            tmp_val += 0.001

        # Значения введённой функции f(x) в каждой точке отрезка [x_min; x_max]
        func_value = [eval(function_name) for x in x_list]
        # Построение сплайна
        splines = spline_interpolation(interpolation_x, interpolation_y)
        # Значения интерполяционного многочлена в каждой точке отрезка [x_min; x_max]
        splines_value = [spline_values(splines, x) for x in x_list]

        # Печать в консоль значений функции в узлах интерполяции
        print('Таблица значений функции {0}: '.format(function_name))
        print('\t\tX\t\t\tY')
        for i in range(len(interpolation_x)):
            print('{0})\t{1:^10}\t{2:^10.4f}'.format(i + 1, interpolation_x[i], interpolation_y[i]))

        # Отрисовка обоих графиков в одной системе координат
        pylab.plot(x_list, func_value)
        pylab.plot(x_list, splines_value)
        # Открытие окна с графиками
        pylab.show()

    elif des == '3':  # Построение по таблице значений сплайна и его графика
        interpolation_x = []  # узлы интерполяции
        interpolation_y = []  # значение функции в узлах
        x_list = []  # Список значений аргументов функции для красивого вывода на графике

        with open('input3.txt', "r") as inp:  # считывание входных данных из файла
            for line in inp.readlines():
                list_m = line.split()
                interpolation_x.append(float(list_m[0]))
                interpolation_y.append(float(list_m[1]))

        # Вычисление списка аргументов, лежащих внутри отрезка [x_min; x_max]
        x_min = interpolation_x[0]  # Минимальный x
        x_max = interpolation_x[len(interpolation_x) - 1]  # Максимальный x

        tmp_val = x_min
        while tmp_val <= x_max:
            x_list.append(tmp_val)
            tmp_val += 0.001

        # Построение сплайна
        splines = spline_interpolation(interpolation_x, interpolation_y)
        # Значения функции в каждой точке отрезка [x_min; x_max]
        splines_value = [spline_values(splines, x) for x in x_list]

        # Печать в консоль значений функции в узлах интерполяции
        print('Таблица значений')
        print('\t\tX\t\t\tY')
        for i in range(len(interpolation_x)):
            print('{0})\t{1:^10}\t{2:^10.4f}'.format(i + 1, interpolation_x[i], interpolation_y[i]))

        # Построение графика и открытие окна отображения
        pylab.plot(x_list, splines_value)
        # Открытие окна с графиком
        pylab.show()

    elif des == '4':  # Построение графика кривой заданной параметрически
        value_t = []                # значения параметра t
        interpolation_x = []        # x-координата кривой в точке ti
        interpolation_y = []        # y-координата кривой в точке ti

        with open("input4.txt", "r") as inp:  # считывание входных данных из файла
            for line in inp.readlines():
                list_m = line.split()
                interpolation_x.append(float(list_m[0]))
                interpolation_y.append(float(list_m[1]))

        #value_t = [i for i in range(len(interpolation_x))]

        value_t.append(0)
        tmp_val = 0
        for i in range(len(interpolation_x) - 1):
            #tmp_val += sqrt((interpolation_x[i + 1] - interpolation_x[i]) ** 2 +
            #                (interpolation_y[i + 1] - interpolation_y[i]) ** 2)
            tmp_val += abs(interpolation_x[i+1] - interpolation_x[i]) + abs(interpolation_y[i+1] - interpolation_y[i])
            value_t.append(tmp_val)

        t_list = []
        tmp_val = 0
        max_v = value_t[len(value_t) - 1]
        while tmp_val <= max_v:
            t_list.append(tmp_val)
            tmp_val += 0.01

        # Вывод  координат точек
        print("Таблица значений кривой:")
        print('\t\tti\t\t\tXi\t\t\tYi')
        for i in range(len(interpolation_x)):
            print('{0})\t{1:^10.4f}\t{2:^10.4f}\t{3:^10.4f}'
                  .format(i + 1, value_t[i], interpolation_x[i], interpolation_y[i]))

        # Вычисление сплайна для x = x(t)
        splines_x = spline_interpolation(value_t, interpolation_x)
        values_x = [spline_values(splines_x, x) for x in t_list]

        # Вычисление сплайна для y = y(t)
        splines_y = spline_interpolation(value_t, interpolation_y)
        values_y = [spline_values(splines_y, y) for y in t_list]

        # Рисование графика
        pylab.plot(t_list, values_x)
        pylab.plot(t_list, values_y)
        # Открытие окна с графиком
        pylab.show()

    else:
        print("Такой программы нет!")


if __name__ == '__main__':
    main()
