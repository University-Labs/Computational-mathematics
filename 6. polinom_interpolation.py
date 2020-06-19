# Интерполяция многочленами
import sys
import pylab
import matplotlib as mlab
from math import *


# Нахождение значения y_val интерполяционного многочлена Лагранжа
# с узлами points, значениями в этих узлах value_points
# в точке x_val
def polinom_lagrange(points, value_points, x_val):
    y_val = 0  # Итоговый результат
    amount_points = len(points)

    for k in range(amount_points):  # Внешняя сумма
        top_multy = 1
        bot_multy = 1
        for j in range(amount_points):  # Произведение в числителе и знаменателе
            if j != k:
                top_multy *= (x_val - points[j])
                bot_multy *= (points[k] - points[j])

        y_val += (top_multy / bot_multy) * value_points[k]
    return y_val


def main():
    print("Введите 1, чтобы программа определяла значение функции в заданной точке")
    print("Введите 2, чтобы программа выводила графики инт. многочлена и функции")
    des = ''  # решение пользователя по режиму работы программы
    while True:
        des = input()
        if des == '1' or des == '2':
            break

    if des == '1':
        interpolation_x = []  # Узлы интерполяции
        interpolation_y = []  # Значение функции в этих узлах
        value_x = 0  # Точка x, в которой нужно найти значение функции
        value_y = 0  # Значение y, которое нужно найти

        # Считывание из файла пар (x,y) и в конце значение x, для которого нужно найти y
        with open('input1.txt', 'r') as inp:
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
            print('{0})\t{1:^10}\t{2:^10}'.format(i+1, interpolation_x[i], interpolation_y[i]))

        # Нахождение значения функции в точке value_x
        value_y = polinom_lagrange(interpolation_x, interpolation_y, value_x)

        print('_____________________\nОтвет: f({0}) = {1}'
              '\n_____________________'.format(value_x, value_y))

    elif des == '2':
        function_name = ""    # Введённая функция
        interpolation_x = []  # Узлы интерполирования
        interpolation_y = []  # Значения функции
        x_list = []  # Список значений аргументов функции для красивого вывода на графике

        # Считывание из файла названия функции и координат x, в которых необходимо её исследовать
        with open('input2.txt', 'r') as inp:
            function_name = inp.readline().replace('\n', '')
            interpolation_x = list(map(float, inp.readline().split()))

        # Нахождение значений функции программно
        for arg in interpolation_x:
            x = arg
            interpolation_y.append(eval(function_name))

        # Вычисление списка аргументов, лежащих внутри отрезка [x_min; x_max]
        x_min = interpolation_x[0]  # Минимальный x
        x_max = interpolation_x[len(interpolation_x) - 1]  # Максимальный x

        tmp_val = x_min
        while tmp_val < x_max:
            x_list.append(tmp_val)
            tmp_val += 0.001

        # Значения введённой функции f(x) в каждой точке отрезка [x_min; x_max]
        func_value = [eval(function_name) for x in x_list]
        # Значения интерполяционного многочлена в каждой точке отрезка [x_min; x_max]
        poli_value = [polinom_lagrange(interpolation_x, interpolation_y, x) for x in x_list]

        # Печать в консоль значений функции в узлах интерполяции
        print('Таблица значений функции {0}: '.format(function_name))
        print('\t\tX\t\t\tY')
        for i in range(len(interpolation_x)):
            print('{0})\t{1:^10}\t{2:^10.4f}'.format(i + 1, interpolation_x[i], interpolation_y[i]))

        # Отрисовка обоих графиков в одной системе координат
        pylab.plot(x_list, func_value)
        pylab.plot(x_list, poli_value)
        # Открытие окна с графиками
        pylab.show()
    else:
        print('Нет такой программы')


if __name__ == '__main__':
    main()
