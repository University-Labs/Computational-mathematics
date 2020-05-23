# Численное интегрирование методом Монте-Карло
import os
from math import *
import pylab
import random


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


# Вычисление первой производной
def find_differ(splines, value_x):
    # определяем какому много члену принадлежит значение value_x
    size = len(splines)
    if value_x <= splines[0][0]:  # значение x слева от минимального сплайна
        curr_spline = splines[0]
    elif value_x >= splines[size - 1][0]:  # значение x справа от максимального сплайна
        curr_spline = splines[size - 1]
    else:  # значение x в одном из внутренних Si
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

    # Вычисляем значение первой производной в заданной точке
    return curr_spline[2] + curr_spline[3] * (value_x - curr_spline[0]) \
           + 0.5 * curr_spline[4] * (value_x - curr_spline[0]) ** 2


# Вычисление интеграла по формуле Симпсона
def simpson(x_spline, y_spline, a, b, epsilon):
    # Вычисление формулой Симпсона
    res_simpson = 0
    h = 1
    cur_t = a
    while cur_t < b:  # Продолжаем вычисления, пока не дойдем до границы промежутка
        if cur_t + h > b:
            h = b - cur_t
        # Вычисление интеграла с шагом h
        in_h = h / 6 * (get_val(cur_t, x_spline, y_spline) + 4 * get_val(cur_t + h / 2, x_spline, y_spline)
                        + get_val(cur_t + h, x_spline, y_spline))
        # Вычисление интеграла с шагом h/2
        in_h2 = (get_val(cur_t, x_spline, y_spline) + 4 * get_val(cur_t + h / 4, x_spline, y_spline)
                 + 2 * get_val(cur_t + h / 2, x_spline, y_spline) + 4 * get_val(cur_t + 3 * h / 4, x_spline, y_spline) +
                 get_val(cur_t + h, x_spline, y_spline)) * h / 12
        if abs(in_h - in_h2) > epsilon:
            h /= 2
        else:
            res_simpson += in_h
            cur_t += h
            h *= 2
    return abs(res_simpson)


# Вычисление площади области, ограниченной кривой
def get_val(t, x_spline, y_spline):
    return spline_values(y_spline, t) * find_differ(x_spline, t)


# Вычисление площади методом Монте-Карло
def monte_carlo(spline_values_x, spline_values_y, amount_try):
    x_min = min(spline_values_x)
    x_max = max(spline_values_x)
    y_min = min(spline_values_y)
    y_max = max(spline_values_y)

    random_points = []  # Добавленные точки
    amount = 0     # Количество точек внутри области

    # Добавление точек
    for i in range(amount_try):
        # Генерация случайной точки
        x = random.uniform(x_min, x_max)
        y = random.uniform(y_min, y_max)
        # Внутри ли области заданная точка
        if is_inside(x, y, spline_values_x, spline_values_y):
            amount += 1
            random_points.append((x, y, True))
        else:
            random_points.append((x, y, False))

    # S = J * 1/amount_try * доля точек внутри области
    j = (abs((x_max - x_min) * (y_max - y_min)))
    s = j * amount / amount_try
    return s, random_points


# Проверка, принадлежит ли точка заданной области
def is_inside(x, y, x_points, y_points):
    result = False
    size = len(x_points)
    j = size - 1
    for i in range(size):
        if ((y_points[i] < y and y_points[j] >= y) or (y_points[j] < y and y_points[i] >= y)) and \
                (x_points[i] + (y - y_points[i]) / (y_points[j] - y_points[i]) * (x_points[j] - x_points[i])) < x:
            result = False if result is True else True
        j = i
    return result


def main():
    x_values = []   # список значений x
    y_values = []   # список значений y
    value_t = []    # список значений параметра t
    eps = 1         # точность вычисления
    amount_try = 1  # количество вычислений интеграла
    integral_simpson = 0
    integral_monte_karlo = 0

    # Чтение данных из файла
    with open("input.txt", "r") as inp:
        next_in = False
        for line in inp.readlines():
            list_m = line.split()
            if len(list_m) == 2:
                x_values.append(float(list_m[0]))
                y_values.append(float(list_m[1]))
            else:
                if next_in is False:
                    eps = float(list_m[0])
                    next_in = True
                else:
                    amount_try = int(list_m[0])
                    break

    # Задаем параметр t
    value_t = [i for i in range(len(x_values))]

    # формирование списка значений параметра для красивого вывода
    t_list = []
    tmp_val = 0
    max_v = value_t[len(value_t) - 1]
    while tmp_val <= max_v:
        t_list.append(tmp_val)
        tmp_val += 0.01

    # Вывод  координат точек
    print("Таблица значений кривой:")
    print('\t\tti\t\t\tXi\t\t\tYi')
    for i in range(len(x_values)):
        print('{0})\t{1:^10.4f}\t{2:^10.4f}\t{3:^10.4f}'
              .format(i + 1, value_t[i], x_values[i], y_values[i]))

    # Вычисление сплайна для x = x(t)
    splines_x = spline_interpolation(value_t, x_values)
    spline_answers_x = [spline_values(splines_x, x) for x in t_list]

    # Вычисление сплайна для y = y(t)
    splines_y = spline_interpolation(value_t, y_values)
    spline_answers_y = [spline_values(splines_y, y) for y in t_list]

    # Вычисление площади по формуле Симпсона
    integral_simpson = simpson(splines_x, splines_y, value_t[0], value_t[len(value_t) - 1], eps)
    print(100 * '_')
    print("По формуле Cимпсона: {0}".format(integral_simpson))

    # Площадь методом Монте-Карло
    integral_monte_karlo, points = monte_carlo(spline_answers_x, spline_answers_y, amount_try)
    print("По формуле Монте-Карло: {0}".format(integral_monte_karlo))
    print("Количество испытаний: {0}".format(amount_try))

    # отрисовка генерируемых точек (внутри - зеленым, снаружи - красным)
    for i in range(len(points)):
        x, y, inside = points[i]
        if i < 1000:
            if inside is True:
                pylab.scatter(x, y, s=5, c="green")
            else:
                pylab.scatter(x, y, s=5, c="red")

    # Рисование графика
    pylab.plot(spline_answers_x, spline_answers_y, "c")
    pylab.grid(True)
    pylab.show()


if __name__ == '__main__':
    main()
