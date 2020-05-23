# Численное дифференцирование (поиск 1, 2 и 3 производных на отрезке)
from math import *
import numpy
import sympy
import random


# Вычисление 1, 2 и 3 производных с помощью интерполяционного многочлена Лагранжа
def find_derivations(x, u):
    derivative1_list = []  # Список из первых производных исходной функции
    derivative2_list = []  # Список из вторых производных исходной функции
    derivative3_list = []  # Список из третих производных исходной функции
    h = abs(x[1] - x[0])  # Шаг сетки аргументов функции
    size = len(x)

    # ПЕРВЫЕ ПРОИЗВОДНЫЕ #
    # Нахождение производной в x0 (из L2)
    derivative1_list.append((-3 * u[0] + 4 * u[1] - u[2]) / (2 * h))
    # Нахождение производной в xi (i = 1,.., n-1)
    for i in range(1, size - 1):
        derivative1_list.append((u[i + 1] - u[i - 1]) / (2 * h))
    # Нахождение производной в xn (формула была выведена из L2)
    derivative1_list.append((u[size - 3] - 4 * u[size - 2] + 3 * u[size - 1]) / (2 * h))

    # ВТОРЫЕ ПРОИЗВОДНЫЕ #
    # Нахождение 2 производной в x0 (из L3)
    derivative2_list.append((u[0] * 2 - u[1] * 5 + u[2] * 4 - u[3]) / (h ** 2))
    # Нахождение 2 производной в xi (i = 1,.., n-1)
    for i in range(1, len(x) - 1):
        derivative2_list.append((u[i - 1] - 2 * u[i] + u[i + 1]) / (h ** 2))
    # Нахождение 2 производной в xn (из L3)
    derivative2_list.append((u[size - 1] * 2 - u[size - 2] * 5 + u[size - 3] * 4 - u[size - 4]) / (h ** 2))

    # ТРЕТЬИ ПРОИЗВОДНЫЕ #z
    # Нахождение 3 производной в x0 и x1 (из L4)
    derivative3_list.append((-5 * u[0] + 18 * u[1] - 24 * u[2] + 14 * u[3] - 3 * u[4]) / (2 * (h ** 3)))
    derivative3_list.append((-3 * u[0] + 10 * u[1] - 12 * u[2] + 6 * u[3] - u[4]) / (2 * (h ** 3)))
    # Нахождение 3 производной в xi (i = 2,.., n-2)
    for i in range(2, len(x) - 2):
        derivative3_list.append((-u[i - 2] + 2 * u[i - 1] - 2 * u[i + 1] + u[i + 2]) / (2 * (h ** 3)))
    # Нахождение 3 производной в xn-1 и xn (из L4)
    derivative3_list.\
        append((3 * u[size - 1] - 10 * u[size - 2] + 12 * u[size - 3] - 6 * u[size - 4] + u[size - 5])
               / (2 * (h ** 3)))
    derivative3_list.\
        append((5 * u[size - 1] - 18 * u[size - 2] + 24 * u[size - 3] - 14 * u[size - 4] + 3 * u[size - 5])
               / (2 * (h ** 3)))

    return derivative1_list, derivative2_list, derivative3_list


def main():
    random.seed()

    function_name = ""   # Имя функции
    function_bord = []   # Границы функции
    argument_value = []  # Значения аргументов
    function_value = []  # Значение функции в этих точках
    step = 0             # Шаг сетки узлов интерполирования
    interference = 0     # Предел возмущений

    with open("input.txt", "r") as inp:
        function_name = inp.readline().replace('\n', '')
        function_bord = list(map(float, (inp.readline()).split()))
        step = float(inp.readline())
        interference = float(inp.readline())
    amount_digits = len(str(int(1 / step)))  # Количество знаков после запятой у шага сетки
    # Заполнение таблицы значений аргументов и функции
    x = function_bord[0]  # Изначально x - самое минимальное значение
    while x <= function_bord[1]:  # Добавляем в таблицу значения, пока не дойдём до x max
        argument_value.append(x)
        function_value.append(eval(function_name))
        amount_digits_x = len(str(x).split('.')[1])  # Количество знаков после запятой у x
        x += step
        # Округляем, чтобы избежать ошибки float
        if amount_digits > amount_digits_x:
            x = round(x, amount_digits)
        else:
            x = round(x, amount_digits_x)

    # Вычисление первой, второй и третьей производной
    first_derivation, second_derivation, third_derivation = find_derivations(argument_value, function_value)

    print('Таблица значений функции и производных для функции {0} на отрезке  [{1}; {2}]:'
          .format(function_name, function_bord[0], function_bord[1]))
    print('(Шаг сетки равен {0})'.format(step))

    # Вычисление точных производных
    x = sympy.Symbol('x')
    df1 = sympy.diff(function_name, x, 1)
    df2 = sympy.diff(function_name, x, 2)
    df3 = sympy.diff(function_name, x, 3)

    print("Точные значения производной:")
    print("\t\t\tx\t\t\t\t\tf(x)\t\t\t\tf'(x)\t\t\t\tf''(x)\t\t\tf'''(x)")
    for i in range(len(argument_value)):
        # Вычисление точных значений производных в точках
        x = argument_value[i]
        df1_x = eval(str(df1))
        df2_x = eval(str(df2))
        df3_x = eval(str(df3))
        print('{0})\t{1:^20}{2:^20.4f}{3:^20.4f}{4:^20.4f}{5:^20.4f}'
              .format(i+1, argument_value[i], function_value[i],
                      df1_x, df2_x, df3_x))
    print("\nВычисленные значения производной:")
    print("\t\t\tx\t\t\t\t\tf(x)\t\t\t\tf'(x)\t\t\t\tf''(x)\t\t\tf'''(x)")
    for i in range(len(argument_value)):
        print('{0})\t{1:^20}{2:^20.4f}{3:^20.4f}{4:^20.4f}{5:^20.4f}'
              .format(i+1, argument_value[i], function_value[i],
                      first_derivation[i], second_derivation[i], third_derivation[i]))

    # Добавляем возмущения к значению исходной функции
    inter_function_value = []
    for i in range(len(function_value)):
        inter_function_value.append(function_value[i] + random.uniform(-interference, interference))
    # Вычисляем производные для таких новых значений
    inter_first_derivation, inter_second_derivation, inter_third_derivation \
        = find_derivations(argument_value, inter_function_value)

    print("\nПроизводные функции с возмущениями (разница в скобках):")
    print("\t\t\tx\t\t\t\t\tf(x)\t\t\t\t\t\tf'(x)\t\t\t\t\t\tf''(x)\t\t\t\t\tf'''(x)")
    for i in range(len(argument_value)):
        print('{0})\t{1:^15}{2:^15.4f}({3:<8.5f})\t{4:^15.4f}({5:<8.4f})\t{6:^15.4f}({7:<8.4f})\t{8:^15.4f}({9:<8.4f})'
              .format(i + 1, argument_value[i], inter_function_value[i], function_value[i] - inter_function_value[i],
                      inter_first_derivation[i], first_derivation[i] - inter_first_derivation[i],
                      inter_second_derivation[i], second_derivation[i] - inter_second_derivation[i],
                      inter_third_derivation[i], third_derivation[i] - inter_third_derivation[i]))


if __name__ == '__main__':
    main()
