# Решение одного нелинейного алгебраического уравнения комбинированным методом хорд и касательных
from math import *
import numpy
import sympy


# Вычисление значения функции по её названию и значению аргумента
def comp_function(func_name, x):
    return eval(func_name)


def main():
    funct_name = ""      # Исходная функция
    first_der = ""      # Первая производная
    second_der = ""     # Вторая производная
    x_list = []         # Список решений уравнения, также хранится число итераций

    # чтение данных из файла
    with open("input.txt", 'r') as inp:
        funct_name = inp.readline().replace('\n', '')
        first_der = inp.readline().replace('\n', '')
        second_der = inp.readline().replace('\n', '')
        a, b = map(float, (inp.readline()).split())
        eps = float(inp.readline())
        h = float(inp.readline())

    # Вывод исходных данных
    print("Исходные данные:\nНайти решения {0} = 0, на промежутке [{1},{2}]".format(funct_name, a, b))
    print("с точностью epsilon = {0}. Шаг = {1}".format(eps, h))
    print("Если известны также её первая f'(x) = {0} и вторая f''(x) = {1} производные".format(first_der, second_der))

    left_border = a         # Левая граница текущего отрезка
    right_border = a + h    # Правая граница текущего отрезка
    # Продолжаем обход промежутка от a до b
    while right_border < b + h:
        if right_border > b:
            right_border = b
        # Если знаки на концах промежутков совпадают, то здесь нет корней
        # (или их четное число и нельзя найти с заданным шагом)
        if comp_function(funct_name, left_border) * comp_function(funct_name, right_border) > eps:
            left_border = right_border
            right_border += h
            continue

        # Иначе на промежутке [left_border, right_border] находится ровно 1 корень
        iter_count = 0
        x1 = left_border
        x2 = right_border
        # Если на границе промежутка имеется нулевое значение, то она является корнем
        if abs(comp_function(funct_name, x1)) < eps:
            x2 = x1
        elif abs(comp_function(funct_name, x2)) < eps:
            x1 = x2
        # Повторяем итерации пока решения не будут достаточно близко к друг другу
        while abs(x1 - x2) > eps:
            iter_count += 1
            x1_tmp = x1
            x2_tmp = x2

            # Если в x1 первая и вторая производные имеют один знак, то уточняем корень в x1
            # методом касательных, а в x2 - методом хорд
            if comp_function(first_der, x1) * comp_function(second_der, x1) > 0:
                x1_tmp = x1 - comp_function(funct_name, x1)/comp_function(first_der, x1)
                x2_tmp = x2 - comp_function(funct_name, x2) * (x2 - x1) /\
                         (comp_function(funct_name, x2) - comp_function(funct_name, x1))
            # Если же в x1 они имеют разные знаки, то уточняем корень в ней
            # методом хорд, а в x2 - методом касательных
            if comp_function(first_der, x1) * comp_function(second_der, x1) < 0:
                x1_tmp = x1 - comp_function(funct_name, x1) * (x2 - x1) / \
                       (comp_function(funct_name, x2) - comp_function(funct_name, x1))
                x2_tmp = x2 - comp_function(funct_name, x2) / comp_function(first_der, x2)
            x1 = x1_tmp
            x2 = x2_tmp
        # Проверим, чтобы не добавлять в список повторяющийся корень
        if len(x_list) > 0:
            if abs(x_list[len(x_list) - 1][0] - ((x1 + x2) / 2)) < eps:
                x_list.pop()
        # Добавляем корень, который находится посередине промежутка
        if abs((x1 + x2) / 2) < eps:
            x_list.append([0, iter_count])
        else:
            x_list.append([(x1 + x2) / 2, iter_count])

        # передвигаемся к следующему отрезку
        left_border = right_border
        right_border += h

    # Вывод корней
    print(100 * '_')
    print("Были найдены следующие корни:")
    if len(x_list) == 0:
        print("Корней нет")
    for x_val in x_list:
        # вычисление невязки
        nevyazka = comp_function(funct_name, x_val[0])
        print("x = {0}, за {1} итераций. Невязка = {2}".format(x_val[0], x_val[1], nevyazka))


if __name__ == '__main__':
    main()
