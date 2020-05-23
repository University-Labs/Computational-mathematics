# Решение систем нелинейных уравнений методом Ньютона
from math import *
from functools import reduce
from copy import copy

start_system = []       # Исходная система
differ_system = []      # Первые производные указанных функций
var_names = []          # имена переменных


# Из этой функции берутся начальные данные
def form_start_system():
    start_system.append("x ** 2 + (y ** 2) / 4 - 1")
    start_system.append("x ** 2 + (y + 0.5) ** 2 - 10")

    differ_system.append(["2 * x", "y / 2"])
    differ_system.append(["2 * x", "2 * y + 1"])

    var_names.append("x")
    var_names.append("y")


# преобразование матрицы A в треугольную
def make_triangle(matrix_a, free_b, is_main_elem = True):
    length = len(matrix_a)
    amount_changing = 0  # количество перестановок строк
    for i in range(0, length):
        # Если встретили на главной диагонали нулевой элемент
        if abs(matrix_a[i][i]) < 10e-16:
            for j in range(i+1, length):
                if matrix_a[j][i] != 0:  # Ищем под ним строку j с не нулём на позиции i
                    matrix_a[j], matrix_a[i] = matrix_a[i], matrix_a[j]
                    free_b[j], free_b[i] = free_b[i], free_b[j]
                    amount_changing += 1  # поменяли строки местами
                    break
            if matrix_a[i][i] < 10e-16:  # матрица вырождена
                return -1

        # Если происходит поиск главного элемента
        if is_main_elem:
            max_index = i
            for j in range(i+1, length):
                if abs(matrix_a[max_index][i]) < abs(matrix_a[j][i]):
                    max_index = j
            if max_index != i:
                amount_changing += 1  # выполнили перестановку строки
                matrix_a[max_index], matrix_a[i] = matrix_a[i], matrix_a[max_index]
                free_b[max_index], free_b[i] = free_b[i], free_b[max_index]

        # Отнимаем от каждой строчки i-ую, умноженную на C
        for j in range(i + 1, length):
            c = matrix_a[j][i] / matrix_a[i][i]
            for k in range(i, length):
                matrix_a[j][k] -= matrix_a[i][k] * c
                if abs(matrix_a[j][k]) < 10e-16:  # чтобы избежать ошибок округления
                    matrix_a[j][k] = 0
            free_b[j] -= free_b[i] * c
            if abs(free_b[j]) < 10e-16:  # чтобы избежать ошибок округления
                free_b[j] = 0
    return amount_changing


# нахождение решений x1,...xm (обратный ход Гаусса)
def find_decision(matrix_a, free_b):
    length = len(matrix_a)
    decision = [0 for _ in range(length)]
    for i in range(length - 1, -1, -1):
        answer = free_b[i]
        for j in range(i + 1, length):
            answer -= decision[j] * matrix_a[i][j]
        answer /= matrix_a[i][i]
        if abs(answer) < 10e-16:  # чтобы избежать ошибок округления
            answer = 0
        decision[i] = answer
    return decision


# Получение матрицы Якоби
def yakobi_val(value_x):
    res_matrix = []
    amount_equality = len(differ_system)

    for i in range(amount_equality):
        temp = list()
        for j in range(amount_equality):
            temp.append(func_val(differ_system[i][j], value_x))
        res_matrix.append(temp)

    return res_matrix


# решение системы методом Ньютона
def newton_method(start_decision, epsilon):
    x_k = start_decision    # начальное приближение
    amount_var = len(x_k)   # количество переменных
    iter_count = 0          # количество итераций

    # Пока не превысим возможное количество итераций (бесконечный цикл),
    # либо пока норма вектора d_x не будет меньше epsilon
    while iter_count < 100000:
        iter_count += 1

        # Вычисление матрицы Якоби для x_k
        matrix_b = yakobi_val(x_k)

        # Переносим F(x_k) вправо
        minus_f = [(-1) * func_val(start_system[i], x_k) for i in range(amount_var)]

        # Решаем СЛАУ B*dx = -F методом Гаусса
        is_correct = make_triangle(matrix_b, minus_f)
        if is_correct == -1:        # матрица вырождена
            print("Матрица вырождена")
            iter_count = 100001
            break
        d_x = find_decision(matrix_b, minus_f)

        if any(map(lambda xi: xi >= 10e15, d_x)):
            print("Метод не сходится")
            iter_count = 100001
            break

        # Фомирование x_(k+1)
        x_k = [d_x[i] + x_k[i] for i in range(amount_var)]

        # Вычисляем норму вектора
        if sqrt((float(reduce(lambda x, y: x ** 2 + y ** 2, d_x)))) < epsilon:
            break
        # Иначе продолжаем поиски решения
    return x_k, iter_count


# Вычисляет значение функции нескольких переменных с заданным значением переменных
def func_val(function_name, values):
    size = len(values)
    for i in range(size):   # Замена наименований переменных на их значения
        function_name = function_name.replace(str(var_names[i]), "(" + str(values[i]) + ")")
    res = eval(str(function_name))
    return res


def main():
    start_decision = []     # Начальное приближение
    epsilon = 1             # Точность

    print("Решить систему:")
    form_start_system()         # Формирование начальной системы программно
    for function_name in start_system:
        print("{0:20} = 0;".format(function_name))
    print('Матрица Якоби имеет вид: ')
    amount_var = len(differ_system[0])
    for func_i in differ_system:
        for i in range(amount_var):
            print("{0:20}".format(func_i[i]), end="\t")
        print()

    # Ввод начального приближения
    print("Начальное приближение:")
    start_decision = list(map(float, input().split()))

    # Ввод точности
    print("Точность:")
    epsilon = float(input())

    # Решение системы нелинейных уравнений методом Ньютона
    result, iter_count = newton_method(start_decision, epsilon)
    if iter_count >= 100000:
        print("Решение не было найдено")
    else:
        print("Решение:")
        for i in range(len(result)):
            print("{0} = {1}".format(var_names[i], result[i]))
        print("Решение было найдено за {0} итераций".format(iter_count))
        print("Невязка")
        for i in range(len(start_system)):
            print("r{0} = {1};".format(i, func_val(start_system[i], result)))


if __name__ == '__main__':
    main()
