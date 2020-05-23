# Решение задачи Коши для системы обыкновенных Дифференциальных уравнений
from math import *
import pylab
import sympy
from copy import deepcopy


start_system = []       # Исходная система
var_names = []          # Имена функций
par_name = []           # Имя параметра
true_solution = []      # Решение задачи


# Из этой функции берутся начальные данные
def form_start_system():
    start_system.append("cos(x) * x ** 2 + 2 * sin(x) * x")

    # Имена функций
    var_names.append("y")

    # Название аргумента функций
    par_name.append("x")


# Вычисляет значение функции нескольких переменных
# с аргументом arg_value и остальными функциями fun_values
def func_val(function_name, arg_value, fun_values):
    size = len(fun_values)
    function_name = function_name.replace(str(par_name[0]), "(" + str(arg_value) + ")")
    for i in range(size):   # Замена наименований переменных на их значения
        function_name = function_name.replace(str(var_names[i]), "(" + str(fun_values[i]) + ")")
    res = eval(str(function_name))
    return res


# Вычисление y_n+1
def get_kq(amount_equalities, cur_x, cur_y, step):
    k1 = []
    k2 = []
    k3 = []
    k4 = []
    # Новые значения y_n+1
    new_y = []

    # Вычисление коэффициентов
    for i in range(amount_equalities):
        k1.append(step * func_val(start_system[i], cur_x, cur_y))
    for i in range(amount_equalities):
        k2.append(step * func_val(start_system[i], cur_x + step/2, [cur_y[j] + k1[j]/2 for j in range(amount_equalities)]))
    for i in range(amount_equalities):
        k3.append(step * func_val(start_system[i], cur_x + step/2, [cur_y[j] + k2[j]/2 for j in range(amount_equalities)]))
    for i in range(amount_equalities):
        k4.append(step * func_val(start_system[i], cur_x + step, [cur_y[j] + k3[j] for j in range(amount_equalities)]))

    # Получение значений функций в следующей точке
    new_y = [cur_y[i] + 1/6 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) for i in range(amount_equalities)]
    return new_y


# решение системы методом Рунге-Кутты с заданным шагом
def fixed_runge_kutta(start_cond, left_b, right_b, h):
    amount_equalities = len(start_system)

    x_list = []     # Список значений параметра функции
    y_list = []     # Список решений системы ДУ

    # Заносим в списки начальное решение
    cur_x = right_b
    x_list.append(cur_x)
    y_list.append([start_cond[i] for i in range(len(start_cond))])

    # До окончания интегрирования по всему промежутку
    h = -h
    while cur_x > left_b:
        # Получение решений в следующей точке
        new_y = get_kq(amount_equalities, cur_x, y_list[-1], h)
        cur_x += h

        # Добавление в список решений
        x_list.append(cur_x)
        y_list.append(new_y)
    return x_list, y_list


# решение системы методом Рунге-Кутты с автоматическим шагом
def auto_runge_kutta(start_cond, left_b, right_b, epsilon):
    # Начальная длина шага
    h = abs(right_b - left_b) / 100
    if h < 10e-15:
        h = abs(right_b - left_b) / 10

    # Граничные длины для шага
    h_max = abs(right_b - left_b) / 5
    h_min = abs(right_b - left_b) / 20000

    amount_equalities = len(start_system)
    epsilon1 = epsilon / 32

    x_list = []     # Список значений параметра функции
    y_list = []     # Список решений системы ДУ

    is_less_prev = False    # уменьшение шага на предыдущей итерации
    is_less = False         # уменьшение шага на текущей итерации

    # Заносим в списки начальное решение
    cur_x = left_b
    x_list.append(cur_x)
    y_list.append([start_cond[i] for i in range(amount_equalities)])
    cur_y = [start_cond[i] for i in range(amount_equalities)]

    k1 = []
    k2 = []
    k3 = []
    k4 = []
    k5 = []
    k6 = []

    # До окончания интегрирования по всему промежутку
    while cur_x < right_b:
        is_less = False
        while True:
            # считаем ошибку по формулам Фельберга
            for i in range(amount_equalities):
                k1.append(h * func_val(start_system[i], cur_x, cur_y))
            for i in range(amount_equalities):
                k2.append(h * func_val(start_system[i], cur_x + h / 4,
                                          [cur_y[j] + (k1[j] / 4) for j in range(amount_equalities)]))
            for i in range(amount_equalities):
                k3.append(h * func_val(start_system[i], cur_x + (3*h)/8, [cur_y[j] - (3 * k1[j] / 32) + (9 * k2[j] / 32)
                                                                          for j in range(amount_equalities)]))
            for i in range(amount_equalities):
                k4.append(h * func_val(start_system[i], cur_x + 12 * h / 13,
                                       [cur_y[j] + 1932 / 2197 * k1[j] - 7200 / 2197 * k2[j] + 7296/2197 * k3[j]
                                        for j in range(amount_equalities)]))

            for i in range(amount_equalities):
                k5.append(h * func_val(start_system[i], cur_x + h,
                                       [cur_y[j] + 439 / 216 * k1[j] - 8 * k2[j] + 3680/513 * k3[j] - 845/4104 * k4[j]
                                        for j in range(amount_equalities)]))

            for i in range(amount_equalities):
                k6.append(h * func_val(start_system[i], cur_x + h/2,
                                       [cur_y[j] - 8/27 * k1[j] + 2 * k2[j] - 3544/2565 * k3[j] + 1859/4104 * k4[j] - 11/40 * k5[j]
                                        for j in range(amount_equalities)]))

            # Вычисление ошибки
            err_eq = []
            for i in range(amount_equalities):
                err_eq.append(k1[i] / 360 - 128/4275 * k3[i] - 2197/75240 * k4[i] + k5[i] / 50 + 2/55 * k6[i])
            err = abs(max(err_eq) / 31)

            # Очистка списков с коэффициентами
            k1.clear()
            k2.clear()
            k3.clear()
            k4.clear()
            k5.clear()
            k6.clear()

            # Если ошибка меньше epsilon или длина шага больше максимальной(меньше минимальной)
            if err <= epsilon or h < h_min or h > h_max:
                if h < h_min:
                    h = h_min
                if h > h_max:
                    h = h_max
                break
            else:
                # иначе продолжаем уменьшать шаг
                h /= 2
                is_less = True

        # Получаем решение на данной  итерации с заданным значением шага
        new_y = get_kq(amount_equalities, cur_x, cur_y, h)
        cur_x += h
        cur_y = deepcopy(new_y)
        # Добавление в список новых решений
        x_list.append(cur_x)
        y_list.append(new_y)

        if not is_less and is_less_prev:  # если на прошлой итерации h < eps1
            h *= 2

        # Если ошибка меньше, чем eps1 ( 0 < eps1 < epsilon )
        if err < epsilon1:
            is_less_prev = True
        else:
            is_less_prev = False
    return x_list, y_list


def main():
    form_start_system()                     # формирование начальных условий
    amount_functions = len(start_system)    # количество функций в системе
    start_conditions = []                   # начальные условия задачи
    h = 0                                   # шаг
    a = b = 0                               # границы интегрирования
    epsilon = 1                             # точность

    print("Система имеет вид:")
    for i in range(amount_functions):
        print("{0}'({1}) = {2}".format(var_names[i], par_name[0], start_system[i]))

    print("Ввод начальных условий:")
    for i in range(amount_functions):
        start_conditions.append(float(input("{0}({1}0) = ".format(var_names[i], par_name[0]))))

    print("Промежуток интегрирования:", end="")
    a, b = map(float, input().split())

    # Если шаг None, то проводим вычисления с автоматическим шагом и заданной точностью
    # иначе ищем решения с шагом h
    print("Шаг:", end="")
    in_str = input()
    if in_str == "None":
        h = -1
        print("Точность:", end="")
        epsilon = float(input())
    else:
        h = float(in_str)

    # решение системы методом Рунге-Кутты
    x_values = []
    y_values = []
    if h == -1:  # с автоматическим шагом
        print("Решение задачи с автоматическиим выбором шага и точностью eps = {0}\n".format(epsilon))
        x_values, y_values = auto_runge_kutta(start_conditions, a, b, epsilon)
    else:  # с заданным шагом
        print("Решение задачи с фиксированным шагом h = {0}\n".format(h))
        x_values, y_values = fixed_runge_kutta(start_conditions, a, b, h)

    # Таблица значений функций
    out = open("output.txt", "w")   # запись в файл
    amount_points = len(x_values)
    cur_str = "%15s" % par_name[0]
    for i in range(amount_functions):
        cur_str += "%15s" % var_names[i]
    print(cur_str)
    out.write(cur_str + "\n")
    for i in range(amount_points):
        cur_str = str("%15f" % x_values[i])
        for j in range(amount_functions):
            cur_str += str("%15f" % y_values[i][j])
        print(cur_str)
        out.write(cur_str + "\n")
    out.close()

    # Рисование графика
    fig, axes = pylab.subplots()
    for i in range(amount_functions):
        axes.plot(x_values, [y_val[i] for y_val in y_values], label="{0}({1})".format(var_names[i], par_name[0]))
    axes.legend(loc='lower right')

    # Задано точное решение системы, тогда выводится его графики
    if len(true_solution) > 0:
        true_y = []
        for i in range(amount_points):
            true_y.append([func_val(true_solution[j], x_values[i], []) for j in range(amount_functions)])

        # отрисовка графика
        fig2, setka = pylab.subplots()
        for i in range(amount_functions):
            setka.plot(x_values, [y[i] for y in true_y],
                       label="{0}({1}) = {2}".format(var_names[i], par_name[0], true_solution[i]))
        setka.legend(loc='lower right')

        # погрешность
        norm = 0
        for i in range(amount_points):
            for j in range(amount_functions):
                razn = abs(true_y[i][j] - y_values[i][j])
                if razn > norm:
                    norm = razn
        print("Норма глобальной погрешности = {0}".format(norm))
    pylab.show()


if __name__ == "__main__":
    main()
