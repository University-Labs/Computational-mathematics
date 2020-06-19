# Метод стрельбы решения краевой задачи
from math import *
import pylab
import sympy
from copy import deepcopy


start_system = []       # Исходная система
var_names = []          # Имена функций
par_name = []           # Имя параметра


# Из этой функции берутся начальные данные
def form_start_system(betta, gamma, delta):
    start_system.append("z")
    start_system.append("delta * y * exp(gamma * betta * (1-y) / (1 + betta * (1-y)))")
    start_system[1] = start_system[1].replace(str("betta"), "(" + str(betta) + ")")
    start_system[1] = start_system[1].replace(str("gamma"), "(" + str(gamma) + ")")
    start_system[1] = start_system[1].replace(str("delta"), "(" + str(delta) + ")")

    # Имена функций
    var_names.append("y")
    var_names.append("z")

    # Название аргумента функций
    par_name.append("t")


# Вычисляет значение функции нескольких переменных
# с аргументом arg_value и остальными функциями fun_values
def func_val(function_name, arg_value, fun_values):
    size = len(fun_values)
    function_name = function_name.replace(str(par_name[0]), "(" + str(arg_value) + ")")
    for i in range(size):   # Замена наименований переменных на их значения
        function_name = function_name.replace(str(var_names[i]), "(" + str(fun_values[i]) + ")")
    res = eval(str(function_name))
    return res


# решение системы методом Рунге-Кутты с автоматическим шагом
def auto_runge_kutta(start_cond, left_b, right_b, epsilon):
    direct = True   # направление интегрирования

    # Начальная длина шага
    h = abs(right_b - left_b) / 100
    if h < 10e-15:
        h = abs(right_b - left_b) / 10

    # необходимо сменить направление интегрирования, в случае если интегрируем от большего числа к меньшему
    if left_b > right_b:
        h = -h
        direct = False

    # Граничные длины для шага
    h_max = abs(right_b - left_b) / 5
    h_min = abs(right_b - left_b) / 20000

    amount_equalities = len(start_system)
    epsilon1 = epsilon / 64

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
    while fin(cur_x, right_b, direct):
        is_less = False
        while True:
            # считаем ошибку
            for i in range(amount_equalities):
                k1.append(h * func_val(start_system[i], cur_x, cur_y))
            for i in range(amount_equalities):
                k2.append(h * func_val(start_system[i], cur_x + h / 4,
                                          [cur_y[j] + 0.25 * k1[j] for j in range(amount_equalities)]))
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
            if err <= epsilon or abs(h) < h_min or abs(h) > h_max:
                if abs(h) < h_min:
                    h = (abs(h) / h) * h_min
                if abs(h) > h_max:
                    h = (abs(h) / h) * h_max
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

    # если вычислили за границей отрезка, то заменяем это значение, чтобы высчитать на его границе
    if direct is True:
        border = right_b
    else:
        border = left_b
    if x_list[-1] > border:
        cur_x -= h
        h = border - x_list[-2]
        new_y = get_kq(amount_equalities, cur_x, cur_y, h)
        x_list.pop()
        y_list.pop()
        x_list.append(left_b)
        y_list.append(new_y)
    return x_list, y_list


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


# определяет окончание условия цикла
def fin(cur_x, right_b, direct):
    if direct is True:
        return cur_x < right_b
    else:
        return cur_x > right_b


# решение краевой задачи методом стрельбы
def decide_task(start_cond, fin_cond, left_b, right_b, epsilon):
    # пусть изначально z(1) = shooting_range
    shooting_range = 0.5
    x_values_r, y_values_r = [], []             # Содержит решение задачи коши при параметре на левой границе отрезка
    x_values_l, y_values_l = [], []             # Содержит решение задачи коши при параметре на правой границе отрезка
    x_values_middle, y_values_middle = [], []   # Содержит решение задачи коши при параметре на середине отрезка

    # будем искать промежуток [-shooting_condition; 0] или [0; +shooting_condition],
    # на котором располагается решение уравнения на параметр shooting_condition
    while shooting_range < 1000:
        # вычисляем решения задач Коши при разных значениях параметра
        x_values_r, y_values_r = auto_runge_kutta([start_cond, shooting_range], left_b, right_b, epsilon)
        x_values_middle, y_values_middle = auto_runge_kutta([start_cond, 0], left_b, right_b, epsilon)
        x_values_l, y_values_l = auto_runge_kutta([start_cond, -shooting_range], left_b, right_b, epsilon)

        # Если на разных концах промежутков значение z(1) лежит по разные стороны от fin_cond,
        # то где-то на них присутствует решение уравнения
        if (y_values_r[-1][1] - fin_cond) * (y_values_middle[-1][1] - fin_cond) <= 0 or\
                (y_values_l[-1][1] - fin_cond) * (y_values_middle[-1][1] - fin_cond) <= 0:
            break
        else:
            # иначе увеличиваем промежуток просмотра
            shooting_range *= 2

    # получили промежуток, на котором лежит решение
    x_values_old, y_values_old = [], []
    x_values_new, y_values_new = [], []
    if y_values_l[-1][1] * y_values_middle[-1][1] <= 0:
        x_values, y_values = x_values_l, y_values_l
        shooting_range *= -1
    elif y_values_r[-1][1] * y_values_middle[-1][1] <= 0:
        x_values, y_values = x_values_r, y_values_r
    else:
        print("Слишком большой промежуток для поиска")
        return [], []

    # Решается уравнение методом секущих
    old_par, new_par = 0, shooting_range
    x_values_new, y_values_new = auto_runge_kutta([start_cond, old_par], left_b, right_b, epsilon)
    while True:
        # решаем задачу Коши для разных значений параметра
        x_values_old, y_values_old = deepcopy(x_values_new), deepcopy(y_values_new)
        x_values_new, y_values_new = auto_runge_kutta([start_cond, new_par], left_b, right_b, epsilon)

        # берём значение функции в точках (z(0) = y'(0)) с разными значениями z(1)
        f0 = y_values_old[-1][1]
        f1 = y_values_new[-1][1]
        new_par, old_par = new_par - f1 / (f1 - f0) * (new_par - old_par), new_par
        # Если разница между вычисленным параметром и начальным условией меньше заданного эпсилон
        # то параметр new_par найден
        if abs(f1 - fin_cond) < epsilon:
            break

    # возвращаем решение задачи Коши
    return auto_runge_kutta([start_cond, new_par], left_b, right_b, epsilon)


def show_task_cond(start_cond, fin_cond, left_b, right_b, epsilon):
    print("Решить следующую краевую задачу:")
    print("Найти функцию y, где:")
    print("y'' = {0}; \t\ty'({1}) = {2}, y({3}) = {4}".format(start_system[1], right_b, fin_cond, left_b, start_cond))
    print("если точность epsilon = {0}".format(epsilon))


def main():
    # формирование начальной системы условий
    betta = float(input("Введите betta: "))
    gamma = float(input("Введите gamma: "))
    delta = float(input("Введите delta: "))
    form_start_system(betta, gamma, delta)

    start_cond = 1      # y(1) = 1
    fin_cond = 0        # y'(0) = 0
    a = 1               # левая граница
    b = 0               # правая граница
    epsilon = float(input("Введите точность: "))

    # печать условий задачи
    show_task_cond(start_cond, fin_cond, a, b, epsilon)

    # решение поставленной задачи
    x_values, y_values = decide_task(start_cond, fin_cond, a, b, epsilon)

    # Рисование графика
    fig, axes = pylab.subplots()
    for i in range(len(start_system)):
        axes.plot(x_values, [y_val[i] for y_val in y_values], label="{0}({1})".format(var_names[i], par_name[0]))
    axes.legend(loc='lower right')

    # Таблица значений функций
    amount_points = len(x_values)
    cur_str = "%15s" % par_name[0]
    for i in range(len(start_system)):
        cur_str += "%15s" % var_names[i]
    print(cur_str)
    for i in range(amount_points):
        cur_str = str("%15f" % x_values[i])
        for j in range(len(start_system)):
            cur_str += str("%15f" % y_values[i][j])
        print(cur_str)

    pylab.show()


if __name__ == "__main__":
    main()
