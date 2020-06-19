#Численное интегрирования методами средних прямоугольников, трапеций и Симпсона
from math import *
import numpy
import sympy


# Вычисление значения функции по её названию и значению аргумента
def compute_function(func_name, x):
    return eval(func_name)


# Вычисление интеграла всеми методами с фиксированным шагом
def integral_fixed_step(func_name, a, b, h):

    amount_digits = len(str(int(1 / h)))  # Количество знаков после запятой у шага сетки

    # Заполняем список узлов квадратурной формулы
    list_x = []
    x = a
    while x <= b:  # Добавляем в список узлов значения, пока не дойдём до b
        list_x.append(x)
        amount_digits_x = len(str(x).split('.')[1])  # Количество знаков после запятой у x
        x += h
        # Округляем, чтобы избежать ошибки float
        if amount_digits > amount_digits_x:
            x = round(x, amount_digits)
        else:
            x = round(x, amount_digits_x)
    # количество узлов квадратурной формулы
    N = len(list_x)

    # переменные для хранения результата вычисления интеграла
    res_rectangles = 0
    res_trapeze = 0
    res_simpson = 0

    # Вычисление интеграла методом средних квадратов
    # res_rectangles = f(x_0.5) + f(x_1.5) + ... + f(x_n-0.5)
    for i in range(1, N):
        res_rectangles += compute_function(func_name, list_x[i] - h / 2)
    res_rectangles *= h

    # Вычисление интеграла методом трапеций
    # res_trapeze = (f(x0)/2 + f(x1) + f(x2) + ... + f(xn-1) + f(xn) / 2) * h
    res_trapeze += compute_function(func_name, list_x[0]) / 2 + compute_function(func_name, list_x[N - 1]) / 2
    for i in range(1, N - 1):
        res_trapeze += compute_function(func_name, list_x[i])
    res_trapeze *= h

    # Вычисление интеграла методом Симпсона
    # res_simpson = (f(x0) + 4f(x1) + 2f(x2) + ... + 2f(xn-2) + 4f(xn-1) + f(xn)) * (h/3)
    cur_koef = 2    # Текущий коэффициент в выражении
    x = list_x[0]
    res_simpson += compute_function(func_name, list_x[0]) + compute_function(func_name, list_x[N-1])
    for i in range(1, N - 1):
        cur_koef = 2 if cur_koef == 4 else 4  # Заменяем коэффициент при переменной на каждом шаге
        res_simpson += compute_function(func_name, list_x[i]) * cur_koef
    res_simpson *= (h / 3)

    return res_rectangles, res_trapeze, res_simpson


# Вычисление интеграла всеми методами с автоматическим шагом
def integral_auto_step(func_name, a, b, epsilon):
    h = 1  # начальное значение шага
    steps_count = 0  # счетчик числа шагов при вычислении интеграла

    # переменные для хранения результата вычисления интеграла
    res_rectangles = 0
    res_trapeze = 0
    res_simpson = 0

    # Вычисление методом средних прямоугольников
    steps_count = 0
    h = 1
    val_x = a
    while val_x < b:  # Продолжаем вычисления, пока не дойдем до границы промежутка
        cur_val = val_x + h  # текущее значение узла
        # Если вышли с таким шагом за границу промежутка
        if cur_val > b:
            cur_val = b
            h = b - val_x
        # Вычисление интеграла In с шагом h
        in_h = compute_function(func_name, cur_val - h/2) * h
        # Вычисление интеграла In с шагом h/2
        in_h2 = compute_function(func_name, cur_val - h/4) * h/2 + compute_function(func_name, cur_val - 3*h/4) * h/2
        # Если ещё не достигли заданной точности, то продолжаем уменьшать шаг
        if abs(in_h - in_h2) > epsilon:
            h /= 2
        else:
            # Иначе можно перейти к  следующему узлу и увеличить шаг
            steps_count += 1
            res_rectangles += in_h
            val_x += h
            h *= 2
    print("Вычисление методом средних прямоугольников в {0} шагов".format(steps_count))

    # Вычисление методом трапеций
    steps_count = 0
    h = 1
    val_x = a
    while val_x < b:  # Продолжаем вычисления, пока не дойдем до границы промежутка
        if val_x + h > b:
            h = b - val_x
        # Вычисление интеграла с шагом h
        in_h = h * (compute_function(func_name, val_x) + compute_function(func_name, val_x + h)) / 2
        # Вычисление интеграла с шагом h/2
        in_h2 = h/4 * (compute_function(func_name, val_x) + 2 * compute_function(func_name, val_x + h/2)
                       + compute_function(func_name, val_x + h))
        # Если ещё не достигли заданной точности, то продолжаем уменьшать шаг
        if abs(in_h - in_h2) > epsilon:
            h /= 2
        else:
            steps_count += 1
            res_trapeze += in_h
            val_x += h
            h *= 2
    print("Вычисление методом трапеций в {0} шагов".format(steps_count))

    # Вычисление формулой Симпсона
    steps_count = 0
    h = 1
    val_x = a
    while val_x < b:  # Продолжаем вычисления, пока не дойдем до границы промежутка
        if val_x + h > b:
            h = b - val_x
        # Вычисление интеграла с шагом h
        in_h = h / 6 * (compute_function(func_name, val_x) + 4 * compute_function(func_name, val_x + h/2)
                        + compute_function(func_name, val_x + h))
        # Вычисление интеграла с шагом h/2
        in_h2 = (compute_function(func_name, val_x) + 4 * compute_function(func_name, val_x + h/4)
                 + 2 * compute_function(func_name, val_x + h/2) + 4 * compute_function(func_name, val_x + 3*h/4) +
                 compute_function(func_name, val_x + h)) * h / 12
        if abs(in_h - in_h2) > epsilon:
            h /= 2
        else:
            steps_count += 1
            res_simpson += in_h
            val_x += h
            h *= 2
    print("Вычисление методом Симпсона в {0} шагов".format(steps_count))
    return res_rectangles, res_trapeze, res_simpson


def main():
    function_name = ""      # наименование функции
    a = 0                   # нижний предел интегрирования
    b = 0                   # верхний предел интегрирования
    h = 0                   # шаг сетки
    eps = 0                   # точность вычислений при автоматическом выборе шага

    # чтение данных из файла
    with open("input.txt", 'r') as inp:
        function_name = inp.readline().replace('\n', '')
        a, b = map(float, (inp.readline()).split())
        h = float(inp.readline())
        eps = float(inp.readline())

    # вывод данных на экран
    print("Вычисление определенного интеграла от {0} на [{1}, {2}] ".format(function_name, a, b))
    x = sympy.Symbol('x')
    print("Точное значение = {0}".format(sympy.integrate(function_name, (x, a, b))))
    print('_' * 100)
    print("С фиксированным шагом h = {0}".format(h))
    # сам процесс вычисления интеграла с фиксированным шагом
    fix_rectangle, fix_trapeze, fix_simpson = integral_fixed_step(function_name, a, b, h)
    print("Вычисление методом средних прямоугольников: {0}".format(fix_rectangle))
    print("Вычисление методом трапеций:                {0}".format(fix_trapeze))
    print("Вычисление методом Симпсона:                {0}".format(fix_simpson))

    print('#' * 100)

    print("С автоматическим выбором шага h, с точностью epsilon = {0}".format(eps))
    # процесс вычисления интеграла с автоматическим выбором шага
    auto_rectangle, auto_trapeze, auto_simpson = integral_auto_step(function_name, a, b, eps)
    print("\nВычисление методом средних прямоугольников: {0}".format(auto_rectangle))
    print("Вычисление методом трапеций:                {0}".format(auto_trapeze))
    print("Вычисление методом Симпсона:                {0}".format(auto_simpson))


if __name__ == '__main__':
    main()
