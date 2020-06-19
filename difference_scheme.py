# Решение краевой задачи с помощью разностной схемы
from math import *
from numpy import linspace
import pylab


def fun_val(func_name, x_value):
    func_name = func_name.replace("X", "(" + str(x_value) + ")")
    return eval(func_name)


# Решение СЛАУ методом прогонки, с заданными коэф-ми a b c d
def slau_solution(a_list, b_list, c_list, d_list):
    length = len(a_list)
    alpha = [0 for _ in range(length + 1)]
    betta = [0 for _ in range(length + 1)]
    answer = [0 for _ in range(length)]

    # прямой ход
    for i in range(length):
        alpha[i + 1] = c_list[i] / (b_list[i] - (alpha[i] * a_list[i]))
        betta[i + 1] = ((betta[i] * a_list[i]) - d_list[i]) / (b_list[i] - (alpha[i] * a_list[i]))
    # обратный ход
    answer[length - 1] = betta[length]
    for i in range(length - 2, -1, -1):
        answer[i] = alpha[i + 1] * answer[i + 1] + betta[i + 1]
    return answer


# Решение с аппроксимацией первого порядка
def first_accuracy(h, x_list):
    a_list = []
    b_list = []
    c_list = []
    d_list = []

    # вычисление  коэффициентов для x_0
    a_list.append(0)
    b_list.append(D1 - h * F1)
    c_list.append(D1)
    d_list.append(h * E1)

    # вычисление коэффициентов для x_i (i=1,...,n-1)
    amount_points = len(x_list)
    for i in range(1, amount_points - 1):
        a_list.append(1 - fun_val(A_x, x_list[i]) * h / 2)
        b_list.append(2 - h ** 2 * fun_val(B_x, x_list[i]))
        c_list.append(1 + fun_val(A_x, x_list[i]) * h / 2)
        d_list.append(h ** 2 * fun_val(C_x, x_list[i]))

    # вычисление коэффициентов для x_n
    a_list.append(-D2)
    b_list.append(-h * F2 - D2)
    c_list.append(0)
    d_list.append(h * E2)

    # решение полученной системы методом прогонки
    y_list = slau_solution(a_list, b_list, c_list, d_list)
    return y_list


# Решение с аппроксимацией второго порядка
def second_accuracy(h, x_list):
    a_list = []
    b_list = []
    c_list = []
    d_list = []

    # вычисление  коэффициентов для x_0
    a_list.append(0)
    b_list.append(D1 - h * F1 + h / 2 * D1 * (fun_val(A_x, a) - h * fun_val(B_x, a)))
    c_list.append(D1 + fun_val(A_x, a) * D1 * h / 2)
    d_list.append(h * E1 + fun_val(C_x, a) * D1 * h ** 2 / 2)

    # вычисление коэффициентов для x_i (i=1,...,n-1)
    amount_points = len(x_list)
    for i in range(1, amount_points - 1):
        a_list.append(1 - fun_val(A_x, x_list[i]) * h / 2)
        b_list.append(2 - h ** 2 * fun_val(B_x, x_list[i]))
        c_list.append(1 + fun_val(A_x, x_list[i]) * h / 2)
        d_list.append(h ** 2 * fun_val(C_x, x_list[i]))

    # вычисление коэффициентов для x_n
    a_list.append(-D2 + fun_val(A_x, b) * D2 * h / 2)
    b_list.append(-h * F2 - D2 + h / 2 * D2 * (fun_val(A_x, b) + h * fun_val(B_x, b)))
    c_list.append(0)
    d_list.append(h * E2 - fun_val(C_x, b) * D2 * h ** 2 / 2)

    # решение полученной системы методом прогонки
    y_list = slau_solution(a_list, b_list, c_list, d_list)
    return y_list


# Вычисление нормы разности между приближенным и точным решениями
def norm(y_list, y_list_true):
    length = len(y_list)
    # вычисляется максимум разности между двумя функциями на всем отрезке интегрирования
    max_val = abs(y_list[0] - y_list_true[0])
    for i in range(1, length):
        cur_val = abs(y_list[i] - y_list_true[i])
        if cur_val > max_val:
            max_val = cur_val
    return max_val


# Промежуток, на котором находится решение
a = 0
b = 1
A_x = "-2"  # Функция при u'
B_x = "0"  # Функция при u
C_x = "exp(X) * (X ** 2 + X - 3)"  # Свободный член
# коэффициенты в первом граничном условии
F1 = 0
D1 = 1
E1 = 2
# коэффициенты во втором граничном условии
F2 = -1
D2 = 1  
E2 = exp(1) * (exp(1) - 3)
true_solution = "exp(X) * (exp(X) - X ** 2 - X + 1)"  # точное решение задачи


def main():
    h = 0
    amount_steps = 1  # Количество частей, на которых разбивается исходный промежуток

    # Вывод условий задачи
    print("Решить краевую задачу:")
    print("u''(X)   +   {0} * u'(X)   +   {1} * u(X)  =  {2}".format(A_x, B_x, C_x))
    print("на промежутке [{0}; {1}]".format(a, b))
    print("при условиях:")
    print("{0} * u(a) + {1} u'(a) = {2}".format(F1, D1, E1))
    print("{0} * u(a) + {1} u'(a) = {2}".format(F2, D2, E2))
    # Ввод данных
    amount_steps = float(input("Количество частей, на которые разбивается отрезок интегрирования: "))

    h = (b - a) / amount_steps  # Шаг
    x_list = []     # список значений аргумента
    cur_x = a
    while cur_x < b + h:    # Формирование списка значений аргумента
        x_list.append(cur_x)
        cur_x += h

    # получение значений функции-решения
    y_true = [fun_val(true_solution, cur_x) for cur_x in x_list]    # точного решения
    y_list1 = first_accuracy(h, x_list)     # приближенного 1-ым порядком
    y_list2 = second_accuracy(h, x_list)    # приближенного 2-ым порядком

    # Рисование графиков в одной системе координат
    fig, axes = pylab.subplots()
    axes.plot(x_list, y_true, label='Точное решение')
    axes.plot(x_list, y_list1, label='O(h)')
    axes.plot(x_list, y_list2, label='O(h ** 2)')
    axes.grid(True)
    axes.legend()

    norm_val = norm(y_list1, y_true)

    # Высчитываем норму разности при разных значениях шага
    h_vals = [0.2, 0.1, 0.05, 0.04, 0.02, 0.01, 0.008, 0.005, 0.004, 0.002, 0.001]      # величины шага
    #h_vals = [0.01, 0.008, 0.005, 0.004, 0.002, 0.001, 0.0008, 0.0005, 0.0004, 0.0002, 0.0001]
    norm_vals = ([], [])    # значения нормы разности
    for h in h_vals:
        # Формирование значений x-ов
        x_list = []  # список значений аргумента
        cur_x = a
        while cur_x < b + (h/10):  # Формирование списка значений аргумента
            x_list.append(cur_x)
            cur_x += h
        # получения значений функции
        y_list1 = first_accuracy(h, x_list)  # приближенного 1-ым порядком
        y_list2 = second_accuracy(h, x_list)  # приближенного 2-ым порядком
        y_true = [fun_val(true_solution, cur_x) for cur_x in x_list]    # точного решения
        # Вычисление нормы разности
        norm_vals[0].append(norm(y_list1, y_true))
        norm_vals[1].append(norm(y_list2, y_true))

    # рисования норм разности в другой системе координат

    figur, axes1 = pylab.subplots()
    axes1.plot(h_vals, norm_vals[0], label="норма разности для O(h)")
    axes1.plot(h_vals, norm_vals[1], label="норма разности для O(h ** 2)")
    axes1.grid(True)
    axes1.legend()

    # Вывод графиков
    pylab.show()


if __name__ == "__main__":
    main()
