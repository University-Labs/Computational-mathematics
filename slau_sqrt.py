# Решение СЛАУ методом квадратного корня
import math
import copy


def replace_columns(matrix, num1, num2):  # переставляет колонки num1 и num2
    length = len(matrix)
    for i in range(length):
        matrix[i][num1], matrix[i][num2] = matrix[i][num2], matrix[i][num1]
def check_epsilon(num):  # проверяет значение на возможную ошибку округления
    if abs(num) < 10e-16:
        return 0
    else:
        return num


def inp_system():  # Ввод матрицы A
    read_f = open("inputA.txt", "r")
    matrix_vars = []
    for line in read_f:
        matrix_vars.append([float(x) for x in line.split()])
    read_f.close()

    # CHECK FOR CORRECT SIZE OF MATRIX
    size = len(matrix_vars)
    for line in matrix_vars:
        if len(line) != size:
            return []
    return matrix_vars


def inp_freedata():  # ввод столбца свободных членов
    read_f = open("inputB.txt", "r")
    matrix_free = []
    for line in read_f:
        matrix_free.append(float(line))
    read_f.close()
    return matrix_free


def check_symmetric(matrix_a):  # Проверка матрицы на симметричность
    length = len(matrix_a)
    for i in range(length):
        for j in range(i+1, length, 1):
            if matrix_a[i][j] != matrix_a[j][i]:  # если элементы. симметричные
                                                  # относительно главной диагонали
                                                  # не равны - matrix_a не симметричная
                return False
    return True


def show_system(matrix_a, free_b):  # печать системы уравнений
    print("Ваша система СЛАУ:")
    length = len(matrix_a)
    for i in range(0, length):
        for j in range(0, length):
            if j != length - 1:
                print("{:10.3f} * x{:1d} + ".format(matrix_a[i][j], j + 1), end="")
            else:
                print("{:10.3f} * x{:1d}".format(matrix_a[i][j], j + 1), end="")
        print(" = {:10.3f}".format(free_b[i]))


def initialize_s_d(matrix_a, free_b, matrix_s, matrix_d, numbers):  # создание и инициализация матриц S и D
    # numbers - массив, хранящий индексы итоговых переменных, меняется, если s[i][i] оказывается нулём
    # и мы переставляем столбцы и строчки матрицы A
    length = len(matrix_a)  # размер матрицы

    # Создание матриц S и D
    for i in range(length):
        matrix_s.append([0 for j in range(length)])
        matrix_d.append([0 for j in range(length)])
        numbers.append(i)
    # Инициализация матриц S и D
    j = 0
    for i in range(length):
        # Вычисляем Sii до тех пор, пока он не будет нулевым, либо матрица не закончится
        while True:
            j += 1
            var_tmp = matrix_a[i][i]
            for k in range(i):  # k = 0; k < i; k++
                var_tmp -= matrix_s[k][i] * matrix_s[k][i] * matrix_d[k][k]
            if abs(var_tmp) < 10e-16 and i + j < length:  # если на главной диагонали возник Sii = 0
                # Переставляем в матрице A столбцы i и j и строчки i и j
                # чтобы на главной диагонали оказался не 0
                numbers[i], numbers[i + j] = numbers[i + j], numbers[i]  # записываем в историю перестановок
                matrix_a[i], matrix_a[i + j] = matrix_a[i + j], matrix_a[i]
                free_b[i], free_b[i + j] = free_b[i + j], free_b[i]
                replace_columns(matrix_a, i, i + j)
                replace_columns(matrix_s, i, i + j)
            else:  # если Sii != 0, то записываем его в матрицу и переходим к дальнейшему заполнению строки
                if var_tmp == 0:
                    return -1
                matrix_d[i][i] = 1 if var_tmp >= 0 else -1  # знак элемента должен совпадать со знаком выражения
                matrix_s[i][i] = math.sqrt(check_epsilon(abs(var_tmp)))
                break
        for j in range(i+1, length):  # Вычисление элементов выше главной диагонали (тк S - верхняя треугольная)
            var_tmp = matrix_a[i][j]
            for k in range(i):
                var_tmp -= matrix_s[k][i] * matrix_s[k][j] * matrix_d[k][k]
            # Проводим дополнительную проверку в случае возможной ошибки округления
            matrix_s[i][j] = check_epsilon(var_tmp / (matrix_s[i][i] * matrix_d[i][i]))
    return 0


def find_decision(matrix_a, free_b, matrix_s, matrix_d):  # нахождение решения СЛАУ
    length = len(matrix_a)
    decisions_z = [0 for _ in range(length)]  # вектор со значениями Zi (т.е. S(T) * Z = b)
    decisions_y = [0 for _ in range(length)]  # вектор со значениями Yi (т.е. D * Y = Z)
    decisions_x = [0 for _ in range(length)]  # вектор с итоговыми решениями Xi (т.е. S * X = Y)

    for k in range(length):  # матрица вырождена
        if check_epsilon(matrix_s[k][k]) == 0:
            return []

    # Находим сначала решения Zi и Yi
    for i in range(length):
        # Zi = (Bi - (Z1*S1i + Z2*S2i + Z3*S3i + ... + Zm*Smi)) / Sii
        decisions_z[i] = free_b[i]
        for j in range(i):
            decisions_z[i] -= decisions_z[j] * matrix_s[j][i]
        # дополнительно проверяем, вдруг число будет СЛИШКОМ маленьким
        decisions_z[i] = check_epsilon(decisions_z[i] / matrix_s[i][i])
        # Yi = Zi / Dii
        decisions_y[i] = decisions_z[i] / matrix_d[i][i]
    # Теперь находим окончательные решения Xi
    for i in range(length - 1, -1, -1):
        # Xi находятся с последней строчки, тк S - верхняя треугольная матрица
        decisions_x[i] = decisions_y[i]
        for j in range(i+1, length):
            decisions_x[i] -= matrix_s[i][j] * decisions_x[j]
        decisions_x[i] = check_epsilon(decisions_x[i] / matrix_s[i][i])
    return decisions_x


def find_reverse(matrix_a, matrix_reverse):
    length = len(matrix_a)

    unit_matrix = []  # единичная матрица
    for i in range(length):
        unit_matrix.append([1 if i == j else 0 for j in range(length)])
        matrix_reverse.append([0 for j in range(length)])

    for i in range(0, length):  # приводим матрицу к диагональному виду
        # Если встретили на главной диагонали нулевой элемент
        if matrix_a[i][i] < 10e-16:
            for j in range(i + 1, length):
                if matrix_a[j][i] != 0:  # Ищем под ним строку j с не нулём на позиции i
                    matrix_a[j], matrix_a[i] = matrix_a[i], matrix_a[j]
                    unit_matrix[j], unit_matrix[i] = unit_matrix[i], unit_matrix[j]
                    break
        # Отнимаем от каждой строчки i-ую, умноженную на C, также и с правой частью
        for j in range(i + 1, length):
            c = matrix_a[j][i] / matrix_a[i][i]
            for k in range(i, length):
                matrix_a[j][k] -= matrix_a[i][k] * c
                if abs(matrix_a[j][k]) < 10e-16:  # чтобы избежать ошибок округления
                    matrix_a[j][k] = 0
            for k in range(length):
                unit_matrix[j][k] -= unit_matrix[i][k] * c
                if abs(unit_matrix[j][k]) < 10e-16:  # чтобы избежать ошибок округления
                    unit_matrix[j][k] = 0

    # Теперь совершаем обратный ход Гаусса
    for i in range(length):
        right_part = [unit_matrix[s][i] for s in range(length)]
        decision = [0 for _ in range(length)]
        for j in range(length - 1, -1, -1):
            answer = right_part[j]
            for f in range(j + 1, length):
                answer -= decision[f] * matrix_a[j][f]
            answer /= matrix_a[j][j]
            answer = check_epsilon(answer)
            decision[j] = answer  # получаем столбец в обратной матрице
        for f in range(length):
            matrix_reverse[f][i] = decision[f]


def conditionality(matrix_a):  # вычисляет число обусловленности для матрицы matrix_a
    cond = 1  # будет храниться число обусловленности ( ||A|| * ||A-1|| )
    maximum = 0  # для нахождения максимума по матрице
    length = len(matrix_a)
    for i in range(length):  # Нахождение ||A||
        summ = 0
        for j in range(length):
            summ += abs(matrix_a[i][j])
        if summ > maximum:
            maximum = summ
    cond *= maximum

    matrix_a_rev = []  # Обратная матрица для A
    find_reverse(matrix_a, matrix_a_rev)

    maximum = 0
    for i in range(length):  # Нахождение ||A-1||
        summ = 0
        for j in range(length):
            summ += abs(matrix_a_rev[i][j])
        if summ > maximum:
            maximum = summ
    cond *= maximum
    return cond


def gilbert():
    print("Число обусловленности для матрицы гильберта")
    for n in range(2, 8):  # нахождение числа обусловленности для n = 2,3,...,7
        matrix_gilbert = []  # матрица гильберта
        matrix_gilbert_rev = []  # обратная к матрице гильберта
        for i in range(n):
            matrix_gilbert.append([1/(i + j + 1) for j in range(n)])  # построение матрицы гильберта для заданного n
        # Нахождение числа обусловленности
        print("При n = {0} число обусловленности равно {1}".format(n, conditionality(matrix_gilbert)))


# Ввод матрицы и столбца свободных членов из файлов
print("Матрица и столбец свободных членов считываются из файла")
A = inp_system()
B = inp_freedata()

# Проверка входных данных на корректность
if len(A) == 0 or len(B) == 0 or len(A) != len(B):
    print("Invalid input")
    exit(1)
length_system = len(A)  # размеры матрицы
show_system(A, B)  # печать на экран системы уравнений
print("\n_________________________________________________________________________________\n")

if check_symmetric(A) is False:  # Если матрица A не симметрична
    print("A не симметрична, решается другим методом")
else:  # матрица A симметрична - продолжаем решать
    S = []  # верхняя треугольная матрица
    D = []  # диагональная матрица, с эл-ми на главной диагонали +- 1
    numeric = []  # будет хранить историю перестановок столбцов матрицы

    if initialize_s_d(copy.deepcopy(A), B, S, D, numeric) == -1:  # заполнение матриц S и D
        print("Матрица вырождена")
        exit(2)
    # Печать получившихся матриц
    print("Матрица S (верхняя треугольная):")
    for line in S:
        for hm in range(length_system):
            print("{0:10.3f}".format(line[hm]), end="\t")
        print()
    print("\n_________________________________________________________________________________\n")
    print("Матрица D:")
    for line in D:
        for hm in range(length_system):
            print("{0}".format(line[hm]), end="\t")
        print()

    X_replaced = find_decision(A, B, S, D)  # поиск решения (переменные могут быть переставлены)

    wr = open("output2.txt", "w")  # файл для записи ответа
    print("\n_________________________________________________________________________________\n")
    print("Ответ:")
    if not X_replaced:  # решение не найдено - матрица вырожденная
        print("Матрица вырожденная, решений нет")
        wr.write("Матрица вырожденная, решений нет\n")
    else:
        # преобразование индексов переменных соответственно проведённым перестановкам столбцов в numeric    n
        X = [0 for _ in range(length_system)]
        for p in range(length_system):
            X[numeric[p]] = X_replaced[p]
        # печать ответа на экран и в файл
        for p in range(0, length_system):
            wr.write("x{0:d} = {1:10.3f}\n".format(p + 1, X[p]))
            print("x{0:d} = {1:10.3f}".format(p + 1, X[p]))
        print("Обратная матрица: ")
        A_rev = []
        find_reverse(copy.deepcopy(A), A_rev)  # нахождение обратной матрицы
        for line in A_rev:
            for hm in range(length_system):
                print("{0:10.3f}".format(line[hm]), end="\t")
            print()
    wr.close()
    print("______________________________________________________________")
    gilbert()  # нахождение чисел обусловленности матрицы гильберта
