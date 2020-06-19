# Частичная проблема собственных чисел (два максимальных по модулю собственных числа, собственное число, ближайшее к заданному)
import math
import sys
from copy import copy, deepcopy


def inp_matrix():  # Ввод матрицы A
    read_f = open("inputA.txt", "r")
    matrix = []
    for line in read_f:
        matrix.append([float(x) for x in line.split()])
    read_f.close()

    # CHECK FOR CORRECT SIZE OF MATRIX
    size = len(matrix)
    for l in matrix:
        if len(l) != size:
            return []
    return matrix


def show(matrix_a):  # Показ матрицы
    print("\nИсходная матрица:")
    for line in matrix_a:
        for el in line:
            print("{0:10.3f}".format(el),end="\t")
        print()


# Скалярное произведение 2х векторов
def scalar_product(a, b):
    # Сумма произведений соответсвующих координат
    sum = 0
    for i in range(len(a)):
        sum += a[i] * b[i]

    return sum


# Нормализация вектора:
# (деление на норму)
def normalize(vector):
    norm = math.sqrt(scalar_product(vector, vector))

    for i in range(len(vector)):
        vector[i] = vector[i] / norm
    return vector


# Умножение матрицы на вектор
def multiply(matrix, vector):
    res = []

    for i in range(len(matrix)):
        sum = 0
        for j in range(len(matrix[0])):
            sum += matrix[i][j] * vector[j]
        res.append(sum)

    return res


def residual_num(matrix, eigenvalue, eigenvector):
    mul = multiply(matrix, eigenvector)

    e = deepcopy(eigenvector)

    for i in range(len(eigenvector)):
        e[i] *= eigenvalue

    for i in range(len(eigenvector)):
            e[i] -= mul[i]

    return e


# Вычисление максимального собственного числа
def compute_eigenvalue(matrix, x, epsilon):
    lambd = 0   # Собственное число
    k = 0       # Число итераций
    next_x = x  # Текущий вектор х

    while True:
        k += 1
        prev_lambd = copy(lambd)  # Собственное число на предыдущей итерации

        # Запоминание х предыдущей итерации
        # и его нормирование
        previos_x = normalize(next_x)

        # Вычисление нового х умножением матрицы на х:
        # x(n+1) = Ax(n)
        next_x = multiply(matrix, previos_x)

        # Вычисление собственного числа на текущей итерации
        lambd = scalar_product(next_x, previos_x) / scalar_product(previos_x, previos_x)

        if abs(lambd - prev_lambd) < epsilon:
            break

    return lambd, previos_x, k


# Нахождение второго максимального собственного числа
def compute_second_eigenvalue(matrix, y, epsilon, e, g):
    lambd = 0   # Собственное число
    k = 0       # Число итераций
    next_y = y  # Текущий вектор у

    while True:
        k += 1
        prev_lambd = copy(lambd)  # Собственное число на предыдущей итерации

        previos_y = copy(next_y)  # Предыдущий вектор у

        # Вычисление нового у, умножением матрицы на у:
        # у(n+1) = Aу(n)
        next_y = multiply(matrix, previos_y)

        # Корректировка у (стр 57)
        numerator = scalar_product(next_y, g)
        denominator = scalar_product(e, g)
        next_y = [next_y[i] - (numerator / denominator * e[i]) for i in range(len(next_y))]

        # Вычисление собственного числа на текущей итерации
        lambd = scalar_product(next_y, previos_y) / scalar_product(previos_y, previos_y)

        if abs(lambd - prev_lambd) < epsilon:
            break

    return lambd, previos_y, k


# Обратные итерации, нахождение собственного числа близкого к shift
def back_run(matrix, x, epsilon, shift):
    lambd = 0       # Собственное число
    k = 0           # Число итераций
    next_x = x      # Текущий вектор х

    while True:
        k += 1
        prev_lambd = copy(lambd)  # Собственное число на предыдущей итерации

        # Запоминание х предыдущей итерации
        # и его нормирование,
        previos_x = normalize(next_x)

        # Следующий х находится решением системы уравнений:
        # Bx(n+1) = x(n)
        next_x = solve_gauss(deepcopy(matrix), deepcopy(previos_x))

        # Вычисление собственного числа на текущей итерации
        # Обратите внимание что тут дробь перевернутая
        numerator = scalar_product(previos_x, previos_x)
        denominator = scalar_product(next_x, previos_x)
        lambd = shift + numerator / denominator

        if abs(lambd - prev_lambd) < epsilon:
            break

    return lambd, previos_x, k


def solve_gauss(matrix, free_members, with_main_elem=True):
    # Размерность матрицы
    size = len(matrix)

    # Единичная матрица
    identity_matrix = [[1 if i == j else 0 for i in range(size)]
                       for j in range(size)]

    # Приведение к треугольному виду
    # k - строка
    # i - столбец
    for i in range(size):
        # Проверка на вырожденность матрицы
        # Если на главной диагонали 0,
        # ищем в столбце строку с не 0 элементом
        # и меняем их местами
        if matrix[i][i] == 0:
            for k in range(i, size):
                if matrix[k][i] != 0:
                    matrix[i], matrix[k] = \
                        deepcopy(matrix[k]), deepcopy(matrix[i])

                    free_members[i], free_members[k] = \
                        copy(free_members[k]), copy(free_members[i])

                    identity_matrix[i], identity_matrix[k] = \
                        deepcopy(identity_matrix[k]), deepcopy(identity_matrix[i])

                    break
            # Если в столбце нет не нулевого элемента,
            # матрица не вырожденная
            else:
                print("Матрица вырождена.")
                return

        # Поиск строки с максимальным элементом в столбце
        if with_main_elem:
            max_elem = abs(matrix[i][i])
            max_row = i
            is_changed = False

            for k in range(i, size):
                if abs(matrix[k][i]) > max_elem:
                    max_elem = abs(matrix[k][i])
                    max_row = k
                    is_changed = True

            # Перестановка строк
            if is_changed:
                matrix[i], matrix[max_row] = \
                    deepcopy(matrix[max_row]), deepcopy(matrix[i])

                free_members[i], free_members[max_row] = \
                    deepcopy(free_members[max_row]), deepcopy(free_members[i])

                identity_matrix[i], identity_matrix[max_row] = \
                    deepcopy(identity_matrix[max_row]), deepcopy(identity_matrix[i])

        # Отнимаем из строки k строку i * на С
        for k in range(i + 1, size):
            c = matrix[k][i] / matrix[i][i]

            for j in range(size):
                matrix[k][j] -= c * matrix[i][j]
                identity_matrix[k][j] -= c * identity_matrix[i][j]

            free_members[k] -= c * free_members[i]

    # Решение
    answer = gaus_back_run(matrix, free_members)

    return answer


# Обратный ход метода Гауса
def gaus_back_run(matrix, free_members):
    size = len(matrix)

    # Вектор решения
    x = [0 for i in range(size)]

    for i in range(size - 1, -1, -1):
        # х = свободный член / коэффицент при х
        x[i] = free_members[i] / matrix[i][i]

        # В каждой строке отнимаем от свободного члена
        # Полученный х * на его коэффицент
        for k in range(i - 1, -1, -1):
            free_members[k] -= matrix[k][i] * x[i]

    return x


def main():
    file_name = "inputParams.txt"
    matrix_a = []

    # Чтение входных данных из файла
    with open(file_name) as file:
        # Начальный вектор х
        x = list(map(float, file.readline().split()))
        # Заданное собственное число
        lambd = float(file.readline())
        # Точность
        epsilon = float(file.readline())

    # Матрица
    matrix_a = inp_matrix()

    amount_digits = len(str(int(1 / epsilon)))

    start_vector = deepcopy(x)

    show(matrix_a)

    # Первое число
    first_eigenval, first_eigenvector, k = compute_eigenvalue(matrix_a, x, epsilon)
    first_eigenvector = normalize(first_eigenvector)

    print(('\n\nПервое максимальное собственное число: {0:.' + str(amount_digits) + 'f}').format(first_eigenval))
    print("Вектор соответствующий этому собственному значению: ")
    for j in range(len(first_eigenvector)):
        if j == len(first_eigenvector) // 2:
            print(('h = \t{0:.' + str(amount_digits) + 'f}').format(first_eigenvector[j]))
        else:
            print(('\t\t{0:.' + str(amount_digits) + 'f}').format(first_eigenvector[j]))
    print("Невязка: \n\t", residual_num(matrix_a, first_eigenval, first_eigenvector))
    print("Число итераций: \t", k)

    # Подготовка к вычислению 2го числа

    # Получаем а транспонированную
    a_trans = [[matrix_a[i][j] for i in range(len(matrix_a))] for j in range(len(matrix_a))]
    # Находим для нее собственный вектор g
    first_eigenval, g, k = compute_eigenvalue(a_trans, x, epsilon)
    g = normalize(g)
    # Вычисляем начальный вектор y (стр 57)
    numerator = scalar_product(x, g)
    denominator = scalar_product(first_eigenvector, g)
    y = [x[i] - (numerator / denominator * first_eigenvector[i])
         for i in range(len(x))]

    # Второе число
    second_eigenval, y_iter, n = compute_second_eigenvalue(matrix_a, y, epsilon, first_eigenvector, g)
    second_eigenvector = normalize(y_iter)

    print(('\n\nВторое максимальное собственное число: {0:.' + str(amount_digits) + 'f}').format(second_eigenval))
    print("Вектор соответствующий этому собственному значению: ")
    for j in range(len(second_eigenvector)):
        if j == len(second_eigenvector) // 2:
            print(('h = \t{0:.' + str(amount_digits) + 'f}').format(second_eigenvector[j]))
        else:
            print(('\t\t{0:.' + str(amount_digits) + 'f}').format(second_eigenvector[j]))
    print("Невязка: \n\t", residual_num(matrix_a, second_eigenval, second_eigenvector))
    print("Число итераций: ", n)

    # Нахождение собственного числа ближайщего к заданному
    # Матрица B = A - lambda*E
    matrix_b = deepcopy(matrix_a)
    for i in range(len(matrix_b)):
        matrix_b[i][i] -= lambd

    eigenval, eigenvector, m = back_run(matrix_b, start_vector, epsilon, lambd)

    print("__________________________________________________________________)")
    print(('\nБлижайщее к {0} собственное число: {1:.' + str(amount_digits) + 'f}').format(lambd, eigenval))
    print("Вектор соответствующий этому собственному значению: ")
    for j in range(len(eigenvector)):
        if j == len(eigenvector) // 2:
            print(('h = \t{0:.' + str(amount_digits) + 'f}').format(eigenvector[j]))
        else:
            print(('\t\t{0:.' + str(amount_digits) + 'f}').format(eigenvector[j]))
    print("Невязка: \n\t", residual_num(matrix_b, eigenval, eigenvector))
    print("Число итераций: ", m)


if __name__ == '__main__':
    main()
