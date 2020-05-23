# Полная проблема собственных чисел
import math
import copy


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


def inp_params():
    with open('inputParams.txt', 'r') as f:
        return float(f.readline())
    return 0.1


def check_simmetr(matrix_a):
    length = len(matrix_a)
    for i in range(length):
        for j in range(i + 1, length, 1):
            if matrix_a[i][j] != matrix_a[j][i]:  # если элементы. симметричные
                # относительно главной диагонали
                # не равны - matrix_a не симметричная
                return False
    return True


def show(matrix_a):  # Показ матрицы
    print("\nИсходная матрица:")
    for line in matrix_a:
        for el in line:
            print("{0:10.3f}".format(el),end="\t")
        print()


def unit_matrix(size):  # Создание единичной матрицы
    u_matrix = []
    for i in range(size):
        u_matrix.append([0 if i != j else 1 for j in range(size)])
    return u_matrix


def get_m(matrix_a, index_k, index_l):  # Возвращает мю (необходимо при вычислении альфа и бетта)
    return 2 * matrix_a[index_k][index_l] / (matrix_a[index_k][index_k] - matrix_a[index_l][index_l])


def get_alpha_betta(matrix_a, index_k, index_l):  # Возвращает alpha для матрицы U
    # Если a[k][k] = a[l][l] ==> альфа = бетта = √1/2
    if matrix_a[index_k][index_k] == matrix_a[index_l][index_l]:
        alpha = betta = math.sqrt(0.5)
        return alpha, betta

    alpha = math.sqrt(0.5 * (1 + 1 / math.sqrt(1 + get_m(matrix_a, index_k, index_l) ** 2)))
    if get_m(matrix_a, index_k, index_l) > 0:
        betta = 1
    else:
        betta = -1
    betta *= math.sqrt(0.5 * (1 - 1 / math.sqrt(1 + get_m(matrix_a, index_k, index_l) ** 2)))
    return alpha, betta


def find_max_el(matrix_a):  # Поиск индексов максимального наддиагонального элемента
    index_k = 1
    index_l = 0
    size = len(matrix_a)
    for i in range(size):                # Проход по матрице
        for j in range(i + 1, size):     # Выше главной диагонали
            if abs(matrix_a[index_k][index_l]) < abs(matrix_a[i][j]):
                index_k, index_l = i, j  # И поиск максимального по модулю значения
    return index_k, index_l


def create_new_matrix(matrix_a, index_k, index_l, alpha, betta):  # Создание новой матрицы A
    matrix_c = copy.deepcopy(matrix_a)

    # Формирование матрицы C = AU
    for i in range(len(matrix_a)):
        matrix_c[i][index_k] = matrix_a[i][index_k] * alpha + matrix_a[i][index_l] * betta
        matrix_c[i][index_l] = -matrix_a[i][index_k] * betta + matrix_a[i][index_l] * alpha

    matrix_b = copy.deepcopy(matrix_c)

    # Формирование матрицы B = UtC
    for i in range(len(matrix_a)):
        matrix_b[index_k][i] = matrix_c[index_k][i] * alpha + matrix_c[index_l][i] * betta
        matrix_b[index_l][i] = -matrix_c[index_k][i] * betta + matrix_c[index_l][i] * alpha
    matrix_b[index_k][index_l] = matrix_b[index_l][index_k] = 0

    return matrix_b


def multiplication_matrix(matrix_1, matrix_2):  # Функция перемножения двух квадратных  матриц
    new_matrix = []  # Матрица - результат умножения
    if len(matrix_1) == len(matrix_2):
        size = len(matrix_1)
        for i in range(size):
            new_matrix.append([0 for j in range(size)])
            for j in range(size):
                summ = 0
                for k in range(size):
                    summ += matrix_1[i][k] * matrix_2[k][j]
                new_matrix[i][j] = summ
        return new_matrix
    else:
        return new_matrix


def sobstv_numbers(matrix_a):  # Нахождение собственных векторов для матрицы A
    matrix_a_copy = copy.deepcopy(matrix_a)
    matrix_un = unit_matrix(len(matrix_a))                      # Создание единичной матрицы
    matrix_d = unit_matrix(len(matrix_a))                       # Матрица, осуществляющая поворот
    while True:
        index_k, index_l = find_max_el(matrix_a)               # Нахождение максимального
        max = matrix_a[index_k][index_l]                       # недиагонального элемента
        if abs(max) < epsilon:                                 # Если все недиагональные элементы нули - вектора найдены
            break
        #  Нахождение Альфа, бетта, и формирование новой матрицы, уже повёрнутой
        alpha, betta = get_alpha_betta(matrix_a, index_k, index_l)
        matrix_a = create_new_matrix(matrix_a, index_k, index_l, alpha, betta)

        #  Заполнение позиций (k,k), (k,l), (l,k) и (l,l) в матрице D
        for i in range(len(matrix_a)):
            matrix_d[i] = [1 if i == j else 0 for j in range(len(matrix_a))]
        matrix_d[index_k][index_k] = matrix_d[index_l][index_l] = alpha
        matrix_d[index_k][index_l] = -betta
        matrix_d[index_l][index_k] = betta

        #  Перемножение матрицы с единичной из предыдущей итерации
        matrix_un = multiplication_matrix(matrix_un, matrix_d)

    amount_digits = len(str(int(1 / epsilon)))
    #  Вывод собственных значений и собственных векторов
    for i in range(len(matrix_a)):
        print(('Собственное значение {0} = {1:.' + str(amount_digits) + 'f}').format(i+1, matrix_a[i][i]))
        print("Вектор соответствующий этому собственному значению: ")
        for j in range(len(matrix_un)):
            if j == len(matrix_un) // 2:
                print(('h = \t{0:.' + str(amount_digits) + 'f}').format(matrix_un[j][i]))
            else:
                print(('\t\t{0:.' + str(amount_digits) + 'f}').format(matrix_un[j][i]))
        print(40 * '_')
    find_nevyazka(matrix_a_copy, matrix_a, matrix_un)


def find_nevyazka(matrix_a, sobstv_numb, sobstv_vec):  # Вычисление Невязки
    # Ah = lh
    size = len(matrix_a)
    nevyazka = []
    for i in range(size):
        nevyazka.append([0 for i in range(size)])

    for i in range(size):
        for j in range(size):
            nevyazka[i][j] = sobstv_numb[i][i] * sobstv_vec[j][i]
            for k in range(size):
                nevyazka[i][j] -= sobstv_vec[k][i] * matrix_a[j][k]

    print(40 * '#')
    print("Невязка: ")
    for i in range(size):
        for j in range(size):
            if j == size // 2:
                print("r{0} = \t {1}".format(i+1, nevyazka[i][j]))
            else:
                print("\t\t {0}".format(nevyazka[i][j]))


def main():
    matrix_A = inp_matrix()                                     # Ввод
    if check_simmetr(matrix_A) is True:                         # Проверка на симметричность
        show(matrix_A)                                          # Вывод на экран
        print(40*'_')
        sobstv_numbers(matrix_A)                                # Нахождение собственных чисел и векторов
    else:
        print("Матрица не является симметричной!")


epsilon = 0.1
if __name__ == '__main__':
    epsilon = inp_params()
    main()
