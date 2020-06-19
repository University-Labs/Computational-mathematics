# Решение СЛАУ методом Гаусса
import copy
import msvcrt


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


def make_triangle(matrix_a, free_b, is_main_elem):  # преобразование матрицы A в треугольную
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
    print("Треугольная матрица имеет вид:")
    show_system(matrix_a, free_b)
    return amount_changing


def find_decision(matrix_a, free_b):  # нахождение решений x1,...xm (обратный ход Гаусса)
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


def show_system(matrix_a, free_b):  # печать системы уравнений
    length = len(matrix_a)
    for i in range(0, length):
        for j in range(0, length):
            if j != length - 1:
                print("{:10.5f} * x{:1d} + ".format(matrix_a[i][j], j + 1), end="")
            else:
                print("{:10.5f} * x{:1d}".format(matrix_a[i][j], j + 1), end="")
        print(" = {:10.5f}".format(free_b[i]))


def find_residual(matrix_a, free_b, decision):  # нахождение величины невязки
    r = []
    length = len(matrix_a)
    for i in range(0, length):
        num = free_b[i]  # изначально полагаем, что она равна правой части
        for j in range(0, length):  # затем отнимаем сумму произведений из левой части
            num -= matrix_a[i][j] * decision[j]
        r.append(num)
    for i in range(length):
        print("Невязка r{} = {}".format(i+1, r[i]))


def find_determinant(matrix_a, change_decision):  # нахождение определителя
    det = 1
    for i in range(len(matrix_a)):  # определитель считается произведением элементов на диагонали
        det *= matrix_a[i][i]
    # Количество перестановок строк
    if change_decision % 2 == 1:
        det = -det
    return det


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
        decisions = find_decision(matrix_a, right_part)  # получаем столбец в обратной матрице
        for j in range(length):
            matrix_reverse[j][i] = decisions[j]


# Ввод матрицы и столбца свободных членов из файлов
print("Матрица и столбец свободных членов считываются из файла")
A = inp_system()
B = inp_freedata()
A_copy = copy.deepcopy(A)  # сохраняем копии начальных матриц
B_copy = copy.deepcopy(B)
# Проверка входных данных на корректность
if len(A) == 0 or len(B) == 0 or len(A) != len(B):
    print("Invalid input")
    exit(1)
length_system = len(A)  # размеры матрицы
print("Ваша система СЛАУ:")
show_system(A, B)  # печать системы уравнений
print("\n_________________________________________________________________________________\n")
print("Как решаем? С выбором главного элемента? Введите y, иначе n")
choose = False  # ответ пользователя по методу решения
ch = ''
while ch != 'y' and ch != 'n':
    ch = input()  # Считывание нажатой клавиши, пока не будет y или n
    if ch == "y":
        print("Система решается с выбором главного элемента")
        choose = True
    elif ch == 'n':
        print("Система решается без выбора главного элемента")
        choose = False

change_amount = make_triangle(A, B, choose)  # преведение к треугольному виду
if change_amount >= 0:  # матрица не вырожденнная
    X = find_decision(A, B)  # поиск решения
    print("_____________________________________________________________________________________")
    print("\nОтвет на задачу: ")
    write_f = open("output.txt", 'w')
    for p in range(0, length_system):
        write_f.write("x{0:d} = {1:10.5f}\n".format(p+1, X[p]))
        print("x{0:d} = {1:10.5f}".format(p+1, X[p]))
    print()
    find_residual(A_copy, B_copy, X)  # нахождение невязки
    determinant = find_determinant(A, change_amount)  # нахождение определителя
    write_f.write("Определитель матрицы A равен |A| = {0:10.5f}\n".format(determinant))
    print("Определитель матрицы A равен |A| = {0:10.5f}".format(determinant))
    A_rev = []
    find_reverse(A_copy, A_rev)  # нахождение обратной матрицы
    print("_______________________________________________________________________________________")
    print("Обратная матрица:")
    write_f.write("Обратная матрица:\n")
    for line in A_rev:
        for el in line:
            write_f.write("{:10.5f}\t".format(el))
            print("{:10.5f}".format(el), end="\t")
        print()
        write_f.write("\n")
    write_f.close()
else:
    write_f = open("output.txt", 'w')
    print("Матрица вырожденная")
    write_f.write("Матрица вырожденная")
    write_f.close()
