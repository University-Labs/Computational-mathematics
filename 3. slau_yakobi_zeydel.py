# Решение СЛАУ методами Якоби и Зейделя
import copy


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


# Метод Зейделя
def Seidel(matrix, freeMembers, solution, epsilon, maxIterCount):
    size = len(matrix[0])       # Размер матрицы
    tempSolution = [0] * size   # Предыдущие решения
    itterationCount = 0         # Кол-во иттераций

    while True:
        # Считаем кол-во иттераций
        itterationCount += 1

        # В tempSolution сохраним старые решения
        for i in range(size):
            tempSolution[i] = solution[i]

        for i in range(size):
            var = 0

            # До i-го подставляем новые решения
            for j in range(i):
                var += matrix[i][j] * solution[j]

            # После i-го подставляем старые решения
            for j in range(i + 1, size):
                var += matrix[i][j] * tempSolution[j]

            # Получаем новое решение (Xi)
            solution[i] = (freeMembers[i] - var) / matrix[i][i]

        # Условие окончания
        if sum((solution[i] - tempSolution[i]) ** 2 for i in range(size)) ** 0.5 < epsilon or \
            itterationCount >= maxIterCount:

            break

    return itterationCount


# Метод Якоби
def Jacobi(matrix, freeMembers, solution, epsilon, maxIterCount):
    size = len(matrix[0])       # Размер матрицы
    tempSolution = [0] * size   # Новые решения
    norm = None                 # Норма
    itterationCount = 0         # Кол-во иттераций

    while True:
        # Считаем кол-во иттераций
        itterationCount += 1

        for i in range(size):
            tempSolution[i] = freeMembers[i]
            for j in range(size):
                if i != j:
                    tempSolution[i] -= matrix[i][j] * solution[j]
            tempSolution[i] /= matrix[i][i]

        # Вычисление нормы (наибольшая разность иксов соседних итераций)
        norm = abs(solution[0] - tempSolution[0])
        for i in range(size):
            if abs(solution[i] - tempSolution[i]) > norm:
                norm = abs(solution[i] - tempSolution[i])

            # Сразу присвоим новые решения
            solution[i] = tempSolution[i]

        # Условие окончания
        if norm < epsilon or itterationCount >= maxIterCount:
            break

    return itterationCount


def main():
    matrix = list()     # Исходная матрица
    freeMembers = list()  # Свободные члены
    epsilon = None  # точность вычисления

    initial_approximation = None  # Начальное приближение
    max_IterCount = None  # Максимальное кол-во иттераций
    # Решения СЛАУ,
    # изначально содержит начальное приближение
    solution = None

    #Ввод данных из файлов
    with open('inputParams.txt', 'r') as f:
        epsilon = float(f.readline())
        max_IterCount = int(f.readline())
        initial_approximation = list(map(float, f.readline().split()))
    matrix = inp_system()
    freeMembers = inp_freedata()

    # Проверка на корректность ввода данных
    size = None
    for line in matrix:
        if size is None:
            size = len(line)
            continue
        elif size != len(line):
            print('Неверное введены данные!')
            return -1
    if len(matrix) != len(freeMembers) :
        print('Неверное введены данные!')
        return -1

    # Вывод введенной матрицы на экран
    print('Введенная матрица: ')
    show_system(matrix, freeMembers)
    amount_digits = len(str(int(1 / epsilon)))
    # Поместим в solution
    # начальное приближение
    solution = copy.copy(initial_approximation)

    # Вывод начального приближения на экран
    print('\nНачальное приближение:')
    for i in range(size):
        print('\tx{0}: {1}'.format(i + 1, initial_approximation[i]))

    iterCount = None  # Количество итераций

    # Решение СЛАУ методом Якоби
    iterCount = Jacobi(matrix, freeMembers, solution, epsilon, max_IterCount)

    # Вывод результатов метода Якоби
    print("_________________________________________________________\n")
    print('Метод Якоби:')
    print()
    print('Результат: ')
    for i in range(size):
        print(('\tx{0}: {1:.' + str(amount_digits) + 'f}').format(i + 1, solution[i]))
    print('\nКол-во иттераций: {0}'.format(iterCount))

    # Поместим в solution
    # начальное приближение
    solution = copy.copy(initial_approximation)

    # Решение СЛАУ методом Зейделя
    iterCount = Seidel(matrix, freeMembers, solution, epsilon, max_IterCount)

    # Вывод результатов метода Зейделя
    print("_________________________________________________________\n")
    print('Метод Зейделя:')
    print()
    print('Результат: ')
    for i in range(size):
        print(('\tx{0}: {1:.' + str(amount_digits) + 'f}').format(i + 1, solution[i]))
    print('\nКол-во иттераций: {0}'.format(iterCount))


if __name__ == '__main__':
    main()

# При увеличении чисел начального приближения кол-во итераций в
# методе Якоби растет гораздо быстрее, чем в методе Зейделя.
# Также при увеличении числа epsilon кол-во иттераций в методе Якоби
# ощутимо растет в сравнении с Зейделем.
