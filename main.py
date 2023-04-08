# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import math
from copy import deepcopy
import time
import matplotlib.pyplot as plt


def forward_substitution(A, b):
    x = []  # dla układu równań Ax = b
    for i in range(len(b)):
        x.append(b[i])
        for j in range(i):
            val = x[i] - (A[i][j] * x[j])
            x[i] = val
        x[i] /= A[i][i]
    return x


def backward_substitution(A, b):
    x = []  # dla układu równań Ax = b
    m = len(b) - 1
    for i in range(len(b)):
        x.append(b[m - i])
        for j in range(i):
            val = x[i] - (A[m - i][m - j] * x[j])
            x[i] = val
        x[i] /= A[m-i][m-i]
    x.reverse()
    return x


def sub_vectors(x1, x2):
    output = deepcopy(x1)
    for i in range(len(x1)):
        output[i] -= x2[i]
    return output


def multiply_matrix_by_vector(A, b):
    multiplied_vector = []  # Mnożenie macierzy przez pionowy wektor
    for i in range(len(b)):
        value = 0
        for j in range(len(A[i])):
            value += A[i][j] * b[j]
        multiplied_vector.append(value)
    return multiplied_vector  # wynikiem jest nowy pionowy wektor


def norm(v):  # obliczanie normy zgodnie ze wzorem ( p = 2 )
    norm_value = 0
    for i in range(len(v)):
        norm_value += (v[i] * v[i])
    return math.sqrt(norm_value)


def jacobi(A, x, b):
    iterations = 0
    x_prev = deepcopy(x)
    start = time.time()
    while True:
        x = []
        for i in range(len(A)):
            x.append(b[i])
            for j in range(i):
                x[i] -= (A[i][j] * x_prev[j])
            for j in range(i + 1, len(A)):
                x[i] -= A[i][j] * x_prev[j]
            x[i] /= A[i][i]
        iterations += 1
        help1 = multiply_matrix_by_vector(A, x)  # obliczanie wektora residuum
        res = sub_vectors(help1, b)
        norm_of_residuum = norm(res)
        #print(norm_of_residuum) # Sprawdzenie do zadania C
        if norm_of_residuum >= pow(10,10):  # Dla macierzy z zadania C, norma z wektora residuum ciągle rośnie
            print("Metoda Jacobiego nie zbiega się")  # Dodane zakończenie wykonywania iteracji
            return
        if norm_of_residuum <= pow(10, -9):
            break  # warunek stopu
        # print(iterations)
        x_prev = deepcopy(x)
    end = time.time()
    alg_time = end - start
    print("Metoda Jacobiego")
    print("Iteracje: " + str(iterations))
    print("Czas: " + str(alg_time))
    print(x)
    print("")
    return alg_time


def gauss_Seidl(A,x,b):
    iterations = 0
    x_prev = deepcopy(x)
    start = time.time()
    while True:
        x = []
        for i in range(len(A)):
            x.append(b[i])
            for j in range(i):
                x[i] -= (A[i][j] * x[j])
            for j in range(i + 1, len(A)):
                x[i] -= A[i][j] * x_prev[j]
            x[i] /= A[i][i]
        iterations += 1
        help1 = multiply_matrix_by_vector(A, x)  # obliczanie wektora residuum
        res = sub_vectors(help1, b)
        norm_of_residuum = norm(res)
        # print(norm_of_residuum) # Sprawdzenie do zadania C
        if norm_of_residuum >= pow(10, 10):  # Dla macierzy z zadania C, norma z wektora residuum ciągle rośnie
            print("Metoda Gaussa-Seidla nie zbiega się")  # Dodane zakończenie wykonywania iteracji
            return
        if norm_of_residuum <= pow(10, -9):
            break  # warunek stopu
        # print(iterations)
        x_prev = deepcopy(x)
    end = time.time()
    alg_time = end - start
    print("Metoda Gaussa-Seidla")
    print("Iteracje: " + str(iterations))
    print("Czas: " + str(alg_time))
    print(x)
    print("")
    #print(x)
    return alg_time


def factorization_LU(A, x, b):
    start = time.time()
    U = deepcopy(A)  # Utworzenie macierzy L i U
    L = [0] * len(A)  # macierz I
    for i in range(len(A)):
        L[i] = [0] * len(A)
    for i in range(len(A)):
        for j in range(len(A)):
            if i == j:
                L[i][j] = 1
    for i in range(len(A)-1):
        for j in range(i+1,len(A)):
            L[j][i] = U[j][i]/U[i][i]
            for k in range(i,len(A)):
                U[j][k] = U[j][k] - L[j][i] * U[i][k]
    #  rozwiązanie układu równań Ly = b
    y = forward_substitution(L, b)
    #  rozwiązanie układu równań Ux = y
    x = backward_substitution(U, y)  # wynik
    end = time.time()
    #  sprawdzenie normy z wektora residuum
    help1 = multiply_matrix_by_vector(A, x)
    res = sub_vectors(help1, b)
    norm_of_residuum = norm(res)
    alg_time = end-start
    print("Faktoryzacja LU")
    print("Norma z residuum: " + str(norm_of_residuum))
    print("Czas: " + str(end - start))
    print(x)
    return alg_time


def exerciseE_times():
    N = [100, 500, 1000, 2000, 3000]
    jacobi_times = []
    gauss_seidl_times = []
    lu_factorization_times = []
    print("Zadanie E: ")
    print("")
    for size in N:
        print("Dla wielkości N = " + str(size) + ":")
        e = 7
        f = 4
        a1 = 5 + e
        a2 = a3 = -1
        A = create_matrix_A(size, a1, a2, a3)
        b = []  # generacja wektora b
        for i in range(size):
            val = math.sin((i + 1) * (f + 1))
            b.append(val)
        print("")  # zliczanie czasów trwania dla danych N:
        x = [1] * size
        jacobi_times.append(jacobi(A, x, b))
        x = [1] * size
        gauss_seidl_times.append(gauss_Seidl(A, x, b))
        x = [1] * size
        lu_factorization_times.append(factorization_LU(A, x, b))
        print("")
    #  Utworzenie wykresu czasów w zależności od N, dla 3 metod
    plt.plot(N, jacobi_times, label='Jacobi')
    plt.plot(N, gauss_seidl_times, color='g', label='Gauss-Seidl')
    plt.plot(N, lu_factorization_times, color='r', label='LU')
    plt.title("Zależność czasu trwania algorytmów od liczby niewiadomych")
    plt.xlabel("Liczba niewiadomych N")
    plt.ylabel("Czas trwania [s]")
    plt.legend()
    plt.show()


def create_matrix_A(N, a1, a2, a3):
    A = [0] * N
    for i in range(N):
        A[i] = [0] * N
    for i in range(N):
        for j in range(N):  # wypełnienie macierzy A zgodnie z instrukcją
            if i == j:
                A[i][j] = a1
            if i + 1 == j or i == j + 1:
                A[i][j] = a2
            if i + 2 == j or i == j + 2:
                A[i][j] = a3
    return A


# Indeks 184751 - e = 7, f = 4
if __name__ == '__main__':
    # Zadanie A
    N = 951  # Generacja macierzy A
    e = 7
    f = 4
    a1 = 5 + e
    a2 = a3 = -1
    A = create_matrix_A(N, a1, a2, a3)

    b = []  # generacja wektora b
    for i in range(N):
        val = math.sin((i + 1) * (f + 1))
        b.append(val)
    # print(b)
    x = [1] * N
    # Zadanie B

    print("Zadanie B:")
    print("")
    jacobi(A, x, b)
    x = [1] * N   # wypełnienie wektora x jedynkami, dla kolejnej metody
    gauss_Seidl(A, x, b)

    # Zadanie C
    print("Zadanie C:")
    a1 = 3
    a2 = a3 = -1
    A = create_matrix_A(N, a1, a2, a3)  #  Utworzenie macierzy dla nowego a1
    print("")
    x = [1] * N
    jacobi(A, x, b)  #  sprawdzenie metody Jacobiego dla nowej macierzy
    x = [1] * N
    gauss_Seidl(A, x, b)  #  sprawdzenie metody Gaussa-Seidla dla nowej macierzy
    x = [1] * N
    print("")
    print("Zadanie D:")
    print("")
    factorization_LU(A, x, b)  # rozwiązanie układu równań za pomocą faktoryzacji LU
    print("")
    exerciseE_times()  # Utworzenie wykresów zależności czasu trwania algorytmów

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
