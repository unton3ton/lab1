# coding: utf8
# - для восприятия кириллицы


import math # подключение библиотеки для мат.функций
from time import clock # для подсчёта времени работы 


def left_rectangle_help(f,a,b,n): # вспомогательная функция
    h = (b-a)/n # шаг
    s=0.0 # кладём нуль с переменную s
    for i in range(n):
        s += f((a*(n-i)+b*i)/n) # вычисляем сумму сторон прямоугольников
    return s*h # вычисляем площадь прямоугольников под кривой


def right_rectangle_help(f,a,b,n): # вспомогательная функция
    h = (b-a)/n # шаг
    s=0.0 # кладём нуль с переменную s
    for i in range(n):
        s += f((a*(n-i-1)+b*(i+1))/n) # вычисляем сумму сторон прямоугольников
    return s*h # вычисляем площадь прямоугольников под кривой


def medium_rectangle_help(f,a,b,n): # вспомогательная функция
    h = (b-a)/n # шаг
    s=0.0 # кладём нуль с переменную s
    for i in range(n):
        s += f((a*(n-i)+b*i)/n +h/2) # вычисляем сумму сторон прямоугольников
    return s*h # вычисляем площадь прямоугольников под кривой


def trapecium_help(f,a,b,n): # вспомогательная функция
    h = (b-a)/n # шаг
    s=0.0 # кладём нуль с переменную s
    for i in range(n):
        s += (f((a*(n-i)+b*i)/n)+f((a*(n-i)+b*i)/n+h))/2 # вычисляем полусумму
        #сторон трапеций
    return s*h # вычисляем площадь трапеций под кривой
    

def simpson_help(f,a,b,n): # вспомогательная функция
    h = (b-a)/n # 
    s=0.0 # 
    for i in range(n):
        s += (f((a*(n-i)+b*i)/n)+4*f((a*(n-i)+b*i)/n+h/2) + f((a*(n-i)+b*i)/n+h))/6 
    return s*h


def integrate(f, int_help, a, b, eps):# основная функция для
    #вычисления интеграла с выбором вспомогательной функции
    # для каждого метода
    n = 1 
    s1 = int_help(f,a,b,n) # считаем площадь для n разбиений
    s2 = int_help(f,a,b,n*2)# считаем площадь для 2n разбиений
    while (abs(s2 - s1) > eps):# условие сравнения с заданной точностью
        #для поиска оптимального разбиения
        n *= 2 # увеличиваем разбиение в 2 раза
        s1, s2 = s2, int_help(f,a,b,n*2) # кладём в s1 значение s2,
        #а в s2 значение с дважды удвоенным разбиением
    return s2 # возвращаем s2


####################################################################################


def polinom_help(f,a,b,n):# вспомогательная функция для
    # нахождения коэффициентов полиномиального разложения 
    h = (b-a)/n # шаг
    A = [[0] * (n + 1) for i in range(n + 1)] # создаём пустую мтр под
    # запись точек разбиения х_i^k, где i,k от 0 до n
    B = [0] * (n + 1)# создаём пустой список под запись столбца
   # свободных членов f(х_i)
    for i in range(n+1):# до n + 1, т.к задаётся степень n, а
       # ей соответствует n+1 точек
        p = a + h * i # разбиение по точкам отрезка от а до b
        B[i] = f(p)# формирование столбца свободных членов
        for j in range(n+1):
            A[i][j] = p ** j # формирование мтр
    return solve(A, B)# возвращает коэффициенты разложения функции 


def polinom(f,a,b,n):# основная функция для интегрирования
    s=0.0
    v = polinom_help(f, a, b, n)# записываем вычисленные коэф разложения
    for k in range(n+1):
        s += (v[k]*(b**(k+1)-a**(k+1)))/(k+1) # вычисляем сумму
       # соответствующую интегралу
    return s # выдаёт значение суммы


def lup(A1):# функция для lup-разложения мтр а1
    A = [[x for x in row] for row in A1]# копирование мтр а1 в а
    n=len(A)# число строк в мтр a
    pi=[i for i in range(n)]# создание массива для мтр перестановок Р
    for k in range(n):
        p=0.0
        k1=k
        for i in range(k,n):# для i=k от 0 до n-1
            if abs(A[i][k])>p:
                p=abs(A[i][k])
                k1=i
        if p==0.0:
            print("error")# выдает ошибку при обнаружении вырожденной мтр
        pi[k],pi[k1]=pi[k1],pi[k]# обмен значений pi[k] и pi[k1]
        A[k],A[k1]=A[k1],A[k]# обмен значений a[k][i] и a[k1][i]
        for i in range(k+1, n):
            A[i][k]=A[i][k]/A[k][k]
            for j in range(k+1,n):
                A[i][j]=A[i][j]-A[i][k]*A[k][j]
    L=[[0]*n for i in range(n)]# создание пустых мтр L,U
    U=[[0]*n for i in range(n)]
    for i in range(n):
        for j in range(n):
            if i>j:
                L[i][j]=A[i][j]
            else:
                U[i][j]=A[i][j]
    for i in range(n):
        L[i][i] = 1.0
    return pi, L, U

def lup_solve(L,U,pi,B):# функция для нахождения корней х с помощью lup
    n=len(L)# число строк в мтр L
    x=[0]*n# создание пустых векторов x,y для записи решения
    y=[0]*n
    for i in range(n):
        s = sum(L[i][j] * y[j] for j in range(i))# суммирование по j до i-1
        y[i] = B[pi[i]] - s# прямая подстановка
    for i in reversed(range(n)): #разворот для прохода от n до 0
        s = sum(U[i][j] * x[j] for j in range(i + 1, n))# суммирование от j=i+1 до n
        x[i] = (y[i] - s) / U[i][i]#  обратная подстановка
    return x# выдает решение х

def solve(A,B):# функция для решения слау
    pi, L, U = lup(A)# вычисление мтр с помощью функции lup() для мтр а
    #print(pi, L, U)
    return lup_solve(L, U, pi, B) # передача вычисленных
              #матриц из lup-разложения в слау


################################################################################


def lagrange(f, a, b, n):
    #  Интегрирование полиномом Лагранжа степени n
    
    # создаём равномерное разбиение отрезка [-1, 1]
    ys = [-1 + 2 * float(i) / n for i in range(n + 1)]# 1/2 = 0, а 1.0/2 = 0.5
    
    # считаем веса:
    c = weights(ys)

    # считаем ответ
    x = lambda y: (b + a) / 2 + y * (b - a) / 2 # отображает разбиение ys на
                                                # промежуток интегрирования
                                                
                                                # lambda создает объект и
                                                # возвращает его в виде результата
    s = sum(c[i] * f(x(y)) for i, y in enumerate(ys))
    #enumerate - ф-ция позволяет
    #раскладывать список на порядковые
    #номера и соответсвующие им элементы
    return s * (b - a)


def weights(ys):
    # Рассчитывает веса
    
    result = [] # создаём пустой список
    for i, y in enumerate(ys):
        # сначала формируем список [yj], j ≠ i
        yj = ys[:i] + ys[i + 1:]

        # определяем коэффициенты числителя
        p = polynomial_from_roots(yj)

        # определяем знаменатель
        value = polynomial_value(p, y)

        # рассчитываем вес и добавляем его в список весов
        result.append(polynomial_integrate(p) / value / 2)  # append добавляет
                                                            # переданный объект
                                                            # в существующий список
    return result


def polynomial_from_roots(ys):
    
    # Вычисляет список коэффициентов многочлена, имеющего корни ys
    
    result = [1]
    for y in ys:
        result = polynomial_multiply(result, y)
    return result


def polynomial_multiply(p, a):
    
     #  Вычисляет коэффициенты многочлена, который получается при умножении
     #  многочлена с коэффициентами p на (x-a)

     #  p -- список коэффициентов [a0, a1, … , an]
    
    n = len(p)
    result = []

    # начинаем формировать результат:
    # 1. добавляем свободный член
    result.append(-a * p[0])

    # 2. добавляем все коэффициенты, кроме самого старшего
    for i in range(1, n):
        result.append(p[i-1] - a * p[i])

    # 3. добавляем старший коэффициент
    result.append(p[-1])
    return result


def polynomial_value(p, x):
    # Вычисляет значение многочлена с коэффициентами p в точке x
    return sum(a * x ** i for i, a in enumerate(p))


def polynomial_integrate(p):
    # Интегрирует многочлен с коэффициентами p на [-1, 1]
    return sum(2 * a / (i + 1) for i, a in enumerate(p) if i % 2 == 0)
    # i % 2 == 0 проверка на чётность (% - остаток)


################################################################################


def legendre(f,a,b):
    s = 0.0
    h = (b-a)/2 # шаг
    n=8 # порядок
    c=[0.101228536290376, 0.222381034453374,
       0.313706645877887, 0.362683783378362,
       0.362683783378362, 0.313706645877887,
       0.222381034453374, 0.101228536290376] # табличные значения весов
    
    x=[-0.960289856497536, -0.796666477413627,
       -0.525532409916329, -0.183434642495650,
       0.183434642495650, 0.525532409916329,
       0.796666477413627, 0.960289856497536] # для соответствующих
                                             # табличных узлов
    for i in range(n-1):
        s += f(((b+a)+(b-a)*x[i])/2)*c[i] # для интегрирования на отрезке [a,b]
    return s*h


#print('left_rectangles',integrate(lambda x: math.sin(x), left_rectangle_help, 0.0, math.pi, 0.001))
#print('right_rectangles',integrate(lambda x: math.sin(x), right_rectangle_help, 0.0, math.pi, 0.001))
#print('medium_rectangles',integrate(lambda x: math.sin(x), medium_rectangle_help, 0.0, math.pi, 0.001))
#print('trapecium',integrate(lambda x: math.sin(x), trapecium_help, 0.0, math.pi, 0.001))
#print('simpson',integrate(lambda x: math.sin(x), simpson_help, 0.0, math.pi, 0.001))
#print('polinom', polinom(lambda x: math.sin(x), 0.0, math.pi, 90))
#print('lagrange', lagrange(lambda x: math.sin(x), 0.0, math.pi, 7))
#print('gauss-legendre', legendre(lambda x: math.sin(x), 0.0, math.pi))


t1 = clock()
res1 = integrate(lambda x: 1/math.sqrt(x-1), left_rectangle_help, 1.0000001, 5.0, 0.001) # math.fabs(x)==|x|
t2 = clock()

t3 = clock()
res2 = integrate(lambda x: 1/math.sqrt(x-1), right_rectangle_help, 1.0, 5.0, 0.001)
t4 = clock()

t5 = clock()
res3 = integrate(lambda x: 1/math.sqrt(x-1), medium_rectangle_help, 1.0, 5.0, 0.001)
t6 = clock()

t7 = clock()
res4 = integrate(lambda x: 1/math.sqrt(x-1), trapecium_help, 1.0000001, 5.0, 0.001)
t8 = clock()

t9 = clock()
res5 = integrate(lambda x: 1/math.sqrt(x-1), simpson_help, 1.0000001, 5.0, 0.001)
t10 = clock()

t11 = clock()
res6 = polinom(lambda x: 1/math.sqrt(x-1), 1.009, 5.0, 128)
t12 = clock()

t13 = clock()
res7 = lagrange(lambda x: 1/math.sqrt(x-1), 1.01, 5.0, 17)
t14 = clock()

t15 = clock()
res8 = legendre(lambda x: 1/math.sqrt(x-1), 1.0000001, 5.0)
t16 = clock()


print 'left_rectangle: %f by %.2f ms' % (res1, (t2-t1) * 1000)
print 'right_rectangle: %f by %.2f ms' % (res2, (t4-t3) * 1000)
print 'medium_rectangle: %f by %.2f ms' % (res3, (t6-t5) * 1000)
print 'trapecium: %f by %.2f ms' % (res4, (t8-t7) * 1000)
print 'simpson: %f by %.2f ms' % (res5, (t10-t9) * 1000)
print 'polinom: %f by %.2f ms' % (res6, (t12-t11) * 1000)
print 'lagrange: %f by %.2f ms' % (res7, (t14-t13) * 1000)
print 'gauss-legendre: %f by %.2f ms' % (res8, (t16-t15) * 1000)
