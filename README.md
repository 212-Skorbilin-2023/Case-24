# Case-24
Программа в папке program решает СЛУ методом отражений (24)

В файлах inoutput.cpp (.h) содержатся функции ввода/вывода матриц, а также функция расчёта правой части СЛУ.
В файлах operations.cpp (.h) содержатся функции произведения матриц (в том числе и на матрицу отражения), вычисления нормы и сама функция решения СЛУ методом отражений.

Чтобы собрать программу необходимо вызвать make файл Makefile.txt
Команда в консоли: make -f Makefile.txt 

Команда в консоли для запуска программы: mycode.exe n m k filename
где n — размер матрицы, m — максимальное количество строк и столбцов матрицы, которое нужно выводить, k — режим работы программы (1 — 4: матрица заполняется по соответствующим формулам, 0: матрица считывается из файла), filename — (нужно только в случае, если k = 0) имя файла, из которого считается матрица.

Для параллельной программы: четвёртый параметр при запуске программы это желаемое количество потоков.

В программе и вектора и матрицы представлены двумерными массивами для универсальной работы некоторых функций.
