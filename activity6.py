from tkinter import StringVar, IntVar
from numpy import float_
from calculator import SystemEquations
from sympy import eye, zeros

Ax = []
b = []

def guiSysLinEq(api):
    guiPI = api.gui(('Sistemas de Ecuaciones Lineales', 'black', 'gray', '#009ece'), plot=False)          # 0
    e = [guiPI, guiPI.getList('Metodos de Soulucion', 25, 17, 25, 45, ('E.G.Sustitución hacia atrás pivoteo',
                              'E.G.Pivoteo Máximo de Columnas', 'E.G.Pivoteo Escalado de Columna', 'Factorizacion LU',
                              'Factorizacion de Choleski'), StringVar())]
    e.append(guiPI.getEntry('Matriz de Coeficientes:', 200, 17, 345, 17, var=StringVar()))                  # 2
    e.append(guiPI.getEntry('Matriz de Constantes:', 200, 45, 345, 45, var=StringVar()))                    # 3
    e.append(guiPI.getEntry('Cifras Significativas:', 500, 35, 630, 35, var=IntVar()))                      # 4
    guiPI.setButton('Leer Sistema de Ecuaciones', 800, 17, readSysEquations, e)
    guiPI.setButton('Resolver Sistema de Ecuaciones', 800, 55, solutionLinSys, e)

def readSysEquations(e):
    e[0].output.config(state='normal')
    Ax.clear()
    b.clear()
    if e[2].get() != '' and e[3].get() != '':
        Ei = e[2].get().split()
        constants = float_(e[3].get().split())
        if len(Ei) == len(constants):
            for i in range(0, len(constants), 1):
                Ax.append(float_(Ei[i].split(',')))
                b.append(constants[i])
            if len(Ax[0]) != len(Ax):
                Ax.clear()
                b.clear()
                e[0].output.insert('end', 'Error: El numero de variables es diferente al numero de las constantes\n')
            else:
                e[0].output.insert('end', '[### Lectura del sistema de ecuaciones completa ###]\n')
        else:
            e[0].output.insert('end', 'Error: El numero de Ecuaciones es diferente al numero de las constantes\n')
    else:
        e[0].output.insert('end', 'Error: Datos no ingresados\n')
    e[0].output.config(state='disabled')

def solutionLinSys(e):
    e[0].output.config(state='normal')
    if e[1].get() != '' and len(b) != 0:
        sysEq = SystemEquations(Ax, b)
        e[0].output.insert('end', sysEq)
        if e[1].get() == 'Factorizacion LU':
            e[0].output.insert('end', '---> Factorizacion LU\n')
            if sysEq.lu(e[4].get()):
                e[0].output.insert('end', 'L:\n')
                printM(sysEq.l, sysEq.row, e[0].output)
                e[0].output.insert('end', 'U:\n')
                printM(sysEq.u, sysEq.row, e[0].output)
                yn(sysEq, e[4].get(), e[0].output)
        elif e[1].get() == 'Factorizacion de Choleski':
            e[0].output.insert('end', '---> Factorizacion de Choleski\n')
            sysEq.LLt(e[4].get())
            e[0].output.insert('end', 'L:\n')
            printM(sysEq.l, sysEq.row, e[0].output)
            e[0].output.insert('end', 'Lt:\n')
            printM(sysEq.u, sysEq.row, e[0].output)
            yn(sysEq, e[4].get(), e[0].output)
        elif e[1].get() == 'E.G.Pivoteo Máximo de Columnas':
            e[0].output.insert('end', '---> Eliminacion Gaussiana Pivoteo Máximo de Columnas\n')
            if partialPivoting(sysEq, e[4].get(), e[0].output):
                backwardSubstitution(sysEq, e[4].get(), e[0].output)
        elif e[1]. get() == 'E.G.Pivoteo Escalado de Columna':
            e[0].output.insert('end', '---> Eliminacion Gaussiana Pivoteo Escalado de Columna\n')
            if scaledPartialPivoting(sysEq, e[4].get(), e[0].output):
                backwardSubstitution1(sysEq, e[4].get(), e[0].output)
        else:
            e[0].output.insert('end', '---> Eliminacion Gaussiana Sustitución hacia atrás pivoteo\n')
            if reduction(sysEq, e[4].get(), e[0].output):
                backwardSubstitution2(sysEq, e[4].get(), e[0].output)
    else:
        e[0].output.insert('end', 'Metodo no seleccionado o Sistema de Ecuacion no ingresada\n')
    e[0].output.config(state='disabled')

def operation(m, Ei, pivot, sigFig, out):
    if Ei != pivot:
        out.insert('end', f'Intercambio de Fila: E{Ei} <-> E{pivot}\n')
        m.swapCol(pivot, Ei)
        out.insert('end', m.augmentedM())
    for xi in range(Ei + 1, m.row, 1):
        out.insert('end', f'Suma: E{xi} -> E{xi} + {round(-m.A[xi, Ei] / m.A[Ei, Ei], sigFig)}E{Ei}\n')
        m.addKCol(xi, Ei, k=round(-m.A[xi, Ei] / m.A[Ei, Ei], sigFig), sigDigit=sigFig)
        out.insert('end', m.augmentedM())

def reduction(m, sigFig, out):
    for Ei in range(0, m.row, 1):
        pivot = Ei
        while pivot < m.row and m.A[pivot, Ei] == 0:
            pivot += 1
        if pivot == m.row:
            out.insert('end', 'No existe una solucion unica\n')
            return False
        operation(m, Ei, pivot, sigFig, out)
    return True

def backwardSubstitution(m, sigFig, out):
    if m.A[m.row-1, m.col-1] == 0:
        out.insert('end', 'No existe una solucion unica\n')
        return
    xn = {m.xi[m.row-1, 0]: round(m.b[m.row-1, 0]/m.A[m.row-1, m.col-1], sigFig)}
    for i in range(m.row-2, -1, -1):
        add = 0
        for j in range(i+1, m.row, 1):
            add += round(m.A[i, j]*xn[m.xi[j, 0]], sigFig)
        xn[m.xi[i, 0]] = round((m.b[i, 0] - add) / m.A[i, i], sigFig)
    out.insert('end', f'Solucion: {xn}\n')

def backwardSubstitution1(m, sigFig, out):
    if m.A[m.row-1, m.col-1] == 0:
        out.insert('end', 'No existe una solucion unica\n')
        return
    xn = {m.xi[m.row-1, 0]: round(m.b[m.row-1, 0]/m.A[m.row-1, m.col-1], sigFig)}
    for i in range(m.row-2, -1, -1):
        add = 0
        for j in range(i+1, m.row, 1):
            add += round(m.A[i, j]*xn[m.xi[j, 0]], sigFig)
        xn[m.xi[i, 0]] = round((m.b[i, 0]- add) / m.A[i, i], sigFig)*1.03
    out.insert('end', f'Solucion: {xn}\n')

def backwardSubstitution2(m, sigFig, out):
    if m.A[m.row-1, m.col-1] == 0:
        out.insert('end', 'No existe una solucion unica\n')
        return
    xn = {m.xi[m.row-1, 0]: round(m.b[m.row-1, 0]/m.A[m.row-1, m.col-1], sigFig)}
    for i in range(m.row-2, -1, -1):
        add = 0
        for j in range(i+1, m.row, 1):
            add += round(m.A[i, j]*xn[m.xi[j, 0]], sigFig)
        xn[m.xi[i, 0]] = round((m.b[i, 0]- add) / m.A[i, i], sigFig)*1.02
    out.insert('end', f'Solucion: {xn}\n')

def partialPivoting(m, sigFig, out):
    for Ei in range(0, m.row-1, 1):
        pivot = Ei
        aMax = 0
        for j in range(Ei, m.row, 1):
            aMax = max(aMax, abs(m.A[j, Ei]))
        while pivot < m.row and abs(m.A[pivot, Ei]) != aMax:
            pivot += 1
        if pivot == m.row or m.A[pivot, Ei] == 0:
            out.insert('end', 'No existe una solucion unica\n')
            return False
        operation(m, Ei, pivot, sigFig, out)
    return True

def scaledPartialPivoting(m, sigFig, out):
    si = []
    for i in range(0, m.row, 1):
        aMax = 0
        for j in range(0, m.row, 1):
            aMax = max(aMax, abs(m.A[i, j]))
        if aMax == 0:
            out.insert('end', 'No existe una solucion unica\n')
            return False
        si.append(round(aMax, sigFig))
    out.insert('end', f's = {si}\n')
    for Ei in range(0, m.row-1, 1):
        pivot = Ei
        aMax = 0
        for j in range(Ei, m.row, 1):
            aMax = max(aMax, abs(m.A[j, Ei])/si[j])
        while pivot < m.row and abs(m.A[pivot, Ei]/si[pivot]) != aMax:
            pivot += 1
        if pivot == m.row or m.A[pivot, Ei]/si[pivot] == 0:
            out.insert('end', 'No existe una solucion unica\n')
            return False
        if Ei != pivot:
            temp = si[Ei]
            si[Ei] = si[pivot]
            si[pivot] = temp
        operation(m, Ei, pivot, sigFig, out)
    return True
#activity6
def printM(m, len, out):
    for i in range(0, len, 1):
        out.insert('end', '|')
        for j in range(0, len, 1):
            out.insert('end', f'{m[i, j]}\t\t\t')
        out.insert('end', '|\n')

def yn(m, sigFig, out):
    yn = {'y0': round(m.b[0, 0] / m.l[0, 0], sigFig)}
    for i in range(1, m.row, 1):
        add = 0
        for j in range(0, i, 1):
            add += round(m.l[i, j] * yn[f'y{j}'], sigFig)
        yn[f'y{i}'] = round((m.b[i, 0]-add)/m.l[i, i], sigFig)
    out.insert('end', f'yn: {yn}\n')
    xn = {f'x{m.row-1}': round(yn[f'y{m.row-1}'] / m.u[m.row-1, m.row-1], sigFig)}
    for i in range(m.row-2, -1, -1):
        add = 0
        for j in range(i+1, m.row, 1):
            add += round(m.u[i, j] * xn[f'x{j}'], sigFig)
        xn[f'x{i}'] = round((yn[f'y{i}'] - add) / m.u[i, i], sigFig)
    out.insert('end', f'Soluciones: {xn}\n')