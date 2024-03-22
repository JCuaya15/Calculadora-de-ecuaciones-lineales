from sympy import sympify, symbols, Eq, diff, Matrix, eye, zeros, sqrt
from numpy import linspace
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from tkinter import Tk, Toplevel, Frame, PhotoImage, Label, Entry, Button, Text, Scrollbar
from tkinter.ttk import Combobox
from functools import partial

class Function():
    def __init__(self, name, expression):
        self.name = name
        self.x = symbols('x')
        self.exp = sympify(expression, convert_xor=True, rational=True, evaluate=False)

    def substitution(self, symbols):  ##symbols = {x: a}
        return self.exp.subs(symbols)

    def evaluate(self, symbols):
        return float(self.exp.subs(symbols))

    def getCodomain(self, domain, symbols):
        cod = []
        for x in domain:
            symbols['x'] = x
            cod.append(self.evaluate(symbols))
        return cod

    def derivate(self):
        return diff(self.exp, self.x)

    def equation(self, a):
        return Eq(self.exp, a)

    def __str__(self):
        return f'{self.name}(x) = {self.exp}'

class Graphic(Figure):
    def __init__(self):
        super().__init__(figsize=(3, 3), dpi=100)
        self.paint = self.add_subplot()
        self.paint.grid(True, linestyle='dotted', color='k')
        self.paint.spines['left'].set_position('center')
        self.paint.spines['left'].set_linestyle('dashed')
        self.paint.spines['bottom'].set_position('center')
        self.paint.spines['bottom'].set_linestyle('dashed')
        self.paint.spines['right'].set_visible(False)
        self.paint.spines['top'].set_visible(False)
        self.tight_layout()

    def reset(self):
        self.paint = self.add_subplot()
        self.paint.grid(True, linestyle='dotted', color='k')
        self.paint.spines['left'].set_position('center')
        self.paint.spines['left'].set_linestyle('dashed')
        self.paint.spines['bottom'].set_position('center')#('data', 0)='zero'
        self.paint.spines['bottom'].set_linestyle('dashed')
        self.paint.spines['right'].set_visible(False)
        self.paint.spines['top'].set_visible(False)
        self.tight_layout()

    def setPlot(self, title, domain, codomain):
        self.paint.plot(domain, codomain, label=title)

    def setFunction(self, fx, symbols, a, b, n):
        domain = linspace(a, b, n)
        codomain = fx.getCodomain(domain, symbols)
        self.paint.plot(domain, codomain, label=f'{fx.name}(x)')

    def setPoint(self, x1, x2, title):
        self.paint.plot(x1, x2, 'o', label=title)

    def domain(self, x1, x2):
        self.paint.plot([x1, x2], [0, 0], '|', label=[x1, x2])

class Gui(Tk):
    def __init__(self, title, title2, colour):
        super().__init__()
        self.config(bg=colour, bd=10, relief='ridge')
        self.iconbitmap('buap.ico')
        self.title(title)
        self.__coverPage = Frame(self)
        self.__coverPage.config(width=300, height=300, bg='white', bd=5, relief='sunken', cursor='xterm')
        self.__coverPage.pack()
        self.__buap = PhotoImage(file='Escudo.png')
        Label(self.__coverPage, text=title2, bg='white', font=('TIMES NEW ROMAN', 15)).place(x=30, y=25)
        Label(self.__coverPage, image=self.__buap, bg='white').place(x=50, y=70)
        self.resizable(False, False)

    def addActivity(self, t, method, colour):
        Button(self, text=t, bg=colour, bd=3, relief='raised', cursor='hand1', command=partial(method, self)).pack()

    def gui(self, entry, plot):#(entry=['title', 'colour1', 'colour2', 'colour1'])
        self.__secondGui = SecondGui(entry[0], entry[1])
        self.__secondGui.setInput(entry[2])
        if plot:
            self.__secondGui.setZonePlot()
        self.__secondGui.setOutput(entry[3])
        return self.__secondGui

class SecondGui(Toplevel):
    def __init__(self, title, colour):
        super().__init__()
        self.config(bg=colour, bd=10, relief='ridge')
        self.iconbitmap('buap.ico')
        self.title(title)

    def setInput(self, colour):
        self.inputColour = colour
        self.input = Frame(self)
        self.input.config(width=683, height=100, bg=colour, bd=5, relief='sunken')
        self.input.pack(side='top', fill='x')

    def setZonePlot(self):
        self.plot = Graphic()
        self.canvas = FigureCanvasTkAgg(self.plot, master=self)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self)
        self.toolbar.pack()
        self.canvas.get_tk_widget().pack(side='left', fill='both', expand=True)
        self.toolbar.update()

    def setOutput(self, colour):
        self.output = Text(self, width=1, height=20, bg=colour, bd=5, relief='sunken', wrap='none')
        self.__vertical = Scrollbar(self, command=self.output.yview)
        self.__vertical.pack(side='right', fill='y')
        self.__horizontal = Scrollbar(self, orient='horizontal', command=self.output.xview)
        self.__horizontal.pack(side='bottom', fill='x')
        self.output.pack(side='left', fill='both', expand=True)
        self.output.config(state='disabled', yscrollcommand=self.__vertical.set, xscrollcommand=self.__horizontal.set)

    def getList(self, txt, x1, y1, x2, y2, data, var):
        Label(self.input, text=txt, bg=self.inputColour, font=('TIMES NEW ROMAN', 11), cursor='xterm').place(x=x1, y=y1)
        list = Combobox(self.input, values=data, textvariable=var)
        list.place(x=x2, y=y2)
        return var

    def getEntry(self, txt, x1, y1, x2, y2, var):
        Label(self.input, text=txt, bg=self.inputColour, font=('TIMES NEW ROMAN', 11), cursor='xterm').place(x=x1, y=y1)
        entry = Entry(self.input, font=('TIMES NEW ROMAN', 11), textvariable=var)
        entry.place(x=x2, y=y2)
        return var

    def setButton(self, t, x1, y1, method, entry):
        Button(self, text=t, bd=3, relief='raised', cursor='hand1', command=partial(method, entry)).place(x=x1, y=y1)

class SystemEquations():
    def __init__(self, A, b):
        self.A = Matrix(A)
        self.b = Matrix(b)
        self.xi = Matrix(symbols(f'x:{len(b)}'))
        (self.row, self.col) = self.A.shape

    def swapCol(self, row1, row2):
        self.A = self.A.elementary_row_op('n<->m', row1=row1, row2=row2)
        self.b = self.b.elementary_row_op('n<->m', row1=row1, row2=row2)

    def mulKCol(self, row, k, sigDigit):
        self.A = self.A.elementary_row_op('n->kn', row=row, k=k)
        self.b = self.b.elementary_row_op('n->kn', row=row, k=k)
        for ai in range(0, self.row, 1):
            self.A[row, ai] = round(self.A[row, ai], sigDigit)
        self.b[row, 0] = round(self.b[row, 0], sigDigit)

    def addKCol(self, i, j, k, sigDigit):
        self.A = self.A.elementary_row_op(op="n->n+km", row=i, row2=j, k=k)
        self.b = self.b.elementary_row_op(op='n->n+km', row=i, row2=j, k=k)
        for ai in range(0, self.row, 1):
            self.A[i, ai] = round(self.A[i, ai], sigDigit)
        self.b[i, 0] = round(self.b[i, 0], sigDigit)

    def significantFigures(self, sigDigit):
        for i in range(0, self.row, 1):
            for j in range(0, self.row, 1):
                self.A[i, j] = round(self.A[i, j], sigDigit)
            self.b[i, 0] = round(self.b[i, 0], sigDigit)

    def augmentedM(self):
        m = ''
        for Ei in range(0, len(self.b), 1):
            m += '|'
            for j in range(0, len(self.b), 1):
                m += f'{self.A[Ei, j]}\t\t\t'
            m += f'|{self.b[Ei, 0]}\n'
        return m

    def lu(self, sigFig):
        self.l = eye(self.row)
        self.u = eye(self.row)
        for j in range(1, self.row, 1):
            self.u[0, j] = self.A[0, j] / self.l[0, 0]
            self.l[j, 0] = self.A[j, 0] / self.u[0, 0]
        for i in range(1, self.row - 1, 1):
            if 0 == self.l[i, i] * self.u[i, i]:
                return False
            summation = 0
            for k in range(0, i, 1):
                summation += round(self.l[i, k] * self.u[k, i], sigFig)
            self.u[i, i] = self.A[i, i] - summation
            for j in range(i + 1, self.row, 1):
                summation = 0
                for k in range(0, i, 1):
                    summation += round(self.l[i, k] * self.u[k, j], sigFig)
                self.u[i, j] = round(round(round(self.A[i, j] - summation, sigFig), sigFig) / self.l[i, i], sigFig)
                summation = 0
                for k in range(0, i, 1):
                    summation += round(self.l[j, k] * self.u[k, i], sigFig)
                self.l[j, i] = round(round(1 / self.u[i, i], sigFig) * round(self.A[j, i] - summation, sigFig), sigFig)
        summation = 0
        for k in range(0, self.row - 1, 1):
            summation += round(self.l[self.row - 1, k] * self.u[k, self.row - 1], sigFig)
        self.u[self.row - 1, self.row - 1] = self.A[self.row - 1, self.row - 1] - summation
        if self.l[self.row - 1, self.row - 1] * self.u[self.row - 1, self.row - 1] == 0:
            return False
        return True

    def LLt(self, sigFig):
        self.l = zeros(self.row)
        self.u = zeros(self.row)
        self.l[0, 0] = round(sqrt(self.A[0, 0]), sigFig)
        self.u[0, 0] = self.l[0, 0]
        for j in range(1, self.row, 1):
            self.l[j, 0] = round(self.A[j, 0]/self.l[0, 0], sigFig)
            self.u[0, j] = self.l[j, 0]
        for i in range(1, self.row-1, 1):
            summation = 0
            for k in range(0, i, 1):
                summation += round(self.l[i, k]**2, sigFig)
            self.l[i, i] = round(sqrt(round(self.A[i, i]-summation, sigFig)), sigFig)
            self.u[i, i] = self.l[i, i]
            print(f'add[{i}, {i}]={summation}')
            for j in range(i+1, self.row, 1):
                summation = 0
                for k in range(0, i, 1):
                    summation += round(self.l[j, k] * self.l[i, k], sigFig)
                self.l[j, i] = round(round(self.A[j, i] - summation, sigFig)/self.l[i, i], sigFig)
                self.u[i, j] = self.l[j, i]
        summation = 0
        for k in range(0, self.row-1, 1):
            summation += round(self.l[self.row-1, k] ** 2, sigFig)
        self.l[self.row-1, self.row-1] = round(sqrt(round(self.A[self.row-1, self.row-1]-summation, sigFig)), sigFig)
        self.u[self.row-1, self.row-1] = self.l[self.row-1, self.row-1]

    def __str__(self):
        sys = ''
        for Ei in range(0, len(self.b), 1):
            sys += f'E{Ei}: '
            for j in range(0, len(self.b), 1):
                sys += f'{self.A[Ei, j]}{self.xi[j, 0]}\t\t\t'
            sys += f'= {self.b[Ei, 0]}\n'
        return sys