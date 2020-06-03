from matplotlib import MatplotlibDeprecationWarning
from mpmath import *
from np import *
from sympy import *
import matplotlib.pyplot as plt
import  numpy
import numpy as np
def get_centroid (poles):
    poles_sum = 0
    for i in range(len(poles)):
        poles_sum = poles_sum + poles[i].real
    centroid = poles_sum / len(poles)
    return centroid

def get_asymptotic_angles(max_q):
    asymptotic_angles = []
    i=0
    for i in range (max_q+1):
        asymptotic_angles.append(((2*i+1)*180)/4)
    return  asymptotic_angles


def get_real_pples(poles):
    real_poles =[]
    for i in range (len(poles)):
        if (np.iscomplex(poles[i]) == false):
            real_poles.append(poles[i])
    return  real_poles


def possible_root_locus(real_poles):
    s,e =(0,0)
    for i in range(1,len(real_poles)):
        if (i%2!=0):
            s,e= (real_poles[i],real_poles[i-1])
    return s,e


def specific_root_locus (ans,possible):
    for i in range(len(ans)):
        if (ans.args[i]<possible[1] and ans.args[i]>possible[0]):
            return ans.args[i]

    return 0


def get_equ(poles):
    s = symbols('s')
    function=1
    for i in range (len(poles)):
        function= function * (s-poles[i])
    function = expand(function)
    return function

def get_Diff(function):
    s = symbols('s')
    dif = diff(function, s)
    dif = simplify(dif)
    return dif

def routh(equ ):
    s =symbols('s')
    coef = Poly(equ, s)
    co_list =coef.coeffs()
    routh_arr = []
    for i in range (2):
        c=i
        row = []
        for j in range(2):
            row.append(co_list[c])
            c= c+2
        routh_arr.append(row)

    k = symbols('k')
    routh_arr[0].append(k)
    routh_arr[1].append(0)

    for i in range(2, 4):
        row = []
        for j in range(2):
            cell = ((routh_arr[i - 1][0]) * (routh_arr[i - 2][j + 1]) - (routh_arr[i - 2][0]) * (routh_arr[i - 1][j + 1])) / (routh_arr[i - 1][0])
            row.append(cell)
        row.append(0)
        routh_arr.append(row)

    K = solveset(routh_arr[len(routh_arr)-1][0], k, domain=S.Reals).args[0]
    equ = routh_arr[len(routh_arr)-2][0] * s**2
    equ = equ+ routh_arr[len(routh_arr)-2][1]
    cross_with_imag = solveset(equ.subs(k, K))
    return cross_with_imag

def departure_angles(poles):
    #to obtain theta1
    x= -poles[2].real-poles[0]
    y= poles[2].imag
    theta1 = 180 - atan2(y,x)*180/3.14

    # to obtain theta2
    x = poles[2].real - poles[1]
    x=abs(x)
    y = poles[2].imag
    theta2 = 180 - atan2(y, x) * 180 / 3.14

    theta3 = 90
    dep_angles = 180-(theta1+theta2+theta3)
    return dep_angles

def plot_poles(poles):
    for i in range (len(poles)):
        plt.scatter(poles[i].real, poles[i].imag, color="red", marker="o",linewidth = '5')

def plot_break_point ():
    x = np.linspace(0, -25, 10)
    y = np.linspace(0, 0, 10)
    plt.plot(x, y, color='black', linewidth='4')
    plt.scatter(breakPoint, 0, color="green", marker="o", linewidth='5')


def draw_Asym(asymptotic_angles):
    l = 50
    for i in range (len(asymptotic_angles)):
        y = l * sin(radians(asymptotic_angles[i]))
        x = l * cos(radians(asymptotic_angles[i]))
        plt.plot([centroid, centroid+x],[0,y],"--", color = 'black')
    plt.scatter(centroid, 0, color="blue", marker="o", linewidth='5')

def arc ():
    x = linspace(-80, -50, 50)
    y = numpy.sqrt((x + 31.25) ** 2 - 2 - (15.9 ** 2))
    plt.plot(x, y,color = 'blue')
    x = linspace (-80,-50,50)
    y = -1 * numpy.sqrt((x+31.25)**2-2-(15.9**2))
    plt.plot(x,y,color = 'blue')

def draw_dep_angles(dep_angles):
    l = 50
    y = l * sin(radians(int(dep_angles)))
    x = l * cos(radians(int(dep_angles)))
    plt.plot([-78, -x-78],[y,10], color = 'black')
    x = l * sin(radians(int(dep_angles+20)))
    y = l * cos(radians(int(dep_angles+20)))
    plt.plot([-50, -x-50], [-10, y], color='black')


def plotCrossWithImagPoints():
        plt.scatter(0, -22.8035085019828, color="green", marker="o", linewidth='5')
        plt.scatter(0, 22.8035085019828, color="green", marker="o", linewidth='5')

def plot_curve():
    y = linspace(-40, 40, 50)
    b = 22.0997
    parabola = 1 + (y / b) ** 2
    x = b * (numpy.sqrt(parabola)) + (-31.25)
    plt.plot(x, y,color = 'green')

try:
    plt.axes().spines['bottom'].set_position(('data',0))
    plt.axes().spines['left'].set_position(('data',0))
except:
    pass
poles = [0,-25, complex(-50,10), complex(-50,-10)]
centroid = get_centroid (poles)      #step1
max_q = len(poles)-1       #No. of poles - No. of zeros - 1
asymptotic_angles = get_asymptotic_angles(max_q)
real_poles =get_real_pples(poles)
equation = get_equ(poles)
diff = get_Diff(equation)
s = Symbol('s')
ans = solveset(diff,s,domain=S.Reals)
possible = possible_root_locus(real_poles)
breakPoint = specific_root_locus (ans,possible)
cross_with_imag = routh(equation)
dep_angles = departure_angles(poles)
draw_Asym(asymptotic_angles)
plot_break_point ()
plot_poles(poles)
plotCrossWithImagPoints()
plot_curve()
arc ()
plt.show()
print("Centroid = ", centroid)
print("Asymptotic_angles = ",asymptotic_angles)
print("equation= ",equation)
print("Differentiation of the equation = ", diff)
print("BreakPoint = ", breakPoint)
print("Cross_with_imaginary_axis= ",cross_with_imag)
print("Departure_angles = ", dep_angles)