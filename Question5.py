import sympy as sym
import control as C
import matplotlib.pyplot as plt
import numpy as np
import math

# step 0. introduce the necessary symbols
M, m, g, x1, x2, x3, x4, F, ell = sym.symbols('M, m, g, x1, x2, x3, x4, F, ell')
phi = 4*m*ell*x4**2*sym.sin(x3) + 4*F - 3*m*g*sym.sin(x3)*sym.cos(x3)
phi /= 4*(M+m) - 3*m*sym.cos(x3)**2

# declaration of psi function
psi = -3*(m*ell*x4**2*sym.sin(x3)*sym.cos(x3) + F*sym.cos(x3) - (M + m)*g*sym.sin(x3))
psi /= (4*(M + m) - 3*m*sym.cos(x3)**2)*ell

# step 1. differentiate phi wrt F, x3 and x4
d_phi_F = phi.diff(F)
d_phi_x3 = phi.diff(x3)
d_phi_x4 = phi.diff(x4)

# psi differentiations
d_psi_F = psi.diff(F)
d_psi_x3 = psi.diff(x3)
d_psi_x4 = psi.diff(x4)

# step 2. substitute F=0, x3=0, x4=0 into the
#         above derivatives
a = d_phi_F.subs([(F, 0), (x3, 0), (x4, 0)])
b = d_phi_x3.subs([(F, 0), (x3, 0), (x4, 0)])
d_phi_x4_eq = d_phi_x4.subs([(F, 0), (x3, 0), (x4, 0)])

d_psi_F_eq = d_psi_F.subs([(F, 0), (x3, 0), (x4, 0)])
d_psi_x3_eq = d_psi_x3.subs([(F, 0), (x3, 0), (x4, 0)])
d_psi_x4_eq = d_psi_x4.subs([(F, 0), (x3, 0), (x4, 0)])

# sym.pprint(d_psi_F.subs([(F, 0), (x3, 0), (x4, 0)]))
# print(sym.latex(d_psi_x3_eq))
# sym.pprint(a)

c = 3/ell/(4*M+m)
d = 3*(M+m)*g/ell/(4*M+m)

# M = 0.3, m = 0.1, g = 9.81, ell = 0.35


def substitute(z):
    substitutions = [(M, 0.3), (m, 0.1), (g, 9.81), (ell, 0.35)]
    return float(z.subs(substitutions))  # this will return <class 'sympy.core.numbers.Float'>,
    # but we want <class 'float'>. We can check with print(type(a_value))


a_value = substitute(a)
b_value = substitute(b)
c_value = substitute(c)
d_value = substitute(d)

s, t = sym.symbols('s, t')
a, b, c, d = sym.symbols('a, b, c, d', real=True, positive=True)  # clarify that numbers are real and positive,
# for clean output
G_theta = -c / (s**2 - d)

num = [-c_value]
den = [1, 0, -d_value]
G_theta = C.TransferFunction(num, den)
# print(G_theta)


# PID controller below
Kp = 10
Ki = 0.1
Kd = 2

G_c = -C.TransferFunction([Kd, Kp, Ki], [1, 0])
G_d = C.feedback(G_theta, G_c)

t_imp, x3_imp = C.impulse_response(G_d)
plt.plot(t_imp, (x3_imp*180)/np.pi)
plt.xlabel('Time (ms)')
plt.ylabel('Angle (deg)')
plt.grid()
plt.show()
