import sympy as sym
import control as C
import matplotlib.pyplot as plt
import numpy as np

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

# define a, b, c, d
a = d_phi_F.subs([(F, 0), (x3, 0), (x4, 0)])
b = d_phi_x3.subs([(F, 0), (x3, 0), (x4, 0)])
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
G_x = (a*s**2 - a*d - b*c) / (s**4 - d*s**2)

num = [-c_value]
den = [1, 0, -d_value]
G_theta = C.TransferFunction(num, den)

# Forced response
t_span = np.linspace(0, 0.2, 20)
F_input = np.sin(100*t_span**2)
t_out, x3_out, _ = C.forced_response(G_theta, t_span, F_input)
plt.plot(t_out, x3_out)
plt.xlabel('Time (s)')
plt.ylabel('Angle (rad)')
plt.grid()
plt.show()
