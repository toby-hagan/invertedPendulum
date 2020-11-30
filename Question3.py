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
'''
# DETERMINE THE IMPULSE RESPONSE OF THE SYSTEM (kick)
# If F(s) = 1, X3(s) = G_theta, x3(t) = inv_lap(X3)
x3_t = sym.inverse_laplace_transform(G_theta, s, t)
sym.pprint(x3_t)

# DETERMINE THE STEP RESPONSE OF THE SYSTEM (push)
# If F(s) = 1, X3(s) = G_theta, x3(t) = inv_lap(X3)
x3_step_t = sym.inverse_laplace_transform(G_theta/s, s, t)
sym.pprint(x3_step_t)

# DETERMINE THE FREQUENCY RESPONSE OF THE SYSTEM (shake)
# If F(s) = w^2/(s^2+w^2), X3(s) = G_theta*w^2/(s^2+w^2)
# x3(t) = inv_lap(X3)

w = sym.symbols('w', real=True, positive=True)
x3_freq_t = sym.inverse_laplace_transform(G_theta*w**2/(s**2 + w**2), s, t)
sym.pprint(x3_freq_t.simplify())

# Impulse response
t_imp, x3_imp = C.impulse_response(G_theta)
plt.plot(t_imp, x3_imp)
plt.xlabel('Time (s)')
plt.ylabel('Angle (rad)')
plt.grid()
plt.show()

# Step response
t_step, x3_step = C.step_response(G_theta)
plt.plot(t_step, x3_step)
plt.xlabel('Time (s)')
plt.ylabel('Angle (rad)')
plt.grid()
plt.show()
'''
# DETERMINE THE IMPULSE RESPONSE OF THE SYSTEM (kick)
# F(s) = 1, X1(s) = G_x, x1(t) = inv_lap(x1)
x1_t = sym.inverse_laplace_transform(G_x, s, t)
print(sym.latex(x1_t))
sym.pprint(G_x)
# unfortunately I seem to be getting a bad result here as it cannot be computed in the
# c.impulse_response()
'''
# DETERMINE THE STEP RESPONSE OF THE SYSTEM (push)
x3_step_t = sym.inverse_laplace_transform(G_x/s, s, t)
sym.pprint(x3_step_t)

# DETERMINE THE FREQUENCY RESPONSE OF THE SYSTEM (shake)
w = sym.symbols('w', real=True, positive=True)
x3_freq_t = sym.inverse_laplace_transform(G_x*w**2/(s**2 + w**2), s, t)
sym.pprint(x3_freq_t.simplify())
'''
t_imp, x3_imp = C.impulse_response(G_x)
plt.plot(t_imp, x3_imp)
plt.xlabel('Time (s)')
plt.ylabel('Angle (rad)')
plt.grid()
plt.show()
