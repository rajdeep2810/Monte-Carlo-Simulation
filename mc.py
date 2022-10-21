import math as m

h = 10**(-35); t = 0
# Counter parameters
R_c = 21                # Cathode radius (in mm)
R_a = 0.0125            # Anode radius (in mm)
V_op = 2400             # Operational Voltage (in V)

l = 0.00788              # mm
# general electron properties
m_e = 9.11*(10**(-31))  # kg
e = 1.602*(10**(-19))   # C

def d_dt(param,t):
    dt = h
    return (param(t + dt) - param(t - dt))/(2*dt)

def d2_dt2(param,t):
    dt = h
    return (param(t + 2*dt) + param(t - 2*dt) - 2*param(t))/((2*dt)**2)

def e_field(r):
    return (V_op)/(r*m.log(R_c/R_a))

def DE_1(m_e, r, theta, t):
    return m_e*(d2_dt2(r,t) + r*(d2_dt2(theta,t)))

def DE_2(m_e, r, theta, t):
    return m_e*(r*(d2_dt2(theta,t)) + 2*d_dt(r,t)*d_dt(theta,t))

def r_next(r, psi, cap_psi):
    return ((r + l*m.sin(psi*m.pi/180)*m.cos(cap_psi*m.pi/180))**2 + (l*m.sin(psi*m.pi/180)*m.sin(cap_psi*m.pi/180))**2)**0.5

def theta_next(theta, r, psi, cap_psi):
    if m.sin(cap_psi) < 0:
        return theta - (m.acos((r + l*m.sin(psi*m.pi/180)*m.cos(cap_psi*m.pi/180))/r_next(r, psi, cap_psi)))*(180/m.pi)
    else:
        return theta + (m.acos((r + l*m.sin(psi*m.pi/180)*m.cos(cap_psi*m.pi/180))/r_next(r, psi, cap_psi)))*(180/m.pi)

def z_next(z, psi):
    return z + l*m.cos(psi*m.pi/180)

def energy_next(r, energy, psi, cap_psi):
    return energy - ((V_op*m.log(r_next(r, psi, cap_psi)/r))/m.log(R_c/R_a))

def v_r(energy, psi, cap_psi):
    return ((2*energy*e/m_e)**0.5)*10**(-3)*(m.sin(psi*m.pi/180)*m.cos(cap_psi*m.pi/180))

def v_theta(energy, psi, cap_psi):
    return ((2*energy*e/m_e)**0.5)*10**(-3)*(m.sin(psi*m.pi/180)*m.sin(cap_psi*m.pi/180))

def v_z(energy, psi):
    return ((2*energy*e/m_e)**0.5)*10**(-3)*(m.cos(psi*m.pi/180))

def v_z_next(energy, psi):
    return v_z(energy, psi)

def v_theta_next(r, energy, psi, cap_psi):
    return v_theta(energy, psi, cap_psi)*(r/r_next(r, psi, cap_psi))

def v_r_next(r, energy, psi, cap_psi):
    if m.cos(cap_psi) < 0:
        return -((2*energy_next(r, energy, psi, cap_psi)/m_e) - ((v_theta_next(r, energy, psi, cap_psi))**2 + (v_z_next(energy, psi))**2))**0.5
    else:
        return ((2*energy_next(r, energy, psi, cap_psi)/m_e) - ((v_theta_next(r, energy, psi, cap_psi))**2 + (v_z_next(energy, psi))**2))**0.5

def cap_psi_next(r, energy, psi, cap_psi):
    if v_theta_next(r, energy, psi, cap_psi) < 0:
        return -(m.acos((v_r_next(r, energy, psi, cap_psi)/(((v_r_next(r, energy, psi, cap_psi))**2 + (v_theta_next(r, energy, psi, cap_psi))**2)**0.5))))*(180/m.pi)
    else:
        return (m.acos((v_r_next(r, energy, psi, cap_psi)/(((v_r_next(r, energy, psi, cap_psi))**2 + (v_theta_next(r, energy, psi, cap_psi))**2)**0.5))))*(180/m.pi)

def psi_next(r, energy, psi, cap_psi):
    return (m.acos((v_z_next(energy, psi)/((2*energy_next(energy, r, psi, cap_psi))/m_e)**0.5)))*(180/m.pi)