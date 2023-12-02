import numpy as np

def constant_accel_wrapper(accel):
    def constant_accel(t0, state, mu):
        v_vec = state[3:]
        norm_v = np.linalg.norm(v_vec)
        a_x = accel * v_vec[0] / norm_v
        a_y = accel * v_vec[1] / norm_v
        a_z = accel * v_vec[2] / norm_v
        
        return [0, 0, 0, a_x, a_y, a_z]

    return constant_accel


def atmospheric_drag_wrapper(R, CD, A_over_m, H0, rho0):
    def atmospheric_drag_exponential(t0, state, mu):

        r_vec = state[:3]
        v_vec = state[3:]
        r = np.linarg.norm(r_vec)
        v = np.linarg.norm(v_vec)

        B = CD * A_over_m
        rho = rho0 * np.exp(-(r - R) / H0)

        return -(1.0 / 2.0) * rho * B * v * v_vec
    
def J2_perturbation_wrapper(J2, R):
    def J2_perturbation(t0, state, mu):

        r = state[:3]
        r_norm = np.linalg.norm(r)

        factor = (3.0 / 2.0) * mu * J2 * (R**2) / (r_norm**5)

        a_x = r[0] * factor * (5.0 * r[2] ** 2 / r_norm**2 - 1)
        a_y = r[1] * factor * (5.0 * r[2] ** 2 / r_norm**2 - 1)
        a_z = r[2] * factor * (5.0 * r[2] ** 2 / r_norm**2 - 3)

        return np.array([0, 0, 0, a_x, a_y, a_z])
    return J2_perturbation

def J3_perturbation_wrapper(J3, R):
    def J3_perturbation(t0, state, mu):

        r = state[:3]
        r_norm = np.linalg.norm(r)

        factor = (1.0 / 2.0) * mu * J3 * (R**3) / (r_norm**5)

        cos_phi = r[2] / r_norm

        a_x = 5.0 * r[0] / r_norm * (7.0 * cos_phi**3 - 3.0 * cos_phi)
        a_y = 5.0 * r[1] / r_norm * (7.0 * cos_phi**3 - 3.0 * cos_phi)
        a_z = 3.0 * (35.0 / 3.0 * cos_phi**4 - 10.0 * cos_phi**2 + 1)
        
        return np.array([0, 0, 0, a_x, a_y, a_z]) * factor
    return J3_perturbation