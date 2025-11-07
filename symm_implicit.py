
import numpy as np
import pandas as pd

# Define constants for the implicit curve
a = 1.025
b = 0.09
a2 = a**2

def F(x, y):
    """The implicit curve equation F(x, y) = 0."""
    return y**2 - 2*b*x*y - a2*x + a2*x**3

def Fx(x, y):
    """Partial derivative of F with respect to x."""
    return -2*b*y - a2 + 3*a2*x**2

def Fy(x, y):
    """Partial derivative of F with respect to y."""
    return 2*y - 2*b*x

def Fxx(x, y):
    """Second partial derivative of F with respect to x."""
    return 6*a2*x

def Fyy(x, y):
    """Second partial derivative of F with respect to y."""
    return 2.0

def Fxy(x, y):
    """Mixed partial derivative of F with respect to x and y."""
    return -2.0*b

def solve_y_for_x(x):
    """Solves the quadratic equation y^2 - (2bx)y - (a^2x - a^2x^3) = 0 for y."""
    discriminant = (2*b*x)**2 - 4 * 1 * (-a2*x + a2*x**3)
    if discriminant < 1e-9: # Use a small tolerance
        discriminant = 0
    sqrt_d = np.sqrt(discriminant)
    y1 = (2*b*x + sqrt_d) / 2.0
    y2 = (2*b*x - sqrt_d) / 2.0
    return y1, y2

def trace_curve_loop(resolution=800):
    """Generates an ordered list of points tracing the curve's loop."""
    # The loop of the curve exists for x in [0, x_root] where x_root is the positive root
    # of the discriminant. Find the root of -a^2*x^2 + b^2*x + a^2 = 0.
    roots = np.roots([-a2, b**2, a2])
    x_max = max(r for r in roots if np.isreal(r))

    x_vals_forward = np.linspace(0, x_max, resolution // 2, endpoint=True)
    
    # Trace the upper branch of the loop
    upper_branch = []
    for x in x_vals_forward:
        y1, _ = solve_y_for_x(x)
        upper_branch.append(np.array([x, y1]))
        
    # Trace the lower branch in reverse to form a continuous path
    lower_branch = []
    for x in reversed(x_vals_forward):
        _, y2 = solve_y_for_x(x)
        lower_branch.append(np.array([x, y2]))

    # Combine branches, removing duplicate endpoint
    return upper_branch + lower_branch[1:]

def compute_tangent_normal_implicit(p):
    """Computes unit tangent and normal vectors at a point p on the implicit curve."""
    grad = np.array([Fx(p[0], p[1]), Fy(p[0], p[1])])
    norm_grad = np.linalg.norm(grad)
    if norm_grad < 1e-10:
        return None, None
    normal = grad / norm_grad
    tangent = np.array([-normal[1], normal[0]])
    return tangent, normal

def compute_curvature_implicit(p):
    """Computes curvature at a point p on the implicit curve."""
    x, y = p[0], p[1]
    f_x, f_y = Fx(x, y), Fy(x, y)
    f_xx, f_yy, f_xy = Fxx(x, y), Fyy(x, y), Fxy(x, y)
    
    numerator = f_xx * f_y**2 - 2 * f_xy * f_x * f_y + f_yy * f_x**2
    denominator = (f_x**2 + f_y**2)**1.5
    
    if abs(denominator) < 1e-10:
        return 0.0
    return -numerator / denominator

def compute_focal_point_implicit(p, R_max=10.0):
    """Computes the center of curvature (focal point)."""
    _, N = compute_tangent_normal_implicit(p)
    if N is None:
        return np.array([np.nan, np.nan])
    kappa = compute_curvature_implicit(p)
    if abs(kappa) < 1e-10:
        return np.array([np.nan, np.nan])
    R = 1.0 / kappa
    if abs(R) > R_max:
        return np.array([np.nan, np.nan])
    return p + N * R

def compute_focal_set(points, R_max=3.0):
    """Computes the focal set (evolute) for a set of curve points."""
    focal_points = []
    for p in points:
        c = compute_focal_point_implicit(p, R_max)
        if not np.isnan(c).any():
            focal_points.append(c)
    return focal_points

def compute_symmetry_set(points, lambda_max=3.0):
    """Computes the symmetry set by finding centers of bitangent circles."""
    centers = []
    n_points = len(points)
    if n_points < 2:
        return []

    threshold_idx = int(0.05 * n_points) # Avoid checking points that are too close

    for i in range(n_points):
        p1 = points[i]
        T1, N1 = compute_tangent_normal_implicit(p1)
        if T1 is None: continue

        prev_check1, prev_check2 = np.nan, np.nan
        prev_p2 = None

        for j in range(i + threshold_idx, n_points):
            p2 = points[j]
            T2, N2 = compute_tangent_normal_implicit(p2)
            if T2 is None: continue

            if np.linalg.norm(T1 - T2) < 1e-5 or np.linalg.norm(T1 + T2) < 1e-5:
                prev_check1, prev_check2, prev_p2 = np.nan, np.nan, p2
                continue
                
            diff, T_sum, T_diff = p1 - p2, T1 + T2, T1 - T2
            eq_check1, eq_check2 = np.dot(diff, T_sum), np.dot(diff, T_diff)

            if prev_p2 is not None:
                for prev_check, check, formula in [(prev_check1, eq_check1, 'sum'), (prev_check2, eq_check2, 'diff')]:
                    if not np.isnan(prev_check) and prev_check * check < 0:
                        alpha = abs(prev_check) / (abs(prev_check) + abs(check))
                        p_interp = prev_p2 + alpha * (p2 - prev_p2)
                        _ , N_interp = compute_tangent_normal_implicit(p_interp)
                        
                        if N_interp is not None:
                            N_vec = N1 + N_interp if formula == 'sum' else N1 - N_interp
                            denom = np.dot(N_vec, N_vec)
                            if abs(denom) > 1e-8:
                                lam = np.dot(p_interp - p1, N_vec) / denom
                                if abs(lam) < lambda_max:
                                    centers.append(p1 + lam * N1)
            
            prev_check1, prev_check2, prev_p2 = eq_check1, eq_check2, p2
            
    return centers

def save_to_csv(filename, curve_points, symmetry_points, focal_points):
    """Saves the computed points to a CSV file."""
    curve_df = pd.DataFrame(curve_points, columns=['x', 'y'])
    curve_df['type'] = 'curve'
    symmetry_df = pd.DataFrame(symmetry_points, columns=['x', 'y'])
    symmetry_df['type'] = 'symmetry'
    focal_df = pd.DataFrame(focal_points, columns=['x', 'y'])
    focal_df['type'] = 'focal'
    df = pd.concat([curve_df, symmetry_df, focal_df], ignore_index=True)
    df[['type', 'x', 'y']].to_csv(filename, index=False)

# --- Main Execution ---
resolution = 3000
lambda_max = 15.0
R_max = 10.0

print("1. Tracing the implicit curve...")
curve_points = trace_curve_loop(resolution)

print("2. Computing focal set...")
focal_points = compute_focal_set(curve_points, R_max)

print("3. Computing symmetry set...")
symmetry_points = compute_symmetry_set(curve_points, lambda_max)

output_filename = "implicit_symmetry_focal_set.csv"
save_to_csv(output_filename, curve_points, symmetry_points, focal_points)

print(f"\nProcessing complete.")
print(f"- Found {len(curve_points)} points on the curve.")
print(f"- Found {len(focal_points)} points in the focal set.")
print(f"- Found {len(symmetry_points)} points in the symmetry set.")
print(f"- Data saved to {output_filename}")

if symmetry_points:
    print("\nFirst 5 symmetry set points:")
    for i in range(min(5, len(symmetry_points))):
        p = symmetry_points[i]
        print(f"  ({p[0]:.4f}, {p[1]:.4f})")
