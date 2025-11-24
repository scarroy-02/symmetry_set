#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>

// --- Configuration ---
const int NUM_FRAMES = 120;        // 4 seconds at 30fps
const int SS_RESOLUTION = 40000;    // Symmetry Set resolution (O(N^2) complexity)
const int CURVE_RESOLUTION = 40000; // Drawing resolution for the blue curve
const double TARGET_A = -0.55;
const double TARGET_B = 0.6;
// ---------------------

struct Point {
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {}
    
    Point operator+(const Point& other) const { return Point(x + other.x, y + other.y); }
    Point operator-(const Point& other) const { return Point(x - other.x, y - other.y); }
    Point operator*(double scalar) const { return Point(x * scalar, y * scalar); }
    double dot(const Point& other) const { return x * other.x + y * other.y; }
    double norm() const { return std::sqrt(x * x + y * y); }
    Point normalize() const {
        double n = norm();
        if (n < 1e-10) return Point(0, 0);
        return Point(x / n, y / n);
    }
};

struct TangentNormalResult {
    Point point, tangent, normal;
    bool valid;
    TangentNormalResult() : valid(false) {}
    TangentNormalResult(const Point& p, const Point& t, const Point& n) 
        : point(p), tangent(t), normal(n), valid(true) {}
};

class ParametricCurve {
public:
    std::function<double(double)> x_func;
    std::function<double(double)> y_func;
    double t_min, t_max;
    
    ParametricCurve(std::function<double(double)> x, std::function<double(double)> y, double tmin, double tmax)
        : x_func(x), y_func(y), t_min(tmin), t_max(tmax) {}
    
    Point evaluate(double t) const { return Point(x_func(t), y_func(t)); }
    
    TangentNormalResult computeTangentNormal(double t, double dt = 1e-5) const {
        double dx = (x_func(t + dt) - x_func(t - dt)) / (2.0 * dt);
        double dy = (y_func(t + dt) - y_func(t - dt)) / (2.0 * dt);
        Point tan(dx, dy);
        if (tan.norm() < 1e-10) return TangentNormalResult();
        Point tan_u = tan.normalize();
        return TangentNormalResult(evaluate(t), tan_u, Point(-tan_u.y, tan_u.x));
    }

    double computeCurvature(double t, double dt = 1e-5) const {
        double dx = (x_func(t + dt) - x_func(t - dt)) / (2 * dt);
        double dy = (y_func(t + dt) - y_func(t - dt)) / (2 * dt);
        double ddx = (x_func(t + dt) - 2*x_func(t) + x_func(t - dt)) / (dt * dt);
        double ddy = (y_func(t + dt) - 2*y_func(t) + y_func(t - dt)) / (dt * dt);
        double num = dx*ddy - dy*ddx;
        double den = std::pow(dx*dx + dy*dy, 1.5);
        if (std::abs(den) < 1e-10) return 0.0;
        return num / den;
    }

    Point computeFocalPoint(double t, double R_max = 10.0) const {
        TangentNormalResult tn = computeTangentNormal(t);
        if (!tn.valid) return Point(NAN, NAN);
        double k = computeCurvature(t);
        if (std::abs(k) < 1e-10) return Point(NAN, NAN);
        double R = 1.0 / k;
        if (std::abs(R) > R_max) return Point(NAN, NAN);
        return tn.point + tn.normal * R;
    }
};

std::vector<Point> computeSymmetrySet(const ParametricCurve& curve, int resolution, double lambda_max = 50.0) {
    std::vector<Point> centers;
    double dt = (curve.t_max - curve.t_min) / resolution;
    std::vector<double> tvals;
    for (int i = 0; i < resolution; i++) tvals.push_back(curve.t_min + i * dt);
    
    // Pre-compute tangent/normals for O(N^2) speedup
    std::vector<TangentNormalResult> cache(resolution);
    for(int i=0; i<resolution; i++) cache[i] = curve.computeTangentNormal(tvals[i]);

    double threshold = 5.0 * dt;

    for (size_t i = 0; i < tvals.size(); i++) {
        if (!cache[i].valid) continue;
        Point p1 = cache[i].point;
        Point T1 = cache[i].tangent;
        Point N1 = cache[i].normal;
        
        double prev_check1 = NAN;
        double prev_check2 = NAN;
        double prev_t2 = NAN;
        int prev_idx = -1;
        
        for (size_t j = i + 1; j < tvals.size(); j++) {
            double t2 = tvals[j];
            if (std::abs(t2 - tvals[i]) < threshold) continue;
            if (!cache[j].valid) continue;

            Point p2 = cache[j].point;
            Point T2 = cache[j].tangent;
            Point N2 = cache[j].normal;

            double T_sum_norm = (T1 + T2).norm(); // parallel check
            double T_diff_norm = (T1 - T2).norm(); // anti-parallel check
            
            if (T_diff_norm < 1e-5 || T_sum_norm < 1e-5) {
                prev_check1 = NAN; prev_check2 = NAN;
                continue;
            }

            Point diff = p1 - p2;
            double eq_check1 = diff.dot(T1 + T2);
            double eq_check2 = diff.dot(T1 - T2);

            if (!std::isnan(prev_check1)) {
                // Check 1 transition
                if (prev_check1 * eq_check1 < 0) {
                    double alpha = std::abs(prev_check1) / (std::abs(prev_check1) + std::abs(eq_check1));
                    
                    // Linear Interpolation of Normal and Point (approximate but fast)
                    Point N_prev = cache[prev_idx].normal;
                    Point p_prev = cache[prev_idx].point;
                    Point N_interp = (N_prev * (1.0 - alpha) + N2 * alpha).normalize();
                    Point p_interp = p_prev * (1.0 - alpha) + p2 * alpha;

                    Point Nsum = N1 + N_interp;
                    double denom = Nsum.dot(Nsum);
                    if (std::abs(denom) > 1e-8) {
                        double lam = (p_interp - p1).dot(Nsum) / denom;
                        if (std::abs(lam) < lambda_max) centers.push_back(p1 + N1 * lam);
                    }
                }
                // Check 2 transition
                if (prev_check2 * eq_check2 < 0) {
                    double alpha = std::abs(prev_check2) / (std::abs(prev_check2) + std::abs(eq_check2));
                    
                    Point N_prev = cache[prev_idx].normal;
                    Point p_prev = cache[prev_idx].point;
                    Point N_interp = (N_prev * (1.0 - alpha) + N2 * alpha).normalize();
                    Point p_interp = p_prev * (1.0 - alpha) + p2 * alpha;

                    Point Nsub = N1 - N_interp;
                    double denom = Nsub.dot(Nsub);
                    if (std::abs(denom) > 1e-8) {
                        double lam = (p_interp - p1).dot(Nsub) / denom;
                        if (std::abs(lam) < lambda_max) centers.push_back(p1 + N1 * lam);
                    }
                }
            }
            prev_check1 = eq_check1;
            prev_check2 = eq_check2;
            prev_idx = j;
        }
    }
    return centers;
}

int main() {
    const double PI = 3.14159265358979323846;
    std::ofstream file("animation_data.csv");
    file << "frame,type,x,y\n";
    file << std::fixed << std::setprecision(6);

    std::cout << "Starting generation (" << NUM_FRAMES << " frames)..." << std::endl;

    for (int frame = 0; frame <= NUM_FRAMES; ++frame) {
        double t_progress = (double)frame / NUM_FRAMES;
        double curr_a = TARGET_A * t_progress;
        double curr_b = TARGET_B * t_progress;

        // User's functions with interpolated a and b
        auto x_func = [curr_a, curr_b](double t) { 
            return std::cos(t) * (1 + curr_a * std::exp(-20 * pow((t - 2),2)) + curr_b * std::exp(-70 * pow((t - 2.4),2))); 
        };
        auto y_func = [curr_a, curr_b](double t) { 
            return 1.8 * std::sin(t) * (1 + curr_a * std::exp(-20 * pow((t - 2),2)) + curr_b * std::exp(-70 * pow((t - 2.4),2))); 
        };

        ParametricCurve curve(x_func, y_func, -PI, PI);

        // 1. Compute Curve
        double dt_draw = (2*PI) / CURVE_RESOLUTION;
        for(int i=0; i<=CURVE_RESOLUTION; ++i) {
            Point p = curve.evaluate(-PI + i*dt_draw);
            file << frame << ",curve," << p.x << "," << p.y << "\n";
        }

        // 2. Compute Symmetry Set
        auto ss_points = computeSymmetrySet(curve, SS_RESOLUTION, 50.0);
        for(const auto& p : ss_points) {
            file << frame << ",symmetry," << p.x << "," << p.y << "\n";
        }

        // 3. Compute Focal Set
        double dt_focal = (2*PI) / CURVE_RESOLUTION;
        for(int i=0; i<=CURVE_RESOLUTION; ++i) {
            Point p = curve.computeFocalPoint(-PI + i*dt_focal, 10.0); // R_max = 10
            if (!std::isnan(p.x)) {
                file << frame << ",focal," << p.x << "," << p.y << "\n";
            }
        }
        
        if (frame % 10 == 0) std::cout << "Frame " << frame << " done." << std::endl;
    }

    std::cout << "Data generation complete: animation_data.csv" << std::endl;
    return 0;
}