#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <functional>

struct Point {
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {}
    
    Point operator+(const Point& other) const {
        return Point(x + other.x, y + other.y);
    }
    
    Point operator-(const Point& other) const {
        return Point(x - other.x, y - other.y);
    }
    
    Point operator*(double scalar) const {
        return Point(x * scalar, y * scalar);
    }
    
    double dot(const Point& other) const {
        return x * other.x + y * other.y;
    }
    
    double norm() const {
        return std::sqrt(x * x + y * y);
    }
    
    Point normalize() const {
        double n = norm();
        if (n < 1e-10) return Point(0, 0);
        return Point(x / n, y / n);
    }
};

struct TangentNormalResult {
    Point point;
    Point tangent;
    Point normal;
    bool valid;
    
    TangentNormalResult() : valid(false) {}
    TangentNormalResult(const Point& p, const Point& t, const Point& n) 
        : point(p), tangent(t), normal(n), valid(true) {}
};

class ParametricCurve {
public:
    std::function<double(double)> x_func;
    std::function<double(double)> y_func;
    double t_min;
    double t_max;
    
    ParametricCurve(std::function<double(double)> x, 
                    std::function<double(double)> y,
                    double tmin, double tmax)
        : x_func(x), y_func(y), t_min(tmin), t_max(tmax) {}
    
    Point evaluate(double t) const {
        return Point(x_func(t), y_func(t));
    }
    
    TangentNormalResult computeTangentNormal(double t, double dt = 1e-5) const {
        double dx_dt = (x_func(t + dt) - x_func(t - dt)) / (2.0 * dt);
        double dy_dt = (y_func(t + dt) - y_func(t - dt)) / (2.0 * dt);

        Point tangent(dx_dt, dy_dt);
        double tangent_norm = tangent.norm();
        
        if (tangent_norm < 1e-10) {
            return TangentNormalResult();
        }
        
        Point tangent_unit = tangent.normalize();
        
        Point normal_unit(-tangent_unit.y, tangent_unit.x);
        
        Point point = evaluate(t);
        
        return TangentNormalResult(point, tangent_unit, normal_unit);
    }
    double computeCurvature(double t, double dt = 1e-5) const {
        double dx = (x_func(t + dt) - x_func(t - dt)) / (2 * dt);
        double dy = (y_func(t + dt) - y_func(t - dt)) / (2 * dt);
        double ddx = (x_func(t + dt) - 2*x_func(t) + x_func(t - dt)) / (dt * dt);
        double ddy = (y_func(t + dt) - 2*y_func(t) + y_func(t - dt)) / (dt * dt);
        
        double numerator = dx*ddy - dy*ddx;
        double denominator = std::pow(dx*dx + dy*dy, 1.5);
        
        if (std::abs(denominator) < 1e-10) return 0.0;
        return numerator / denominator;
    }

    Point computeFocalPoint(double t, double R_max = 10.0) const {
        TangentNormalResult tn = computeTangentNormal(t);
        if (!tn.valid) return Point(NAN, NAN);
        
        double kappa = computeCurvature(t);
        if (std::abs(kappa) < 1e-10) return Point(NAN, NAN);
        
        double R = 1.0 / kappa;
        if (std::abs(R) > R_max) return Point(NAN, NAN);  // skip points that are too far
        
        return tn.point + tn.normal * R;
    }
};

std::vector<Point> computeSymmetrySet(const ParametricCurve& curve, 
                                      int resolution = 400,
                                      double lambda_max = 50.0) {
    std::vector<Point> centers;
    
    double dt = (curve.t_max - curve.t_min) / resolution;
    std::vector<double> tvals;
    for (int i = 0; i < resolution; i++) {
        tvals.push_back(curve.t_min + i * dt);
    }
    
    double threshold = 5.0 * dt;
    
    for (size_t i = 0; i < tvals.size(); i++) {
        double t1 = tvals[i];
        TangentNormalResult result1 = curve.computeTangentNormal(t1);
        if (!result1.valid) continue;
        
        Point p1 = result1.point;
        Point T1 = result1.tangent;
        Point N1 = result1.normal;
        
        double prev_check1 = NAN;
        double prev_check2 = NAN;
        double prev_t2 = NAN;
        
        for (size_t j = i + 1; j < tvals.size(); j++) {
            double t2 = tvals[j];
            
            if (std::abs(t2 - t1) < threshold) continue;
            
            TangentNormalResult result2 = curve.computeTangentNormal(t2);
            if (!result2.valid) continue;
            
            Point p2 = result2.point;
            Point T2 = result2.tangent;
            Point N2 = result2.normal;
            
            // Skip if tangents are parallel or anti-parallel
            double T_diff = (T1 - T2).norm();
            double T_sum = (T1 + T2).norm();
            
            if (T_diff < 1e-5 || T_sum < 1e-5) {
                prev_check1 = NAN;
                prev_check2 = NAN;
                continue;
            }
            
            // Compute tangency conditions
            Point diff = p1 - p2;
            Point T_sum_vec = T1 + T2;
            Point T_diff_vec = T1 - T2;
            
            double eq_check1 = diff.dot(T_sum_vec);
            double eq_check2 = diff.dot(T_diff_vec);
            
            // Check for sign changes
            if (!std::isnan(prev_check1)) {
                // Check if eq_check1 changed sign
                if (prev_check1 * eq_check1 < 0) {
                    // Interpolate to find zero crossing
                    double alpha = std::abs(prev_check1) / (std::abs(prev_check1) + std::abs(eq_check1));
                    double t_interp = prev_t2 + alpha * (t2 - prev_t2);
                    
                    TangentNormalResult result_interp = curve.computeTangentNormal(t_interp);
                    if (result_interp.valid) {
                        Point p_interp = result_interp.point;
                        Point N_interp = result_interp.normal;
                        
                        // Use the sum formula for center calculation
                        Point Nsum = N1 + N_interp;
                        double denom = Nsum.dot(Nsum);
                        
                        if (std::abs(denom) > 1e-8) {
                            double lam = (p_interp - p1).dot(Nsum) / denom;
                            
                            if (std::abs(lam) < lambda_max) {
                                Point center = p1 + N1 * lam;
                                centers.push_back(center);
                            }
                        }
                    }
                }
                
                // Check if eq_check2 changed sign
                if (prev_check2 * eq_check2 < 0) {
                    // Interpolate to find zero crossing
                    double alpha = std::abs(prev_check2) / (std::abs(prev_check2) + std::abs(eq_check2));
                    double t_interp = prev_t2 + alpha * (t2 - prev_t2);
                    
                    TangentNormalResult result_interp = curve.computeTangentNormal(t_interp);
                    if (result_interp.valid) {
                        Point p_interp = result_interp.point;
                        Point N_interp = result_interp.normal;
                        
                        // Use the difference formula for center calculation
                        Point Nsub = N1 - N_interp;
                        double denom = Nsub.dot(Nsub);
                        
                        if (std::abs(denom) > 1e-8) {
                            double lam = (p_interp - p1).dot(Nsub) / denom;
                            
                            if (std::abs(lam) < lambda_max) {
                                Point center = p1 + N1 * lam;
                                centers.push_back(center);
                            }
                        }
                    }
                }
            }
            
            // Store current values for next iteration
            prev_check1 = eq_check1;
            prev_check2 = eq_check2;
            prev_t2 = t2;
        }
    }
    
    return centers;
}

std::vector<Point> computeFocalSet(const ParametricCurve& curve, int resolution = 400, double R_max = 10.0) {
    std::vector<Point> focalPoints;
    double dt = (curve.t_max - curve.t_min) / resolution;
    for (int i = 0; i <= resolution; i++) {
        double t = curve.t_min + i * dt;
        Point c = curve.computeFocalPoint(t, R_max);
        if (!std::isnan(c.x) && !std::isnan(c.y))
            focalPoints.push_back(c);
    }
    return focalPoints;
}

void saveToCsv(const std::string& filename, 
               const std::vector<Point>& curve_points,
               const std::vector<Point>& symmetry_points,
               const std::vector<Point>& focal_points) {
    std::ofstream file(filename);
    if (!file.is_open()) return;
    
    file << "type,x,y\n";
    for (const auto& p : curve_points) file << "curve," << p.x << "," << p.y << "\n";
    for (const auto& p : symmetry_points) file << "symmetry," << p.x << "," << p.y << "\n";
    for (const auto& p : focal_points) file << "focal," << p.x << "," << p.y << "\n";
    file.close();
}

int main() {
    const double PI = 3.14159265358979323846;
    
    // Define the ellipse curve
    double a = -0.55;
    double b = 0.6;
    // x[t_] := Cos[t] (1 + -0.55 Exp[-20 (t - 2)^2] + 0.6 Exp[-70 (t - 2.4)^2]);
    // y[t_] := 1.8 Sin[t] (1 + -0.55 Exp[-20 (t - 2)^2] + 0.6 Exp[-70 (t - 2.4)^2]);
    auto x_func = [a,b](double t) { return std::cos(t) * (1 + a * std::exp(-20 * pow((t - 2),2)) + b * std::exp(-70 * pow((t - 2.4),2))); };
    auto y_func = [a,b](double t) { return 1.8 * std::sin(t) * (1 + a * std::exp(-20 * pow((t - 2),2)) + b * std::exp(-70 * pow((t - 2.4),2))); };
    
    ParametricCurve curve(x_func, y_func, - PI, PI);
    
    // Compute symmetry set
    std::cout << "Computing symmetry set..." << std::endl;
    int resolution = 200000;
    double lambda_max = 10.0;
    std::vector<Point> symmetry_centers = computeSymmetrySet(curve, resolution, lambda_max);
    
    // Generate curve points for plotting
    std::vector<Point> curve_points;
    double dt = (curve.t_max - curve.t_min) / resolution;
    for (int i = 0; i <= resolution; i++) {
        double t = curve.t_min + i * dt;
        curve_points.push_back(curve.evaluate(t));
    }

    std::cout << "Computing focal set..." << std::endl;
    double R_max = 50000.0;
    std::vector<Point> focal_points = computeFocalSet(curve, resolution, R_max);
    
    // Save to CSV
    saveToCsv("monodromyfig4_v3.csv", curve_points, symmetry_centers, focal_points);
    std::cout << "Done! Found " << focal_points.size() << " focal set points." << std::endl;
    std::cout << "Done! Found " << symmetry_centers.size() << " symmetry set points." << std::endl;
    
    return 0;
}