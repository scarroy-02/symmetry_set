#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

struct Point {
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {}
    
    Point operator+(const Point& other) const { return Point(x + other.x, y + other.y); }
    Point operator-(const Point& other) const { return Point(x - other.x, y - other.y); }
    Point operator*(double scalar) const { return Point(x * scalar, y * scalar); }
    Point operator/(double scalar) const { return Point(x / scalar, y / scalar); }
    
    double dot(const Point& other) const { return x * other.x + y * other.y; }
    double norm() const { return std::sqrt(x * x + y * y); }
    
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

std::vector<Point> loadPointsFromCsv(const std::string& filename) {
    std::vector<Point> points;
    std::ifstream file(filename);
    std::string line;
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return points;
    }

    // Skip header
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;
        double x, y;

        std::vector<std::string> tokens;
        while(std::getline(ss, cell, ',')) {
            tokens.push_back(cell);
        }

        if(tokens.size() >= 2) {
            try {
                int x_idx = (tokens.size() == 3) ? 1 : 0;
                x = std::stod(tokens[x_idx]);
                y = std::stod(tokens[x_idx+1]);
                points.emplace_back(x, y);
            } catch (...) { continue; }
        }
    }
    return points;
}

class DiscreteCurve {
public:
    std::vector<Point> points;
    
    DiscreteCurve(const std::vector<Point>& pts) : points(pts) {}
    
    size_t size() const { return points.size(); }
    
    Point getPoint(int i) const {
        if (i < 0) i = 0;
        if (i >= points.size()) i = points.size() - 1;
        return points[i];
    }

    TangentNormalResult computeTangentNormal(int i) const {
        if (i <= 0 || i >= points.size() - 1) return TangentNormalResult();

        Point p_prev = points[i - 1];
        Point p_next = points[i + 1];
        Point p_curr = points[i];

        Point tangent = p_next - p_prev; 
        
        if (tangent.norm() < 1e-10) return TangentNormalResult();
        
        Point t_unit = tangent.normalize();
        Point n_unit(-t_unit.y, t_unit.x); 
        
        return TangentNormalResult(p_curr, t_unit, n_unit);
    }

    double computeCurvature(int i) const {
        if (i <= 0 || i >= points.size() - 1) return 0.0;

        Point p_prev = points[i - 1];
        Point p_curr = points[i];
        Point p_next = points[i + 1];

        Point v = (p_next - p_prev) * 0.5;

        Point a = (p_next - (p_curr * 2.0)) + p_prev;

        double numerator = v.x * a.y - v.y * a.x;
        double denom = std::pow(v.x * v.x + v.y * v.y, 1.5);

        if (std::abs(denom) < 1e-10) return 0.0;
        
        return numerator / denom;
    }

    Point computeFocalPoint(int i, double R_max = 10.0) const {
        TangentNormalResult tn = computeTangentNormal(i);
        if (!tn.valid) return Point(NAN, NAN);
        
        double kappa = computeCurvature(i);
        if (std::abs(kappa) < 1e-5) return Point(NAN, NAN); // Flat line
        
        double R = 1.0 / kappa;
        
        // infinite curvature or flat spots
        if (std::abs(R) > R_max) return Point(NAN, NAN); 
        
        return tn.point + tn.normal * R;
    }
};

std::vector<Point> computeSymmetrySet(const DiscreteCurve& curve, 
                                      double lambda_max = 50.0) {
    std::vector<Point> centers;
    int n = curve.size();

    double min_chord_dist = 0.1; 

    double radius_tolerance = 0.05; 

    for (int i = 1; i < n - 1; i += 1) {
        
        TangentNormalResult res1 = curve.computeTangentNormal(i);
        if (!res1.valid) continue;
        
        Point p1 = res1.point;
        Point T1 = res1.tangent;
        Point N1 = res1.normal;
        
        double prev_check1 = NAN;
        double prev_check2 = NAN;
        for (int j = i + 1; j < n - 1; j += 1) {
            
            // 1. Fast Distance Check
            Point p_temp = curve.getPoint(j);
            double dist_sq = (p1.x - p_temp.x)*(p1.x - p_temp.x) + (p1.y - p_temp.y)*(p1.y - p_temp.y);
            if (dist_sq < min_chord_dist * min_chord_dist) continue;

            TangentNormalResult res2 = curve.computeTangentNormal(j);
            if (!res2.valid) continue;
            
            Point p2 = res2.point;
            Point T2 = res2.tangent;
            Point N2 = res2.normal;
            Point diff = p1 - p2;

            Point T_sum_vec = T1 + T2;
            Point T_diff_vec = T1 - T2;

            double eq_check1 = diff.dot(T_sum_vec);
            double eq_check2 = diff.dot(T_diff_vec);
            
            if (!std::isnan(prev_check1) && (prev_check1 * eq_check1 < 0)) {
                Point Nsum = N1 + N2; 
                double denom = Nsum.dot(Nsum);
                if (std::abs(denom) > 1e-5) {
                    double lam = (p2 - p1).dot(Nsum) / denom;
                    if (std::abs(lam) < lambda_max) {
                        Point C = p1 + N1 * lam;
                        double r1 = std::abs(lam);
                        double r2 = (C - p2).norm();
                        if (std::abs(r1 - r2) < (r1 * radius_tolerance)) {
                            centers.push_back(C);
                        }
                    }
                }
            }
            
            if (!std::isnan(prev_check2) && (prev_check2 * eq_check2 < 0)) {
                Point Nsub = N1 - N2;
                double denom = Nsub.dot(Nsub);
                if (std::abs(denom) > 1e-5) {
                    double lam = (p2 - p1).dot(Nsub) / denom;
                    if (std::abs(lam) < lambda_max) {
                        Point C = p1 + N1 * lam;
                        double r1 = std::abs(lam);
                        double r2 = (C - p2).norm();
                        if (std::abs(r1 - r2) < (r1 * radius_tolerance)) {
                            centers.push_back(C);
                        }
                    }
                }
            }

            prev_check1 = eq_check1;
            prev_check2 = eq_check2;
        }
    }
    return centers;
}

std::vector<Point> computeFocalSet(const DiscreteCurve& curve, double R_max = 10.0) {
    std::vector<Point> focalPoints;
    for (int i = 0; i < curve.size(); i++) {
        Point c = curve.computeFocalPoint(i, R_max);
        if (!std::isnan(c.x) && !std::isnan(c.y))
            focalPoints.push_back(c);
    }
    return focalPoints;
}

void saveToCsv(const std::string& filename, 
               const std::vector<Point>& symmetry_points,
               const std::vector<Point>& focal_points) {
    std::ofstream file(filename);
    if (!file.is_open()) return;
    
    file << "type,x,y\n";
    for (const auto& p : symmetry_points) file << "symmetry," << p.x << "," << p.y << "\n";
    for (const auto& p : focal_points) file << "focal," << p.x << "," << p.y << "\n";
    file.close();
}

int main() {
    std::string input_file = "spiral_points.csv";
    std::cout << "Loading points from " << input_file << "..." << std::endl;
    
    std::vector<Point> points = loadPointsFromCsv(input_file);
    
    if (points.size() < 10) {
        std::cout << "Error: Not enough points loaded. Make sure the file exists." << std::endl;
        return 1;
    }
    
    std::cout << "Loaded " << points.size() << " points." << std::endl;
    
    DiscreteCurve curve(points);
    
    std::cout << "Computing focal set..." << std::endl;
    std::vector<Point> focal_points = computeFocalSet(curve, 50.0); // R_max
    
    std::cout << "Computing symmetry set ..." << std::endl;
    std::vector<Point> symmetry_centers = computeSymmetrySet(curve, 10); // lambda_max
    
    // Save output
    std::string output_file = "analysis_output.csv";
    saveToCsv(output_file, symmetry_centers, focal_points);
    
    std::cout << "Done!" << std::endl;
    std::cout << "Focal Set Points: " << focal_points.size() << std::endl;
    std::cout << "Symmetry Set Points: " << symmetry_centers.size() << std::endl;
    std::cout << "Saved to " << output_file << std::endl;
    
    return 0;
}