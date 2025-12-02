#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <set>

// --- GUDHI INCLUDES ---
#include <gudhi/Matrix.h>
#include <gudhi/persistence_matrix_options.h>

using namespace Gudhi;
using namespace Gudhi::persistence_matrix;

// --- 1. CONFIGURATION ---

// We use a Vector-based column for robustness
struct VineyardOptions : Default_options<Column_types::VECTOR, true> {
    static const bool has_vine_update = true;
    static const bool has_column_pairings = true;
    static const bool is_of_boundary_type = false; 
    
    // CRITICAL: Allows us to swap columns by index
    static const Column_indexation_types column_indexation_type = Column_indexation_types::POSITION;
    
    static const bool has_map_column_container = false; 
};

using VineyardMatrix = Matrix<VineyardOptions>;

// --- 2. DATA STRUCTURES ---

struct Point {
    double x, y;
};

struct Simplex {
    int id;                 // Global ID
    int dim;                // 0 or 1
    double current_val;     // Filtration Value
    std::vector<int> boundary_ids; // For safety checks
};

// --- 3. CSV PARSERS (UPDATED FOR YOUR DATA) ---

std::vector<Point> read_potato_curve(const std::string& filename) {
    std::vector<Point> points;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << std::endl;
        return points;
    }
    
    std::string line;
    std::getline(file, line); // Skip Header: Type,X,Y
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string segment;
        std::vector<std::string> row;
        
        while (std::getline(ss, segment, ',')) {
            row.push_back(segment);
        }

        // Format: Type, X, Y
        // We need row[0] == "curve", X is row[1], Y is row[2]
        if (row.size() >= 3 && row[0] == "curve") { 
            try {
                double x = std::stod(row[1]);
                double y = std::stod(row[2]);
                points.push_back({x, y});
            } catch (...) { continue; }
        }
    }
    return points;
}

std::vector<Point> read_centers(const std::string& filename) {
    std::vector<Point> points;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << std::endl;
        return points;
    }

    std::string line;
    std::getline(file, line); // Skip Header: "X","Y"
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string segment;
        std::vector<std::string> row;
        
        while (std::getline(ss, segment, ',')) {
            row.push_back(segment);
        }

        // Format: X, Y
        if(row.size() >= 2) {
            try {
                double x = std::stod(row[0]);
                double y = std::stod(row[1]);
                points.push_back({x, y});
            } catch (...) { continue; }
        }
    }
    return points;
}

double dist(Point p1, Point p2) {
    return std::sqrt(std::pow(p1.x - p2.x, 2) + std::pow(p1.y - p2.y, 2));
}

// Check if simplex A is a face of simplex B
// Used to prevent swapping a Vertex past its Edge
bool is_face_of(const Simplex& A, const Simplex& B) {
    if (A.dim >= B.dim) return false; 
    for (int face_id : B.boundary_ids) {
        if (face_id == A.id) return true;
    }
    return false;
}

// --- 4. MAIN ---

int main() {
    auto curve_points = read_potato_curve("plot_data/potato_all.csv");
    auto centers = read_centers("plot_data/potato_circle_points.csv");
    int n = curve_points.size();

    std::cout << "Loaded " << n << " curve points and " << centers.size() << " center frames.\n";

    if (n < 2 || centers.empty()) {
        std::cerr << "Error: Insufficient data.\n";
        return 1;
    }

    // --- INITIALIZATION ---

    std::vector<Simplex> filtration;
    filtration.reserve(2 * n);

    // 1. Create Vertices (IDs 0 to n-1)
    for(int i=0; i<n; ++i) {
        filtration.push_back({i, 0, 0.0, {}});
    }
    
    // 2. Create Edges (IDs n to 2n-1)
    // Connects point i to i+1 (wrapping around)
    for(int i=0; i<n; ++i) {
        int u = i;
        int v = (i + 1) % n;
        std::vector<int> b = {u, v};
        filtration.push_back({n + i, 1, 0.0, b});
    }

    // 3. Compute Initial Values (Frame 0)
    Point start_center = centers[0];
    for(auto& s : filtration) {
        if (s.dim == 0) {
            s.current_val = dist(curve_points[s.id], start_center);
        } else {
            // Edge value is max of endpoints
            double d1 = dist(curve_points[s.boundary_ids[0]], start_center);
            double d2 = dist(curve_points[s.boundary_ids[1]], start_center);
            s.current_val = std::max(d1, d2);
        }
    }

    // 4. Initial Sort
    // Order: By Value, then by Dimension (Vertex < Edge)
    std::sort(filtration.begin(), filtration.end(), [](const Simplex& a, const Simplex& b) {
        if (std::abs(a.current_val - b.current_val) > 1e-9) 
            return a.current_val < b.current_val;
        return a.dim < b.dim; 
    });

    // 5. Build Map: Global ID -> Initial Sorted Rank
    std::vector<int> id_to_rank(2 * n);
    for(size_t i=0; i<filtration.size(); ++i) {
        id_to_rank[filtration[i].id] = i;
    }

    // 6. Build Matrix Boundaries
    // The matrix works on indices [0..2n-1]. 
    // We must translate Global IDs to these indices.
    std::vector<std::vector<unsigned>> boundaries;
    boundaries.reserve(2 * n);
    
    for(const auto& s : filtration) {
        std::vector<unsigned> b;
        if (s.dim == 1) {
            for(int global_id : s.boundary_ids) {
                b.push_back(id_to_rank[global_id]);
            }
            std::sort(b.begin(), b.end());
        }
        boundaries.push_back(b);
    }

    VineyardMatrix mat(boundaries);
    std::cout << "Matrix initialized successfully.\n";

    // --- VINEYARD LOOP ---

    // To verify output, we can print H0 birth/death
    // H0 features usually: 1 infinite (connected component), many noise (short life)
    
    for (size_t t = 1; t < centers.size(); ++t) {
        Point c = centers[t];

        // A. Update Filtration Values
        for(auto& s : filtration) {
            if (s.dim == 0) {
                s.current_val = dist(curve_points[s.id], c);
            } else {
                double d1 = dist(curve_points[s.boundary_ids[0]], c);
                double d2 = dist(curve_points[s.boundary_ids[1]], c);
                s.current_val = std::max(d1, d2);
            }
        }

        // B. Bubble Sort with Topology Guard
        bool swapped = true;
        while(swapped) {
            swapped = false;
            for(size_t i = 0; i < filtration.size() - 1; ++i) {
                
                Simplex& s1 = filtration[i];
                Simplex& s2 = filtration[i+1];

                // Check if they want to swap (s1 is heavier than s2)
                if (s1.current_val > s2.current_val) {
                    
                    // --- SAFETY CHECK ---
                    // If s1 is a Vertex of s2 (Edge), they CANNOT swap.
                    // This happens if numerical noise makes Vertex > Edge.
                    if (is_face_of(s1, s2)) {
                        // Correct the violation: Edge must be >= Vertex
                        s2.current_val = s1.current_val;
                        continue; // Skip swap
                    }

                    // --- PERFORM SWAP ---
                    mat.vine_swap(i);
                    std::swap(filtration[i], filtration[i+1]);
                    swapped = true;
                }
            }
        }

        // C. Print Barcode (Sample)
        if (t % 10 == 0) {
            std::cout << "--- Frame " << t << " ---\n";
            for (unsigned i = 0; i < mat.get_number_of_columns(); ++i) {
                auto& col = mat.get_column(i);
                if (col.is_paired()) {
                    unsigned pair_idx = col.get_paired_chain_index();
                    if (i < pair_idx) { // Birth
                        double birth = filtration[i].current_val;
                        double death = filtration[pair_idx].current_val;
                        if (death - birth > 0.001) { // Filter tiny noise
                             std::cout << "  Dim " << filtration[i].dim 
                                       << " [" << birth << ", " << death << ")\n";
                        }
                    }
                } else {
                    // Essential feature (Infinite life)
                    // Note: In Chain matrix, sometimes essential features are unpaired in specific ways,
                    // but generally for H0 they appear as unpaired columns or columns paired with -1 logic depending on version.
                    // For R=DV, loops are usually born at 'i' and never die.
                    double birth = filtration[i].current_val;
                    std::cout << "  Dim " << filtration[i].dim << " [" << birth << ", Inf)\n";
                }
            }
        }
    }

    std::cout << "Done.\n";
    return 0;
}