#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>
#include <utility>

#define PI 3.14159265358979323846

// Coordinate conversion functions
std::pair<double, double> polar_to_cartesian(double r, double theta) {
    double x = r * cos(theta);
    double y = r * sin(theta);
    return { x, y };
}

std::pair<double, double> cartesian_to_polar(double x, double y) {
    double r = sqrt(x * x + y * y);
    double theta = atan2(y, x);
    return { r, theta };
}

int main() {
    // Sample points for conversion verification
    std::vector<std::pair<double, double>> polar_points = { {5, PI / 4}, {10, PI / 3}, {7, PI / 6} };

    // 2D Conversion Verification
    for (const auto& p : polar_points) {
        double r = p.first;
        double theta = p.second;

        // Use std::pair to hold the returned values
        std::pair<double, double> cartesian_coords = polar_to_cartesian(r, theta);
        std::pair<double, double> polar_coords = cartesian_to_polar(cartesian_coords.first, cartesian_coords.second);

        // Access elements through .first and .second
        std::cout << "Original (r, theta): (" << r << ", " << theta << ") -> Cartesian (x, y): ("
            << cartesian_coords.first << ", " << cartesian_coords.second << ") -> Reconverted (r, theta): ("
            << polar_coords.first << ", " << polar_coords.second << ")\n";
    }

    return 0;
}
