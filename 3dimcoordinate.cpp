#include <iostream>
#include <cmath>         // For math functions like sin, cos, sqrt, atan2, acos
#include <vector>       // For std::vector
#include <tuple>        // For std::tuple
#include <utility>      // For std::pair
#include <iomanip>      // For std::setprecision

#define PI 3.14159265358979323846

// Coordinate conversion functions
std::tuple<double, double, double> spherical_to_cartesian(double rho, double theta, double phi) {
    double x = rho * sin(phi) * cos(theta);
    double y = rho * sin(phi) * sin(theta);
    double z = rho * cos(phi);
    return { x, y, z };
}

std::tuple<double, double, double> cartesian_to_spherical(double x, double y, double z) {
    double rho = sqrt(x * x + y * y + z * z);
    double theta = atan2(y, x);
    double phi = acos(z / rho);
    return { rho, theta, phi };
}

int main() {
    // Sample points for conversion verification
    std::vector<std::tuple<double, double, double>> spherical_points = {
        {5, PI / 4, PI / 3},
        {10, PI / 3, PI / 6},
        {7, PI / 6, PI / 4}
    };

    // 3D Conversion Verification
    for (const auto& p : spherical_points) {
        double rho, theta, phi;
        std::tie(rho, theta, phi) = p;  // Unpack spherical coordinates

        // Convert spherical to Cartesian
        double x, y, z;  // Declare Cartesian coordinates
        std::tie(x, y, z) = spherical_to_cartesian(rho, theta, phi);

        // Convert back from Cartesian to spherical
        double rho_new, theta_new, phi_new;  // Declare new spherical coordinates
        std::tie(rho_new, theta_new, phi_new) = cartesian_to_spherical(x, y, z);

        // Output the results
        std::cout << std::fixed << std::setprecision(4); // Format output to 4 decimal places
        std::cout << "Original (rho, theta, phi): (" << rho << ", " << theta << ", " << phi << ") -> Cartesian (x, y, z): ("
            << x << ", " << y << ", " << z << ") -> Reconverted (rho, theta, phi): ("
            << rho_new << ", " << theta_new << ", " << phi_new << ")\n";
    }

    return 0;
}

![image](https://github.com/user-attachments/assets/675e14d9-cd43-4b88-96f4-3b7bc61b0870)
