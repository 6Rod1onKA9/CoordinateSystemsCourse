Implementation of Coordinate System Transformations 

Objective: To familiarize oneself with different coordinate systems (Cartesian, polar, and spherical) and gain practical skills in transitioning between them. To determine the computational efficiency of distance calculations in these coordinate systems through benchmarking.

Solution:

Tasks
1) Transition Between Coordinate Systems: 
Two-Dimensional Space: Cartesian and Polar Coordinate Systems 
-	Define the coordinates of several points in the polar coordinate system. 
-	Convert these coordinates to the Cartesian coordinate system. 
-	Perform the reverse conversion from the Cartesian coordinate system to the polar coordinate system. 
-	Verify the correctness of the calculations by ensuring that the initial coordinates match the ones obtained after the reverse conversion. 

Three-Dimensional Space: Cartesian and Spherical Coordinate Systems 
-	Define the coordinates of several points in the spherical coordinate system. 
-	Convert these coordinates to the Cartesian coordinate system. 
-	Perform the reverse conversion from the Cartesian coordinate system to the spherical coordinate system. 
-	Verify the correctness of the calculations by ensuring that the initial coordinates match the ones obtained after the reverse conversion. 

Code listing and result:

2D

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

![image](https://github.com/user-attachments/assets/8c9447a5-8b7c-4e97-b5f4-30d6c787281c)

3D

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

![image](https://github.com/user-attachments/assets/f02a50c6-6144-4aff-89a2-a3fdb1da19b1)
