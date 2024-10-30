#include <iostream>
#include <cmath>        // Required for mathematical functions
#include <vector>      // For std::vector
#include <tuple>       // For std::tuple
#include <iomanip>     // For std::setprecision

#define PI 3.14159265358979323846

// Function to convert spherical coordinates to Cartesian coordinates
std::tuple<double, double, double> spherical_to_cartesian(double rho, double theta, double phi) {
    double x = rho * sin(phi) * cos(theta);
    double y = rho * sin(phi) * sin(theta);
    double z = rho * cos(phi);
    return { x, y, z };
}

// Function to calculate the straight-line distance in Cartesian coordinates
double cartesian_distance(double x1, double y1, double z1, double x2, double y2, double z2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2));
}

// Function to compute the distance between two points in spherical coordinates
double spherical_distance(double rho1, double theta1, double phi1, double rho2, double theta2, double phi2) {
    // Convert spherical to Cartesian
    double x1, y1, z1;  // Declare Cartesian coordinates for first point
    std::tie(x1, y1, z1) = spherical_to_cartesian(rho1, theta1, phi1);

    double x2, y2, z2;  // Declare Cartesian coordinates for second point
    std::tie(x2, y2, z2) = spherical_to_cartesian(rho2, theta2, phi2);

    // Calculate and return the distance
    return cartesian_distance(x1, y1, z1, x2, y2, z2);
}

int main() {
    // Example points in spherical coordinates
    double rho1 = 5.0, theta1 = PI / 4, phi1 = PI / 3; // Point 1
    double rho2 = 10.0, theta2 = PI / 3, phi2 = PI / 6; // Point 2

    // Calculate distance between the two points in spherical coordinates
    double distance = spherical_distance(rho1, theta1, phi1, rho2, theta2, phi2);

    // Output the result
    std::cout << std::fixed << std::setprecision(4); // Format output to 4 decimal places
    std::cout << "Distance between points: " << distance << "\n";

    return 0;
}

![image](https://github.com/user-attachments/assets/f0a2609d-ac2a-46a4-b66b-5419e93e3091)
