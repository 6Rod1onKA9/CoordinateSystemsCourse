#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>
#include <iomanip>
#include <chrono> // For timing
#include <cstdlib> // For random number generation

#define PI 3.14159265358979323846

// Function prototypes
std::tuple<double, double, double> spherical_to_cartesian(double rho, double theta, double phi);
double cartesian_distance_3d(double x1, double y1, double z1, double x2, double y2, double z2);
double great_circle_distance(double rho1, double theta1, double phi1, double rho2, double theta2, double phi2);
double spherical_distance(double rho1, double theta1, double phi1, double rho2, double theta2, double phi2);
void generate_random_points(int n, std::vector<std::tuple<double, double, double>>& points_spherical);

// Function to convert spherical coordinates to Cartesian coordinates
std::tuple<double, double, double> spherical_to_cartesian(double rho, double theta, double phi) {
    double x = rho * sin(phi) * cos(theta);
    double y = rho * sin(phi) * sin(theta);
    double z = rho * cos(phi);
    return { x, y, z };
}

// Function to calculate the straight-line distance in Cartesian coordinates
double cartesian_distance_3d(double x1, double y1, double z1, double x2, double y2, double z2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2));
}

// Function to compute the distance between two points in spherical coordinates
double spherical_distance(double rho1, double theta1, double phi1, double rho2, double theta2, double phi2) {
    // Convert spherical to Cartesian
    double x1, y1, z1;
    std::tie(x1, y1, z1) = spherical_to_cartesian(rho1, theta1, phi1);

    double x2, y2, z2;
    std::tie(x2, y2, z2) = spherical_to_cartesian(rho2, theta2, phi2);

    // Calculate and return the distance
    return cartesian_distance_3d(x1, y1, z1, x2, y2, z2);
}

// Function to calculate the great-circle distance
double great_circle_distance(double rho1, double theta1, double phi1, double rho2, double theta2, double phi2) {
    double delta_theta = theta2 - theta1;
    double delta_phi = phi2 - phi1;

    // Using the haversine formula
    double a = sin(delta_phi / 2) * sin(delta_phi / 2) +
        cos(phi1) * cos(phi2) *
        sin(delta_theta / 2) * sin(delta_theta / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));

    return rho1 * c; // Assuming both points are on the same sphere
}

// Function to generate random spherical points
void generate_random_points(int n, std::vector<std::tuple<double, double, double>>& points_spherical) {
    for (int i = 0; i < n; ++i) {
        double rho = static_cast<double>(rand() % 100 + 1); // Random radius between 1 and 100
        double theta = static_cast<double>(rand()) / RAND_MAX * 2 * PI; // Random theta between 0 and 2π
        double phi = static_cast<double>(rand()) / RAND_MAX * PI; // Random phi between 0 and π
        points_spherical.emplace_back(rho, theta, phi);
    }
}

int main() {
    int num_points = 10000; // Number of random points
    std::vector<std::tuple<double, double, double>> points_spherical;
    generate_random_points(num_points, points_spherical);

    // Benchmarking distance calculations
    auto start = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < points_spherical.size() - 1; ++i) {
        double rho1, theta1, phi1;
        std::tie(rho1, theta1, phi1) = points_spherical[i];

        double rho2, theta2, phi2;
        std::tie(rho2, theta2, phi2) = points_spherical[i + 1];

        double dis = spherical_distance(rho1, theta1, phi1, rho2, theta2, phi2);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Time taken for spherical distance calculation: " << elapsed.count() << " seconds\n";

    start = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < points_spherical.size() - 1; ++i) {
        double rho1, theta1, phi1;
        std::tie(rho1, theta1, phi1) = points_spherical[i];

        double rho2, theta2, phi2;
        std::tie(rho2, theta2, phi2) = points_spherical[i + 1];

        double dis = great_circle_distance(rho1, theta1, phi1, rho2, theta2, phi2);
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Time taken for great-circle distance calculation: " << elapsed.count() << " seconds\n";

    return 0;
}
