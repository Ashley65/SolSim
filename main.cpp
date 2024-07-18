#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <thread>
#include <mutex>
#include <omp.h>



const double G = 6.67430e-11; // Gravitational constant
const double AU = 1.495978707e11; // 1 meter is equal to this many astronomical units
const double c = 299792458; // Speed of light in meters per second


struct CelestialBody {
    std::string name;
    double mass;       // in kilograms
    double radius;     // in meters
    std::vector<double> position;  // {x, y, z} in meters
    std::vector<double> velocity;  // {vx, vy, vz} in meters per second
    std::vector<double> acceleration; // {ax, ay, az} in meters per second squared
};

std::mutex mtx;

void computeGravitationalForce(CelestialBody &body1, CelestialBody &body2){
    std::vector<double> direction(3);
    double distance = 0.0;

    // Calculate the distance between the bodies and the direction vector

    for (int i = 0; i <3; ++i){
        direction[i] = body2.position[i] - body1.position[i];
        distance += direction[i] * direction[i];
    }
    distance = std::sqrt(distance);

    // Calculate the force magnitude
    double forceMagnitude = G * body1.mass * body2.mass / (distance * distance);

    // Relativistic correction
    double relativisticCorrection = 1 + (3 * G * (body1.mass + body2.mass)) / (c * c * distance);

    // Normalize the direction vector and calculate acceleration
    for (int i = 0; i < 3; ++i){
        double force = forceMagnitude * relativisticCorrection * direction[i] / distance;
        mtx.lock();
        body1.acceleration[i] += force / body1.mass;
        body2.acceleration[i] -= force / body2.mass; // Equal and opposite force
        mtx.unlock();
    }
}

void updateBodiesRK4(std::vector<CelestialBody> &bodies, double timestep) {
    struct State {
        std::vector<double> position;
        std::vector<double> velocity;
    };

    std::vector<std::thread> threads;

    for (size_t i = 0; i < bodies.size(); ++i) {
        threads.push_back(std::thread([&, i]() {
            CelestialBody &body = bodies[i];
            State k1, k2, k3, k4;
            k1.position = body.velocity;
            k1.velocity = body.acceleration;

            std::vector<double> originalPosition = body.position;
            std::vector<double> originalVelocity = body.velocity;

            for (int j = 0; j < 3; ++j) {
                body.position[j] = originalPosition[j] + 0.5 * k1.position[j] * timestep;
                body.velocity[j] = originalVelocity[j] + 0.5 * k1.velocity[j] * timestep;
            }

            // Reset accelerations to zero before calculating new forces
            std::fill(body.acceleration.begin(), body.acceleration.end(), 0.0);

            for (size_t j = 0; j < bodies.size(); ++j) {
                if (i != j) {
                    computeGravitationalForce(body, bodies[j]);
                }
            }
            k2.position = body.velocity;
            k2.velocity = body.acceleration;

            for (int j = 0; j < 3; ++j) {
                body.position[j] = originalPosition[j] + 0.5 * k2.position[j] * timestep;
                body.velocity[j] = originalVelocity[j] + 0.5 * k2.velocity[j] * timestep;
            }

            // Reset accelerations to zero before calculating new forces
            std::fill(body.acceleration.begin(), body.acceleration.end(), 0.0);

            for (size_t j = 0; j < bodies.size(); ++j) {
                if (i != j) {
                    computeGravitationalForce(body, bodies[j]);
                }
            }
            k3.position = body.velocity;
            k3.velocity = body.acceleration;

            for (int j = 0; j < 3; ++j) {
                body.position[j] = originalPosition[j] + k3.position[j] * timestep;
                body.velocity[j] = originalVelocity[j] + k3.velocity[j] * timestep;
            }

            // Reset accelerations to zero before calculating new forces
            std::fill(body.acceleration.begin(), body.acceleration.end(), 0.0);

            for (size_t j = 0; j < bodies.size(); ++j) {
                if (i != j) {
                    computeGravitationalForce(body, bodies[j]);
                }
            }
            k4.position = body.velocity;
            k4.velocity = body.acceleration;

            for (int j = 0; j < 3; ++j) {
                body.position[j] = originalPosition[j] + (k1.position[j] + 2 * k2.position[j] + 2 * k3.position[j] + k4.position[j]) * timestep / 6.0;
                body.velocity[j] = originalVelocity[j] + (k1.velocity[j] + 2 * k2.velocity[j] + 2 * k3.velocity[j] + k4.velocity[j]) * timestep / 6.0;
            }
        }));
    }

    for (auto &t : threads) {
        if (t.joinable()) {
            t.join();
        }
    }
}

void convertToAU(std::vector<CelestialBody> &bodies) {
    std::vector<std::thread> threads;

    for (size_t i = 0; i < bodies.size(); ++i) {
        threads.push_back(std::thread([&, i]() {
            for (int j = 0; j < 3; ++j) {
                bodies[i].position[j] /= AU;
            }
        }));
    }

    for (auto &t : threads) {
        if (t.joinable()) {
            t.join();
        }
    }
}

void writeDataToFile(const std::vector<CelestialBody> &bodies,std::ofstream &file, int day) {
    file << "Day " << day << "\n";
    for (const auto &body : bodies) {
        file << body.name << " position: ("
             << std::fixed << std::setprecision(6) << body.position[0] << " AU, "
             << body.position[1] << " AU, " << body.position[2] << " AU) ";
        file << "velocity: ("
             << body.velocity[0] << " m/s, " << body.velocity[1] << " m/s, " << body.velocity[2] << " m/s) ";
        file << "acceleration: ("
             << body.acceleration[0] << " m/s^2, " << body.acceleration[1] << " m/s^2, " << body.acceleration[2] << " m/s^2) ";
        file << "\n";
    }
    file << "\n";
}



// Creating Asteroid Belt

std::vector<CelestialBody> CreateAsteroidBelt(int numAsteroids){

    std::vector<CelestialBody> asteroids;
    double minDistance = 2.1 * AU;
    double maxDistance = 3.3 * AU;

    srand(static_cast<unsigned>(time(0)));

    for (int i = 0; i < numAsteroids; ++i) {
        double distance = minDistance + static_cast<double>(rand()) / RAND_MAX * (maxDistance - minDistance);
        double angle = static_cast<double>(rand()) / RAND_MAX * 2 * M_PI;
        double speed = std::sqrt(G * 1.989e30 / distance); // Circular orbit speed

        CelestialBody asteroid = {
                "Asteroid_" + std::to_string(i),
                1.0e15, // Approximate mass of a small asteroid in kg
                500.0, // Approximate radius of a small asteroid in meters
                {distance * std::cos(angle), distance * std::sin(angle), 0}, // Position in meters
                {-speed * std::sin(angle), speed * std::cos(angle), 0}, // Velocity in m/s
                {0, 0, 0} // Initial acceleration
        };

        asteroids.push_back(asteroid);
    }

    return asteroids;

}





double computeTotalEnergy(const std::vector<CelestialBody> &bodies) {
    double totalEnergy = 0.0;

    // Kinetic energy
    for (const auto &body : bodies) {
        double kineticEnergy = 0.5 * body.mass * (body.velocity[0] * body.velocity[0] +
                                                  body.velocity[1] * body.velocity[1] +
                                                  body.velocity[2] * body.velocity[2]);
        totalEnergy += kineticEnergy;
    }

    // Potential energy
    for (size_t i = 0; i < bodies.size(); ++i) {
        for (size_t j = i + 1; j < bodies.size(); ++j) {
            double distance = std::sqrt((bodies[i].position[0] - bodies[j].position[0]) * (bodies[i].position[0] - bodies[j].position[0]) +
                                        (bodies[i].position[1] - bodies[j].position[1]) * (bodies[i].position[1] - bodies[j].position[1]) +
                                        (bodies[i].position[2] - bodies[j].position[2]) * (bodies[i].position[2] - bodies[j].position[2]));
            double potentialEnergy = -G * bodies[i].mass * bodies[j].mass / distance;
            totalEnergy += potentialEnergy;
        }
    }

    return totalEnergy;
}





int main() {

    CelestialBody sun = {"Sun", 1.989e30, 6.96342e8, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    CelestialBody earth = {"Earth", 5.972e24, 6.371e6, {1.496e11, 0, 0}, {0, 29783, 0}, {0, 0, 0}};
    CelestialBody mars = {"Mars", 6.4171e23, 3.3895e6, {2.279e11, 0, 0}, {0, 24007, 0}, {0, 0, 0}};
    CelestialBody jupiter = {"Jupiter", 1.8982e27, 6.9911e7, {7.785e11, 0, 0}, {0, 13070, 0}, {0, 0, 0}};

    // Earth's Moon: Luna
    CelestialBody luna = {"Luna", 7.342e22, 1.737e6, {1.496e11 + 3.844e8, 0, 0}, {0, 29783 + 1022, 0}, {0, 0, 0}};

    // Jupiter's Moons
    CelestialBody io = {"Io", 8.9319e22, 1.8216e6, {7.785e11 + 4.217e8, 0, 0}, {0, 13070 + 17325, 0}, {0, 0, 0}};
    CelestialBody europa = {"Europa", 4.7998e22, 1.5608e6, {7.785e11 + 6.711e8, 0, 0}, {0, 13070 + 13740, 0}, {0, 0, 0}};
    CelestialBody ganymede = {"Ganymede", 1.4819e23, 2.6341e6, {7.785e11 + 1.0704e9, 0, 0}, {0, 13070 + 10870, 0}, {0, 0, 0}};
    CelestialBody callisto = {"Callisto", 1.0759e23, 2.4103e6, {7.785e11 + 1.8827e9, 0, 0}, {0, 13070 + 8204, 0}, {0, 0, 0}};

    std::vector<CelestialBody> bodies = {sun, earth, mars,
                                         jupiter, luna, io, europa, ganymede, callisto};

    // Create an asteroid belt

    std::vector<CelestialBody> asteroids = CreateAsteroidBelt(1000);
    bodies.insert(bodies.end(), asteroids.begin(), asteroids.end());

    const double timestep = 86400; // 1 day in seconds
    const int numSteps = 365 * 1; // Simulate for 10 years

    std::ofstream outputFile("simulation_results.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not open the file for writing.\n";
        return 1;
    }

    for (int step = 0; step < numSteps; ++step) {
        // Reset accelerations to zero before calculating new forces
        std::vector<std::thread> resetThreads;
        for (size_t i = 0; i < bodies.size(); ++i) {
            resetThreads.push_back(std::thread([&, i]() {
                std::fill(bodies[i].acceleration.begin(), bodies[i].acceleration.end(), 0.0);
            }));
        }
        for (auto &t : resetThreads) {
            if (t.joinable()) {
                t.join();
            }
        }

        // Compute gravitational forces
        std::vector<std::thread> forceThreads;
        for (size_t i = 0; i < bodies.size(); ++i) {
            forceThreads.push_back(std::thread([&, i]() {
                for (size_t j = i + 1; j < bodies.size(); ++j) {
                    computeGravitationalForce(bodies[i], bodies[j]);
                }
            }));
        }
        for (auto &t : forceThreads) {
            if (t.joinable()) {
                t.join();
            }
        }

        // Update positions and velocities using RK4 method
        updateBodiesRK4(bodies, timestep);

        // Convert positions to AU
        convertToAU(bodies);

        // Write the positions of the bodies to the file
        writeDataToFile(bodies, outputFile, step + 1);

        // Output the positions of the bodies in AU (for example, the Earth's position)
        std::cout << "Day " << step + 1 << ": ";
        for (const auto &body : bodies) {
            std::cout << body.name << " position: ("
                      << std::fixed << std::setprecision(6) << body.position[0] << " AU, "
                      << body.position[1] << " AU, " << body.position[2] << " AU) ";
            std::cout << "velocity: ("
                      << body.velocity[0] << " m/s, " << body.velocity[1] << " m/s, " << body.velocity[2] << " m/s) ";

            std::cout << "\n";

        }
        std::cout << "\n";

        // Convert positions back to meters for the next calculation
        for (auto &body : bodies) {
            for (int i = 0; i < 3; ++i) {
                body.position[i] *= AU;
            }
        }
    }
    outputFile.close();
    return 0;
}

