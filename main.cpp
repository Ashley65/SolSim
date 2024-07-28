#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <thread>
#include <mutex>
#include <cstdlib>

// RKF45 Coefficients
const double a2 = 1.0/4.0, a3 = 3.0/8.0, a4 = 12.0/13.0, a5 = 1.0, a6 = 1.0/2.0;
const double b21 = 1.0/4.0;
const double b31 = 3.0/32.0, b32 = 9.0/32.0;
const double b41 = 1932.0/2197.0, b42 = -7200.0/2197.0, b43 = 7296.0/2197.0;
const double b51 = 439.0/216.0, b52 = -8.0, b53 = 3680.0/513.0, b54 = -845.0/4104.0;
const double b61 = -8.0/27.0, b62 = 2.0, b63 = -3544.0/2565.0, b64 = 1859.0/4104.0, b65 = -11.0/40.0;
const double r1 = 1.0/360.0, r3 = -128.0/4275.0, r4 = -2197.0/75240.0, r5 = 1.0/50.0, r6 = 2.0/55.0;
const double c1 = 16.0/135.0, c3 = 6656.0/12825.0, c4 = 28561.0/56430.0, c5 = -9.0/50.0, c6 = 2.0/55.0;
const double d1 = 1.0/360.0, d3 = -128.0/4275.0, d4 = -2197.0/75240.0, d5 = 1.0/50.0, d6 = 2.0/55.0;

// Adaptive Timestep
const double tolerance = 1.0e-6;
const double safetyFactor = 0.9;
const double dtMin = 1.0e-6;
const double dtMax = 1.0e-1;


// Physical constants
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



void computeGravitationalForce(CelestialBody &body1, CelestialBody &body2 ,std::mutex &mtx) {

    std::vector<double> force(3, 0.0);
    std::vector<double> r(3, 0.0);
    double distanceSquared = 0.0;

    for (int i = 0; i < 3; ++i){
        r[i] = body2.position[i] - body1.position[i];
        distanceSquared += r[i] * r[i];
    }

    double distance = std::sqrt(distanceSquared);
    double forceMagnitude = (G * body1.mass * body2.mass) / (distanceSquared * distance);

    for (int i = 0; i < 3; ++i){
        force[i] = forceMagnitude * r[i];
    }

    std::lock_guard<std::mutex> lock(mtx);
    for (int i = 0; i < 3; ++i){
        body1.acceleration[i] += force[i] / body1.mass;
        body2.acceleration[i] -= force[i] / body2.mass;
    }



//    std::vector<double> direction(3);
//    double distance = 0.0;
//
//    // Calculate the distance between the bodies and the direction vector
//
//    for (int i = 0; i <3; ++i){
//        direction[i] = body2.position[i] - body1.position[i];
//        distance += direction[i] * direction[i];
//    }
//    distance = std::sqrt(distance);
//
//    // Calculate the force magnitude
//    double forceMagnitude = G * body1.mass * body2.mass / (distance * distance);
//
//    // Relativistic correction
//    double relativisticCorrection = 1 + (3 * G * (body1.mass + body2.mass)) / (c * c * distance);
//
//    // Normalize the direction vector and calculate acceleration
//    for (int i = 0; i < 3; ++i){
//        double force = forceMagnitude * relativisticCorrection * direction[i] / distance;
//        mtx.lock();
//        body1.acceleration[i] += force / body1.mass;
//        body2.acceleration[i] -= force / body2.mass; // Equal and opposite force
//        mtx.unlock();
//    }
}

void updateBodiesRK4(std::vector<CelestialBody> &bodies, double timestep) {
    struct State {
        std::vector<double> position;
        std::vector<double> velocity;
    };
    std::mutex mtx;

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
                    computeGravitationalForce(body, bodies[j], mtx);
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
                    computeGravitationalForce(body, bodies[j], mtx);
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
                    computeGravitationalForce(body, bodies[j], mtx);
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

void rkf45Step(CelestialBody &body, double dt, std::vector<double> &k1, std::vector<double> &k2, std::vector<double> &k3, std::vector<double> &k4, std::vector<double> &k5, std::vector<double> &k6, std::vector<double> &yTemp) {

    for (int i = 0; i < 3; ++i){
        k1[i] = dt * body.velocity[i];
        k2[i] = dt * (body.velocity[i] + b21 * k1[i]);
        k3[i] = dt * (body.velocity[i] + b31 * k1[i] + b32 * k2[i]);
        k4[i] = dt * (body.velocity[i] + b41 * k1[i] + b42 * k2[i] + b43 * k3[i]);
        k5[i] = dt * (body.velocity[i] + b51 * k1[i] + b52 * k2[i] + b53 * k3[i] + b54 * k4[i]);
        k6[i] = dt * (body.velocity[i] + b61 * k1[i] + b62 * k2[i] + b63 * k3[i] + b64 * k4[i] + b65 * k5[i]);
        yTemp[i] = body.position[i] + c1 *  k1[i] + c3 * k3[i] + c4 * k4[i] + c5 * k5[i] + c6 * k6[i];
    }
}

double estimateError(const std::vector<double> &k1, std::vector<double> &k2, std::vector<double> &k3, std::vector<double> &k4, std::vector<double> &k5, std::vector<double> &k6){
    double error = 0.0;
    for (int i = 0; i <3; ++i){
        error += std::pow(d1 * k1[i] + d3 * k3[i] + d4 * k4[i] + d5 * k5[i] + d6 * k6[i], 2);
    }
    return std::sqrt(error);
}

void updateBodyRKF45(CelestialBody &body, double dt, const std::vector<double> &k1, const std::vector<double> &k2, const std::vector<double> &k3, const std::vector<double> &k4, const std::vector<double> &k5, const std::vector<double> &k6) {
    for (int i = 0; i < 3; ++i) {
        body.position[i] += 25.0/216.0 * k1[i] + 1408.0/2565.0 * k3[i] + 2197.0/4104.0 * k4[i] - k5[i] / 5.0;
        body.velocity[i] += r1 * k1[i] + r3 * k3[i] + r4 * k4[i] + r5 * k5[i] + r6 * k6[i];
    }
}




void convertToAU(std::vector<CelestialBody> &bodies) {
    for (auto &body : bodies) {
        for (int i = 0; i < 3; ++i) {
            body.position[i] /= AU;
        }
    }
}

void writeDataToFile(const std::vector<CelestialBody> &bodies,std::ofstream &file,int step) {
    for (const auto &body : bodies) {
        file << step << "," << body.name << ",";
        for (const auto &pos : body.position) {
            file << pos << ",";
        }
        file << "\n";
    }
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

void simulate(std::vector<CelestialBody> &bodies, double timestep, int numSteps, const std::string &outputFile){
    std::ofstream file(outputFile);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the file for writing.\n";
        return;
    }
    std::mutex mtx;

    for (int step = 0; step < numSteps; ++step) {
        // Reset accelerations to zero before calculating new forces
        std::vector<std::thread> resetThreads;
        for (auto &body : bodies) {
            resetThreads.emplace_back([&]() {
                std::fill(body.acceleration.begin(), body.acceleration.end(), 0.0);
            });
        }
        for (auto &t: resetThreads) {
            if (t.joinable()) {
                t.join();
            }
        }
        // Compute gravitational forces
        std::vector<std::thread> forceThreads;
        for (size_t i = 0; i < bodies.size(); ++i) {
            forceThreads.emplace_back([&, i]() {
                for (size_t j = i + 1; j < bodies.size(); ++j) {
                    computeGravitationalForce(bodies[i], bodies[j], mtx);
                }
            });
        }
        for (auto &t: forceThreads) {
            if (t.joinable()) {
                t.join();
            }
        }

        //RKF45 integration with adaptive timestep
        std::vector<std::thread> rkfThreads;
        for (auto &body: bodies) {
            rkfThreads.emplace_back([&body, &timestep, &mtx]() {
                std::vector<double> k1(3), k2(3), k3(3), k4(3), k5(3), k6(3), yTemp(3);
                double currentStep = timestep;
                double error = 0.0;

                do{
                    rkf45Step(body, currentStep, k1, k2, k3, k4, k5, k6, yTemp);
                    error = estimateError(k1, k2, k3, k4, k5, k6);

                    if (error > tolerance){
                        currentStep *= safetyFactor * std::pow(tolerance / error, 0.25);
                        if (currentStep < dtMin){
                            currentStep = dtMin;
                        }
                    }
                } while ( error > tolerance);
                for (int i= 0; i < 3; ++i) {
                    body.position[i] = yTemp[i];
                    body.velocity[i] += r1 * k1[i] + r3 * k3[i] + r4 * k4[i] + r5 * k5[i] + r6 * k6[i];
                }
                std::lock_guard<std::mutex> lock(mtx);
                timestep = currentStep * std::pow(tolerance / error, 0.2);
                if (timestep > dtMax) {
                    timestep = dtMax;
                }
            });
        }
        for (auto &t: rkfThreads) {
            if (t.joinable()) {
                t.join();
            }

        }

        convertToAU(bodies);

        // Write to fil
        writeDataToFile(bodies, file, step + 1);

        // Convert positions back to meters
        for (auto &body : bodies) {
            for (int i = 0; i < 3; ++i) {
                body.position[i] *= AU;
            }
        }

        file.flush();
    }
    file.close();

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


    const double timestep = 3600; // 1 hour in seconds
    const int numSteps = 365 * 24 * 1; // Simulate for 10 years

    simulate(bodies, timestep, numSteps, "simulation_results.txt");

    // Run the Python visualization script
    int result = system("python3 ../script/vis.py"); // Use "python" or "python3" depending on your system setup

    if (result != 0) {
        std::cerr << "Error: Could not run the Python script.\n";
    }

    return 0;
}

