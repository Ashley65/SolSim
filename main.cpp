#include <iostream>
#include <string>
#include <vector>
#include <cmath>

const double G = 6.67430e-11; // Gravitational constant
const double metersToAstronomicalUnits = 1.495978707e11; // 1 meter is equal to this many astronomical units
const double c = 299792458; // Speed of light in meters per second


struct CelestialBody {
    std::string name;
    double mass;       // in kilograms
    double radius;     // in meters
    std::vector<double> position;  // {x, y, z} in meters
    std::vector<double> velocity;  // {vx, vy, vz} in meters per second
    std::vector<double> acceleration; // {ax, ay, az} in meters per second squared
};



void convertPositionToAU(CelestialBody& body) {
    for (int i = 0; i < body.position.size(); ++i) {
        body.position[i] /= metersToAstronomicalUnits;
    }
}



void computeGravitationalForce(CelestialBody &body1, CelestialBody &body2) {
    std::vector<double> direction(3);
    double distance = 0.0;

    // Calculate the distance between the bodies and the direction vector
    for (int i = 0; i < 3; ++i) {
        direction[i] = body2.position[i] - body1.position[i];
        distance += direction[i] * direction[i];
    }
    distance = std::sqrt(distance);

    // Calculate the force magnitude
    double forceMagnitude = G * body1.mass * body2.mass / (distance * distance);

    // Relativistic correction
    double relativisticCorrection = 1 +(3 * G * (body1.mass + body2.mass)) / (c * c * distance);


    // Normalize the direction vector and calculate acceleration
    for (int i = 0; i < 3; ++i) {
        double force = forceMagnitude * relativisticCorrection * direction[i] / distance;
        body1.acceleration[i] += force / body1.mass;
        body2.acceleration[i] -= force / body2.mass; // Equal and opposite force
    }
}

void updateBodies(std::vector<CelestialBody> &bodies, double timestep) {
    for (auto &body : bodies) {
        for (int i = 0; i < 3; ++i) {
            body.velocity[i] += body.acceleration[i] * timestep;
            body.position[i] += body.velocity[i] * timestep;
            body.acceleration[i] = 0; // Reset acceleration for the next timestep
        }
    }
}
void updateBodiesEuler(std::vector<CelestialBody> &bodies, double timestep) {
    for (auto &body : bodies) {
        for (int i = 0; i < 3; ++i) {
            body.velocity[i] += body.acceleration[i] * timestep;
            body.position[i] += body.velocity[i] * timestep;
            body.acceleration[i] = 0; // Reset acceleration for the next timestep
        }
    }
}

void updateBodiesVerlet(std::vector<CelestialBody> &bodies, double timestep) {
    static std::vector<std::vector<double>> previousPositions(bodies.size(), std::vector<double>(3));

    for (size_t i = 0; i < bodies.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            double newPosition = 2 * bodies[i].position[j] - previousPositions[i][j] + bodies[i].acceleration[j] * timestep * timestep;
            previousPositions[i][j] = bodies[i].position[j];
            bodies[i].position[j] = newPosition;
            bodies[i].velocity[j] = (newPosition - previousPositions[i][j]) / (2 * timestep);
            bodies[i].acceleration[j] = 0; // Reset acceleration for the next timestep
        }
    }
}

void updateBodiesVelocityVerlet(std::vector<CelestialBody> &bodies, double timestep) {
    for (auto &body : bodies) {
        for (int i = 0; i < 3; ++i) {
            body.position[i] += body.velocity[i] * timestep + 0.5 * body.acceleration[i] * timestep * timestep;
            body.velocity[i] += 0.5 * body.acceleration[i] * timestep;
        }
    }

    // After computing the new positions, compute new accelerations based on the updated positions
    for (size_t i = 0; i < bodies.size(); ++i) {
        for (size_t j = i + 1; j < bodies.size(); ++j) {
            computeGravitationalForce(bodies[i], bodies[j]);
        }
    }

    for (auto &body : bodies) {
        for (int i = 0; i < 3; ++i) {
            body.velocity[i] += 0.5 * body.acceleration[i] * timestep;
            body.acceleration[i] = 0; // Reset acceleration for the next timestep
        }
    }
}

void updateBodiesRK4(std::vector<CelestialBody> &bodies, double timestep) {
    struct State {
        std::vector<double> position;
        std::vector<double> velocity;
    };

    for (auto &body : bodies) {
        State k1, k2, k3, k4;
        k1.position = body.velocity;
        k1.velocity = body.acceleration;

        std::vector<double> originalPosition = body.position;
        std::vector<double> originalVelocity = body.velocity;

        for (int i = 0; i < 3; ++i) {
            body.position[i] = originalPosition[i] + 0.5 * k1.position[i] * timestep;
            body.velocity[i] = originalVelocity[i] + 0.5 * k1.velocity[i] * timestep;
        }

        for (size_t i = 0; i < bodies.size(); ++i) {
            for (size_t j = i + 1; j < bodies.size(); ++j) {
                computeGravitationalForce(bodies[i], bodies[j]);
            }
        }
        k2.position = body.velocity;
        k2.velocity = body.acceleration;

        for (int i = 0; i < 3; ++i) {
            body.position[i] = originalPosition[i] + 0.5 * k2.position[i] * timestep;
            body.velocity[i] = originalVelocity[i] + 0.5 * k2.velocity[i] * timestep;
        }

        for (size_t i = 0; i < bodies.size(); ++i) {
            for (size_t j = i + 1; j < bodies.size(); ++j) {
                computeGravitationalForce(bodies[i], bodies[j]);
            }
        }
        k3.position = body.velocity;
        k3.velocity = body.acceleration;

        for (int i = 0; i < 3; ++i) {
            body.position[i] = originalPosition[i] + k3.position[i] * timestep;
            body.velocity[i] = originalVelocity[i] + k3.velocity[i] * timestep;
        }

        for (size_t i = 0; i < bodies.size(); ++i) {
            for (size_t j = i + 1; j < bodies.size(); ++j) {
                computeGravitationalForce(bodies[i], bodies[j]);
            }
        }
        k4.position = body.velocity;
        k4.velocity = body.acceleration;

        for (int i = 0; i < 3; ++i) {
            body.position[i] = originalPosition[i] + (k1.position[i] + 2 * k2.position[i] + 2 * k3.position[i] + k4.position[i]) * timestep / 6.0;
            body.velocity[i] = originalVelocity[i] + (k1.velocity[i] + 2 * k2.velocity[i] + 2 * k3.velocity[i] + k4.velocity[i]) * timestep / 6.0;
        }
    }
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
    std::vector<CelestialBody> bodies = {sun, earth};

    const double timestep = 86400; // 1 day in seconds
    const int numSteps = 365 * 10; // Simulate for 10 years

    std::vector<double> energiesEuler, energiesVerlet, energiesVelocityVerlet, energiesRK4;

    // Simulate with Euler method
    for (int step = 0; step < numSteps; ++step) {
        for (size_t i = 0; i < bodies.size(); ++i) {
            for (size_t j = i + 1; j < bodies.size(); ++j) {
                computeGravitationalForce(bodies[i], bodies[j]);
            }
        }
        updateBodiesEuler(bodies, timestep);
        energiesEuler.push_back(computeTotalEnergy(bodies));
    }

    // Reset bodies
    bodies[0] = sun;
    bodies[1] = earth;

    // Simulate with Verlet method
    for (int step = 0; step < numSteps; ++step) {
        for (size_t i = 0; i < bodies.size(); ++i) {
            for (size_t j = i + 1; j < bodies.size(); ++j) {
                computeGravitationalForce(bodies[i], bodies[j]);
            }
        }
        updateBodiesVerlet(bodies, timestep);
        energiesVerlet.push_back(computeTotalEnergy(bodies));
    }


    // Reset bodies
    bodies[0] = sun;
    bodies[1] = earth;

    // Simulate with Velocity Verlet method
    for (int step = 0; step < numSteps; ++step) {
        for (size_t i = 0; i < bodies.size(); ++i) {
            for (size_t j = i + 1; j < bodies.size(); ++j) {
                computeGravitationalForce(bodies[i], bodies[j]);
            }
        }
        updateBodiesVelocityVerlet(bodies, timestep);
        energiesVelocityVerlet.push_back(computeTotalEnergy(bodies));
    }

    // Reset bodies
    bodies[0] = sun;
    bodies[1] = earth;

    // Simulate with RK4 method
    for (int step = 0; step < numSteps; ++step) {
        for (size_t i = 0; i < bodies.size(); ++i) {
            for (size_t j = i + 1; j < bodies.size(); ++j) {
                computeGravitationalForce(bodies[i], bodies[j]);
            }
        }
        updateBodiesRK4(bodies, timestep);
        energiesRK4.push_back(computeTotalEnergy(bodies));
    }

    // Output energy results for comparison
    std::cout << "Step\tEuler\t\tVerlet\t\tVelocity Verlet\tRK4\n";
    for (int step = 0; step < numSteps; ++step) {
        std::cout << step << "\t" << energiesEuler[step] << "\t" << energiesVerlet[step] << "\t"
                  << energiesVelocityVerlet[step] << "\t" << energiesRK4[step] << "\n";
    }

    return 0;
}

