#include "core/random_generators.h"

RandomGenerators::RandomGenerators(unsigned int seed)
    : m_rng(seed), m_unit_dist(0.0, 1.0), m_mars_gen(seed) {
}

double RandomGenerators::uniform() {
    return m_unit_dist(m_rng);
}

double RandomGenerators::uniform(double a, double b) {
    std::uniform_real_distribution<double> dist(a, b);
    return dist(m_rng);
}

double RandomGenerators::normal(double mean, double stddev) {
    std::normal_distribution<double> dist(mean, stddev);
    return dist(m_rng);
}

double RandomGenerators::gaussian() {
    /// TODO: Later change this to the std::normal_distribution
    return m_mars_gen.gaussian();  // LAMMPS pimd/langevin random numbers
}

int RandomGenerators::integer(int min, int max) {
    std::uniform_int_distribution<int> dist(min, max);
    return dist(m_rng);
}

std::mt19937& RandomGenerators::engine() {
    return m_rng;
}
