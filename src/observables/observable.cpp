#include "observables.h"
#include "simulation.h"

Observable::Observable(int out_freq, const std::string& out_unit) : m_out_freq(out_freq), m_out_unit(out_unit)
{
}

void Observable::initializeLabel(const std::string& label)
{
    quantities.insert({ label, 0.0 });
}

void Observable::initialize(const std::vector<std::string>& labels)
{
    for (const std::string& label : labels) {
        initializeLabel(label);
    }
}

void Observable::resetValues() {
    for (auto it = quantities.begin(); it != quantities.end(); ++it) {
        it.value() = 0.0;
    }
}