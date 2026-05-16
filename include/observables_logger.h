#pragma once

#include <memory>
#include <vector>
#include <fstream>

class Observable;

class ObservablesLogger {
public:
    ObservablesLogger(
        const std::string& filename, 
        int this_bead,
        long frequency,
        const std::vector<std::shared_ptr<Observable>>& observables
    );
    ~ObservablesLogger();

    bool isLoggingStep(long step) const { return step % m_frequency == 0; }

    void log(long step);

private:
    std::ofstream m_file;
    int m_this_bead;
    long m_frequency;
    std::vector<std::shared_ptr<Observable>> m_observables;

    void writeTimeStep(long step);
    void writeObservables();
};