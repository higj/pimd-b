#pragma once

#include <memory>
#include <vector>
#include <fstream>

class Observable;

class ObservablesLogger {
public:
    ObservablesLogger(const std::string& filename, int bead, const std::vector<std::shared_ptr<Observable>>& observables);
    ~ObservablesLogger();

    void log(int step);

private:
    std::ofstream m_file;
    int m_bead;
    std::vector<std::shared_ptr<Observable>> m_observables;
};