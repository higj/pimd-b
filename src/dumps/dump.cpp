#include "dumps/dump.h"

Dump::Dump(int out_freq, const std::string& out_unit) : m_out_freq(out_freq), m_out_unit(out_unit) {
}

Dump::~Dump() {
    if (m_out_file.is_open()) {
        m_out_file.close();
    }
}