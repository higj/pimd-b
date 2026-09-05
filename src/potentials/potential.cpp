#include "potentials/potential.h"

Potential::Potential(bool is_free) : m_is_free(is_free), tail_correction(0.0) {}