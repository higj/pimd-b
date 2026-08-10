#pragma once

struct RpmdContext {
    bool enabled;          // Is RPMD mode active?
    int num_runs;          // Number of independent RPMD simulations to perform
    int nvt_discard_frac;  // Fraction of the NVT trajectory to discard before sampling starts
};