#pragma once

struct BeadContext {
    int nbeads;
    int natoms;
    int this_bead;

    [[nodiscard]] bool isFirst() const { return this_bead == 0; }
    [[nodiscard]] bool isLast() const { return this_bead == nbeads - 1; }
    [[nodiscard]] bool isFirstOrLast() const { return isFirst() || isLast(); }
};