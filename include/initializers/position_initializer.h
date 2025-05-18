#pragma once

#include "common.h"

class PositionInitializer {
public:
    explicit PositionInitializer(const std::shared_ptr<dVec>& coord, double box_size)
        : m_coord(coord), m_box_size(box_size), m_natoms(coord->len()) {}
    virtual ~PositionInitializer() = default;
    virtual void initialize() = 0;

protected:
    std::shared_ptr<dVec> m_coord;
    double m_box_size;
    int m_natoms;
};