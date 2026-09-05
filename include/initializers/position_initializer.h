#pragma once

#include "common.h"
#include "contexts/box_context.h"

#include <memory>

class PositionInitializer {
public:
    explicit PositionInitializer(const std::shared_ptr<VecArray>& coord, const BoxContext& box_ctx)
        : m_coord(coord), m_box_ctx(box_ctx), m_natoms(coord->len()) {}
    virtual ~PositionInitializer() = default;
    virtual void initialize() = 0;

protected:
    std::shared_ptr<VecArray> m_coord;
    BoxContext m_box_ctx;
    int m_natoms;
};