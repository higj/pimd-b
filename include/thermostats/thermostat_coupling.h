#pragma once

#include "common.h"

#include <memory>

class NormalModes;

// Classes to support coupling of the thermostat to Cartesian coords or normal modes of distinguishable ring polymers
class Coupling {
public:
    explicit Coupling(const std::shared_ptr<dVec>& momenta);
    virtual ~Coupling() = default;

    virtual void mpiCommunication() = 0;
    virtual double getMomentumForCalc(const int ptcl_idx, const int axis) = 0;
    virtual double& getMomentumForUpdate(const int ptcl_idx, const int axis) = 0;
    virtual void updateCoupledMomenta() = 0;
protected:
    std::shared_ptr<dVec> m_momenta;
};

class CartesianCoupling final : public Coupling {
public:
    explicit CartesianCoupling(const std::shared_ptr<dVec>& momenta);

    void mpiCommunication() override;
    double getMomentumForCalc(const int ptcl_idx, const int axis) override;
    double& getMomentumForUpdate(const int ptcl_idx, const int axis) override;
    void updateCoupledMomenta() override;    
};

class NormalModesCoupling final : public Coupling {
public:
    NormalModesCoupling(const std::shared_ptr<dVec>& momenta, const std::shared_ptr<NormalModes>& normal_modes, int this_bead);

    void mpiCommunication() override;
    double getMomentumForCalc(const int ptcl_idx, const int axis) override;
    double& getMomentumForUpdate(const int ptcl_idx, const int axis) override;
    void updateCoupledMomenta() override;
private:
    std::shared_ptr<NormalModes> m_normal_modes;
    int m_this_bead;
};
