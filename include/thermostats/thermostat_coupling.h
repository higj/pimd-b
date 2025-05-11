#include "common.h"
#pragma once

class NormalModes;

// Classes to support coupling of the thermostat to Cartesian coords or normal modes of distinguishable ring polymers
class Coupling {
public:
    explicit Coupling(dVec& _momenta);
    virtual ~Coupling() = default;

    virtual void mpiCommunication() = 0;
    virtual double getMomentumForCalc(const int ptcl_idx, const int axis) = 0;
    virtual double& getMomentumForUpdate(const int ptcl_idx, const int axis) = 0;
    virtual void updateCoupledMomenta() = 0;
protected:
    dVec& momenta;
};

class CartesianCoupling : public Coupling {
public:
    CartesianCoupling(dVec& _momenta);
    virtual ~CartesianCoupling() = default;

    void mpiCommunication() override;
    double getMomentumForCalc(const int ptcl_idx, const int axis) override;
    double& getMomentumForUpdate(const int ptcl_idx, const int axis) override;
    void updateCoupledMomenta() override;    
};

class NMCoupling : public Coupling {
public:
    NMCoupling(dVec& _momenta, NormalModes& normal_modes, int this_bead);
    virtual ~NMCoupling() = default;

    void mpiCommunication() override;
    double getMomentumForCalc(const int ptcl_idx, const int axis) override;
    double& getMomentumForUpdate(const int ptcl_idx, const int axis) override;
    void updateCoupledMomenta() override;
private:
    NormalModes& normal_modes;
    int this_bead;
};
