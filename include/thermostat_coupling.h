#pragma once

class Simulation;

// Classes to support coupling of the thermostat to Cartesian coords or normal modes of distinguishable ring polymers
class Coupling {
public:
    explicit Coupling(Simulation& _sim);
    virtual ~Coupling() = default;

    virtual void mpiCommunication() = 0;
    virtual double getMomentumForCalc(const int ptcl_idx, const int axis) = 0;
    virtual double& getMomentumForUpdate(const int ptcl_idx, const int axis) = 0;
    virtual void updateCoupledMomenta() = 0;
protected:
    Simulation& sim;
};

class CartesianCoupling : public Coupling {
public:
    CartesianCoupling(Simulation& _sim);
    virtual ~CartesianCoupling() = default;

    void mpiCommunication() override;
    double getMomentumForCalc(const int ptcl_idx, const int axis) override;
    double& getMomentumForUpdate(const int ptcl_idx, const int axis) override;
    void updateCoupledMomenta() override;    
};

class NMCoupling : public Coupling {
public:
    NMCoupling(Simulation& _sim);
    virtual ~NMCoupling() = default;

    void mpiCommunication() override;
    double getMomentumForCalc(const int ptcl_idx, const int axis) override;
    double& getMomentumForUpdate(const int ptcl_idx, const int axis) override;
    void updateCoupledMomenta() override;
};
