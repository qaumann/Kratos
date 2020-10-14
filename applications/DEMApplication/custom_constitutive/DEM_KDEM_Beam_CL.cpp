// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "DEM_KDEM_Beam_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEM_Beam::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEM_Beam(*this));
        return p_clone;
    }

    void DEM_KDEM_Beam::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        KRATOS_INFO("DEM") << "Assigning DEM_KDEM_Beam to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        this->Check(pProp);
    }

    void DEM_KDEM_Beam::Check(Properties::Pointer pProp) const {
        DEMContinuumConstitutiveLaw::Check(pProp);

        if(!pProp->Has(ROTATIONAL_MOMENT_COEFFICIENT)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable ROTATIONAL_MOMENT_COEFFICIENT should be present in the properties when using DEM_KDEM. 1.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(ROTATIONAL_MOMENT_COEFFICIENT) = 1.0;
        }
    }

    void DEM_KDEM_Beam::CalculateElasticConstants(double& kn_el,
                                                  double& kt_el_0,
                                                  double& kt_el_1,
                                                  double initial_dist,
                                                  double equiv_young,
                                                  double equiv_poisson,
                                                  double calculation_area,
                                                  SphericContinuumParticle* element1,
                                                  SphericContinuumParticle* element2) {

        kn_el = equiv_young * calculation_area / initial_dist;

        const double Inertia_Ix = 0.5 * (element1->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_X] + element2->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_X]);
        const double Inertia_Iy = 0.5 * (element1->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_Y] + element2->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_Y]);

        // kt_el_0 = 3.0 * equiv_young * Inertia_Iy / (calculation_area * initial_dist);
        // kt_el_1 = 3.0 * equiv_young * Inertia_Ix / (calculation_area * initial_dist);

        // 2013 Bourrier. Discrete modeling of granular soils reinforcement by plant roots FOR CYLINDERS!!! CHECK!!!
        kt_el_0 = 12.0 * equiv_young * Inertia_Iy / (initial_dist * initial_dist * initial_dist);
        kt_el_1 = 12.0 * equiv_young * Inertia_Ix / (initial_dist * initial_dist * initial_dist);
    }

    void DEM_KDEM_Beam::CalculateViscoDampingCoeff(double& equiv_visco_damp_coeff_normal,
                                                   double& equiv_visco_damp_coeff_tangential_0,
                                                   double& equiv_visco_damp_coeff_tangential_1,
                                                   SphericContinuumParticle* element1,
                                                   SphericContinuumParticle* element2,
                                                   const double kn_el,
                                                   const double kt_el_0,
                                                   const double kt_el_1) {

        const double equiv_mass  = std::max(element1->GetMass(), element2->GetMass());
        const double aux_mass    = 0.5 * (element1->GetProperties()[BEAM_MASS] / equiv_mass);
        const double equiv_gamma = 0.5 * (element1->GetProperties()[DAMPING_GAMMA] + element2->GetProperties()[DAMPING_GAMMA]);

        equiv_visco_damp_coeff_normal       = 2.0 * equiv_gamma * aux_mass * sqrt(equiv_mass * kn_el  );
        equiv_visco_damp_coeff_tangential_0 = 2.0 * equiv_gamma * aux_mass * sqrt(equiv_mass * kt_el_0);
        equiv_visco_damp_coeff_tangential_1 = 2.0 * equiv_gamma * aux_mass * sqrt(equiv_mass * kt_el_1);
    }

    void DEM_KDEM_Beam::CalculateForces(const ProcessInfo& r_process_info,
                                        double OldLocalElasticContactForce[3],
                                        double LocalElasticContactForce[3],
                                        double LocalElasticExtraContactForce[3],
                                        double LocalCoordSystem[3][3],
                                        double LocalDeltDisp[3],
                                        const double kn_el,
                                        const double kt_el_0,
                                        const double kt_el_1,
                                        double& contact_sigma,
                                        double& contact_tau,
                                        double& failure_criterion_state,
                                        double equiv_young,
                                        double equiv_shear,
                                        double indentation,
                                        double calculation_area,
                                        double& acumulated_damage,
                                        SphericContinuumParticle* element1,
                                        SphericContinuumParticle* element2,
                                        int i_neighbour_count,
                                        int time_steps,
                                        bool& sliding,
                                        double &equiv_visco_damp_coeff_normal,
                                        double &equiv_visco_damp_coeff_tangential_0,
                                        double &equiv_visco_damp_coeff_tangential_1,
                                        double LocalRelVel[3],
                                        double ViscoDampingLocalContactForce[3]) {

        CalculateNormalForces(LocalElasticContactForce,
                              kn_el,
                              equiv_young,
                              indentation,
                              calculation_area,
                              acumulated_damage,
                              element1,
                              element2,
                              i_neighbour_count,
                              time_steps,
                              r_process_info);

        CalculateTangentialForces(OldLocalElasticContactForce,
                                  LocalElasticContactForce,
                                  LocalElasticExtraContactForce,
                                  LocalCoordSystem,
                                  LocalDeltDisp,
                                  kt_el_0,
                                  kt_el_1,
                                  equiv_young,
                                  contact_sigma,
                                  contact_tau,
                                  indentation,
                                  calculation_area,
                                  failure_criterion_state,
                                  element1,
                                  element2,
                                  i_neighbour_count,
                                  sliding,
                                  r_process_info);

        CalculateViscoDampingCoeff(equiv_visco_damp_coeff_normal,
                                   equiv_visco_damp_coeff_tangential_0,
                                   equiv_visco_damp_coeff_tangential_1,
                                   element1,
                                   element2,
                                   kn_el,
                                   kt_el_0,
                                   kt_el_1);

        CalculateViscoDamping(LocalRelVel,
                              ViscoDampingLocalContactForce,
                              indentation,
                              equiv_visco_damp_coeff_normal,
                              equiv_visco_damp_coeff_tangential_0,
                              equiv_visco_damp_coeff_tangential_1,
                              sliding,
                              element1->mIniNeighbourFailureId[i_neighbour_count]);
    }

    void DEM_KDEM_Beam::CalculateNormalForces(double LocalElasticContactForce[3],
                                              const double kn_el,
                                              double equiv_young,
                                              double indentation,
                                              double calculation_area,
                                              double& acumulated_damage,
                                              SphericContinuumParticle* element1,
                                              SphericContinuumParticle* element2,
                                              int i_neighbour_count,
                                              int time_steps,
                                              const ProcessInfo& r_process_info) {

            LocalElasticContactForce[2] = kn_el * indentation;
    }

    void DEM_KDEM_Beam::CalculateTangentialForces(double OldLocalElasticContactForce[3],
                                                  double LocalElasticContactForce[3],
                                                  double LocalElasticExtraContactForce[3],
                                                  double LocalCoordSystem[3][3],
                                                  double LocalDeltDisp[3],
                                                  const double kt_el_0,
                                                  const double kt_el_1,
                                                  const double equiv_young,
                                                  double& contact_sigma,
                                                  double& contact_tau,
                                                  double indentation,
                                                  double calculation_area,
                                                  double& failure_criterion_state,
                                                  SphericContinuumParticle* element1,
                                                  SphericContinuumParticle* element2,
                                                  int i_neighbour_count,
                                                  bool& sliding,
                                                  const ProcessInfo& r_process_info) {

        LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - kt_el_0 * LocalDeltDisp[0]; // 0: first tangential
        LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - kt_el_1 * LocalDeltDisp[1]; // 1: second tangential
    }

    void DEM_KDEM_Beam::CalculateViscoDamping(double LocalRelVel[3],
                                              double ViscoDampingLocalContactForce[3],
                                              double indentation,
                                              double equiv_visco_damp_coeff_normal,
                                              double equiv_visco_damp_coeff_tangential_0,
                                              double equiv_visco_damp_coeff_tangential_1,
                                              bool& sliding,
                                              int failure_id) {

        ViscoDampingLocalContactForce[2] = -equiv_visco_damp_coeff_normal       * LocalRelVel[2];
        ViscoDampingLocalContactForce[0] = -equiv_visco_damp_coeff_tangential_0 * LocalRelVel[0];
        ViscoDampingLocalContactForce[1] = -equiv_visco_damp_coeff_tangential_1 * LocalRelVel[1];
    }

    void DEM_KDEM_Beam::ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                         SphericContinuumParticle* neighbor,
                                                         double equiv_young,
                                                         double distance,
                                                         double calculation_area,
                                                         double LocalCoordSystem[3][3],
                                                         double ElasticLocalRotationalMoment[3],
                                                         double ViscoLocalRotationalMoment[3],
                                                         double equiv_poisson,
                                                         double indentation) {

        KRATOS_TRY
        double LocalDeltaRotatedAngle[3]    = {0.0};
        double LocalDeltaAngularVelocity[3] = {0.0};

        array_1d<double, 3> GlobalDeltaRotatedAngle;
        noalias(GlobalDeltaRotatedAngle) = element->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE) - neighbor->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);
        array_1d<double, 3> GlobalDeltaAngularVelocity;
        noalias(GlobalDeltaAngularVelocity) = element->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY) - neighbor->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);

        GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalDeltaRotatedAngle, LocalDeltaRotatedAngle);
        GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalDeltaAngularVelocity, LocalDeltaAngularVelocity);

        const double MomentOfInertiaX = std::max(element->GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[1], neighbor->GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[1]);
        const double MomentOfInertiaY = std::max(element->GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[2], neighbor->GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[2]);
        const double MomentOfInertiaZ = std::max(element->GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[0], neighbor->GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[0]);

        double norm_distance = (element->GetRadius() + neighbor->GetRadius()) / distance; // If spheres are not tangent the Damping coefficient, DeltaRotatedAngle and DeltaAngularVelocity have to be normalized

        const double auxX = 0.5 * element->GetProperties()[BEAM_PRINCIPAL_MOMENTS_OF_INERTIA_Y] * element->GetProperties()[BEAM_MASS] * norm_distance / MomentOfInertiaX;
        const double auxY = 0.5 * element->GetProperties()[BEAM_PRINCIPAL_MOMENTS_OF_INERTIA_Z] * element->GetProperties()[BEAM_MASS] * norm_distance / MomentOfInertiaY;
        const double auxZ = 0.5 * element->GetProperties()[BEAM_PRINCIPAL_MOMENTS_OF_INERTIA_X] * element->GetProperties()[BEAM_MASS] * norm_distance / MomentOfInertiaZ;

        const double equiv_shear   = equiv_young / (2.0 * (1 + equiv_poisson));

        const double Inertia_Ix = 0.5 * (element->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_X] + neighbor->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_X]);
        const double Inertia_Iy = 0.5 * (element->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_Y] + neighbor->GetProperties()[BEAM_PLANAR_MOMENT_OF_INERTIA_Y]);
        const double Inertia_J  = Inertia_Ix + Inertia_Iy;

        const double my_gamma    = element->GetProperties()[DAMPING_GAMMA];
        const double other_gamma = neighbor->GetProperties()[DAMPING_GAMMA];
        const double equiv_gamma = 0.5 * (my_gamma + other_gamma);

        const double k_rot_x = equiv_young * Inertia_Ix * norm_distance / distance;
        const double k_rot_y = equiv_young * Inertia_Iy * norm_distance / distance;
        const double k_tor   = equiv_shear * Inertia_J  / distance;

        const double visc_param_rot_x = 2.0 * equiv_gamma * auxX * sqrt(MomentOfInertiaX * k_rot_x);
        const double visc_param_rot_y = 2.0 * equiv_gamma * auxY * sqrt(MomentOfInertiaY * k_rot_y);
        const double visc_param_tor   = 2.0 * equiv_gamma * auxZ * sqrt(MomentOfInertiaZ * k_tor  );

        ElasticLocalRotationalMoment[0] = -k_rot_x * LocalDeltaRotatedAngle[0];
        ElasticLocalRotationalMoment[1] = -k_rot_y * LocalDeltaRotatedAngle[1];
        ElasticLocalRotationalMoment[2] = -k_tor   * LocalDeltaRotatedAngle[2];

        ViscoLocalRotationalMoment[0] = -visc_param_rot_x * LocalDeltaAngularVelocity[0];
        ViscoLocalRotationalMoment[1] = -visc_param_rot_y * LocalDeltaAngularVelocity[1];
        ViscoLocalRotationalMoment[2] = -visc_param_tor   * LocalDeltaAngularVelocity[2];

        KRATOS_CATCH("")
    }//ComputeParticleRotationalMoments

} // namespace Kratos
