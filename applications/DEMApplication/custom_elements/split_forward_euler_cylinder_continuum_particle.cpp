// System includes
#include <string>
#include <iostream>
#include <cmath>

#include <fstream>

// External includes

// Project includes
#include "split_forward_euler_cylinder_continuum_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "custom_utilities/discrete_particle_configure.h"
#include "custom_strategies/schemes/glued_to_wall_scheme.h"


namespace Kratos
{
// using namespace GeometryFunctions;

SplitForwardEulerCylinderContinuumParticle::SplitForwardEulerCylinderContinuumParticle()
    : SplitForwardEulerSphericContinuumParticle() {
}

SplitForwardEulerCylinderContinuumParticle::SplitForwardEulerCylinderContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry)
    : SplitForwardEulerSphericContinuumParticle(NewId, pGeometry) {
}

SplitForwardEulerCylinderContinuumParticle::SplitForwardEulerCylinderContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : SplitForwardEulerSphericContinuumParticle(NewId, pGeometry, pProperties) {
}

SplitForwardEulerCylinderContinuumParticle::SplitForwardEulerCylinderContinuumParticle(IndexType NewId, NodesArrayType const& ThisNodes)
    : SplitForwardEulerSphericContinuumParticle(NewId, ThisNodes) {
}

Element::Pointer SplitForwardEulerCylinderContinuumParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SplitForwardEulerCylinderContinuumParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

/// Destructor.
SplitForwardEulerCylinderContinuumParticle::~SplitForwardEulerCylinderContinuumParticle(){
}

// SplitForwardEulerCylinderContinuumParticle& SplitForwardEulerCylinderContinuumParticle::operator=(const SplitForwardEulerCylinderContinuumParticle& rOther) {
//     SplitForwardEulerSphericContinuumParticle::operator=(rOther);

//     return *this;
// }

void SplitForwardEulerCylinderContinuumParticle::ContactAreaWeighting() //MISMI 10: POOYAN this could be done by calculating on the bars. not looking at the neighbous of my neighbours.
{

    double alpha = 1.0;
    double circle_perimeter = 2*Globals::Pi * GetRadius();
    double total_equiv_perimeter = 0.0;
    unsigned int continuous_initial_neighbours_size = mContinuumInitialNeighborsSize;
    Vector& cont_ini_neigh_area = GetValue(NEIGHBOURS_CONTACT_AREAS);

    for (unsigned int i = 0; i < continuous_initial_neighbours_size; i++) {
        SphericParticle* ini_cont_neighbour_iterator = mNeighbourElements[i];
        double other_radius     = ini_cont_neighbour_iterator->GetInteractionRadius();
        double area = mContinuumConstitutiveLawArray[i]->CalculateContactArea(GetRadius(), other_radius, cont_ini_neigh_area); //This call fills the vector of areas only if the Constitutive Law wants.
        total_equiv_perimeter += area;
    } //for every neighbour

    if (continuous_initial_neighbours_size >= 4) {

        if (!IsSkin()) {
            AuxiliaryFunctions::CalculateAlphaFactor2D(continuous_initial_neighbours_size, circle_perimeter, total_equiv_perimeter, alpha);
            for (unsigned int i = 0; i < cont_ini_neigh_area.size(); i++) {
                cont_ini_neigh_area[i] = alpha*cont_ini_neigh_area[i];
            } //for every neighbour
        }
        else { //skin sphere
            for (unsigned int i = 0; i < cont_ini_neigh_area.size(); i++) {
                alpha            = 1.30*(1.10266)*(circle_perimeter/total_equiv_perimeter)*((double(continuous_initial_neighbours_size))/6); // 6 is mean coordination number.
                cont_ini_neigh_area[i] = alpha*cont_ini_neigh_area[i];
            }//loop on cont neighs
        }//skin particles.
    }//if 3 neighbours or more.
}//Contact Area Weighting

double SplitForwardEulerCylinderContinuumParticle::CalculateVolume() {
    return Globals::Pi * GetRadius() * GetRadius();
}

double SplitForwardEulerCylinderContinuumParticle::CalculateMomentOfInertia() {
    return 0.5 * GetMass() * GetRadius() * GetRadius();
}

void SplitForwardEulerCylinderContinuumParticle::FinalizeStressTensor(const ProcessInfo& r_process_info, double& rRepresentative_Volume){

    KRATOS_TRY
    SphericParticle::FinalizeStressTensor(r_process_info, rRepresentative_Volume);

    if (!r_process_info[IMPOSED_Z_STRAIN_OPTION]) return;

    double z_strain_value = r_process_info[IMPOSED_Z_STRAIN_VALUE];
    double myYoung = GetYoung();
    double myPoisson = GetPoisson();

    // (*mStressTensor)(2,2) += E*z_displacement - poisson*(sigma_xx + sigma_yy);
    (*mStressTensor)(2, 2) = myYoung*z_strain_value + myPoisson*((*mStressTensor)(0, 0) + (*mStressTensor)(1, 1));

    KRATOS_CATCH("")
}

void SplitForwardEulerCylinderContinuumParticle::AddContributionToRepresentativeVolume(const double distance,
                                                                    const double radius_sum,
                                                                    const double contact_area) {

    KRATOS_TRY

    double gap = distance - radius_sum;
    double real_distance = GetInteractionRadius() + 0.5 * gap;
    double& rRepresentative_Volume = this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME);
    rRepresentative_Volume += 0.5 * (real_distance * contact_area);

    KRATOS_CATCH("")
}

double SplitForwardEulerCylinderContinuumParticle::CalculateStiffnessNorm(const double& r_normal_stiffness, const double& r_tangential_stiffness) {
    return std::sqrt(r_normal_stiffness*r_normal_stiffness+r_tangential_stiffness*r_tangential_stiffness);
}

}  // namespace Kratos.
