//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//
//

#if !defined(KRATOS_SPLIT_FORWARD_EULER_CYLINDER_CONTINUUM_PARTICLE_H_INCLUDED )
#define  KRATOS_SPLIT_FORWARD_EULER_CYLINDER_CONTINUUM_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "split_forward_euler_spheric_continuum_particle.h"
#include "../custom_utilities/AuxiliaryFunctions.h"
#include "../custom_constitutive/DEM_discontinuum_constitutive_law.h"
#include "../custom_conditions/RigidFace.h"
#include "../custom_conditions/dem_wall.h"
#include "../custom_strategies/schemes/dem_integration_scheme.h"
#include "includes/kratos_export_api.h"
#include "../custom_utilities/properties_proxies.h"
#include "includes/kratos_flags.h"

namespace Kratos
{

class DEMWall;

class KRATOS_API(DEM_APPLICATION) SplitForwardEulerCylinderContinuumParticle : public SplitForwardEulerSphericContinuumParticle
{

public:

/// Pointer definition of SplitForwardEulerCylinderContinuumParticle
KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SplitForwardEulerCylinderContinuumParticle);

/// Default constructor.
SplitForwardEulerCylinderContinuumParticle( IndexType NewId, GeometryType::Pointer pGeometry );
SplitForwardEulerCylinderContinuumParticle( IndexType NewId, NodesArrayType const& ThisNodes);
SplitForwardEulerCylinderContinuumParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

/// Destructor.
virtual ~SplitForwardEulerCylinderContinuumParticle();


/// Turn back information as a string.
virtual std::string Info() const override
{
std::stringstream buffer;
buffer << "SplitForwardEulerCylinderContinuumParticle" ;
return buffer.str();
}

/// Print information about this object.
virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "SplitForwardEulerCylinderContinuumParticle";}

void ContactAreaWeighting() override;

void FinalizeStressTensor(const ProcessInfo& r_process_info, double& rRepresentative_Volume) override;

protected:

SplitForwardEulerCylinderContinuumParticle();

double CalculateVolume() override;
double CalculateMomentOfInertia() override;
void AddContributionToRepresentativeVolume(const double distance, const double radius_sum, const double contact_area) override ;

double CalculateStiffnessNorm(const double& r_normal_stiffness, const double& r_tangential_stiffness) override;

private:

friend class Serializer;

virtual void save(Serializer& rSerializer) const override
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SplitForwardEulerSphericContinuumParticle );

    // public members

    // protected members
}

virtual void load(Serializer& rSerializer) override
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SplitForwardEulerSphericContinuumParticle );

    // public members

    // protected members
}

// SplitForwardEulerCylinderContinuumParticle& operator=(const SplitForwardEulerCylinderContinuumParticle& rOther);


}; // Class SplitForwardEulerCylinderContinuumParticle

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
            SplitForwardEulerCylinderContinuumParticle& rThis){ return rIStream;}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
            const SplitForwardEulerCylinderContinuumParticle& rThis)
{
rThis.PrintInfo(rOStream);
rOStream << std::endl;
rThis.PrintData(rOStream);

return rOStream;
}

}  // namespace Kratos.

#endif // KRATOS_SPLIT_FORWARD_EULER_CYLINDER_CONTINUUM_PARTICLE_H_INCLUDED  defined
