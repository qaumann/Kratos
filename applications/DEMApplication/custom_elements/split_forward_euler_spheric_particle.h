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

#if !defined(KRATOS_SPLIT_FORWARD_EULER_SPHERIC_PARTICLE_H_INCLUDED )
#define  KRATOS_SPLIT_FORWARD_EULER_SPHERIC_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "spheric_particle.h"
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

class KRATOS_API(DEM_APPLICATION) SplitForwardEulerSphericParticle : public SphericParticle
{

public:

/// Pointer definition of SplitForwardEulerSphericParticle
KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SplitForwardEulerSphericParticle);

/// Default constructor.
SplitForwardEulerSphericParticle();
SplitForwardEulerSphericParticle( IndexType NewId, GeometryType::Pointer pGeometry );
SplitForwardEulerSphericParticle( IndexType NewId, NodesArrayType const& ThisNodes);
SplitForwardEulerSphericParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

/// Destructor.
virtual ~SplitForwardEulerSphericParticle();

SplitForwardEulerSphericParticle& operator=(const SplitForwardEulerSphericParticle& rOther);

typedef typename SphericParticle::ParticleDataBuffer ParticleDataBuffer;

// typedef std::unique_ptr<ParticleDataBuffer> BufferPointerType;

void Initialize(const ProcessInfo& r_process_info) override;
void CalculateRightHandSide(const ProcessInfo& r_process_info, double dt, const array_1d<double,3>& gravity) override;

/// Turn back information as a string.
virtual std::string Info() const override
{
std::stringstream buffer;
buffer << "SplitForwardEulerSphericParticle" ;
return buffer.str();
}

/// Print information about this object.
virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "SplitForwardEulerSphericParticle";}

protected:

void EvaluateBallToBallForcesForPositiveIndentiations(SphericParticle::ParticleDataBuffer & data_buffer,
                                                        const ProcessInfo& r_process_info,
                                                        double LocalElasticContactForce[3],
                                                        double DeltDisp[3],
                                                        double LocalDeltDisp[3],
                                                        double RelVel[3],
                                                        const double indentation,
                                                        double ViscoDampingLocalContactForce[3],
                                                        double& cohesive_force,
                                                        SphericParticle* p_neighbour_element,
                                                        bool& sliding,
                                                        double LocalCoordSystem[3][3],
                                                        double OldLocalCoordSystem[3][3],
                                                        array_1d<double, 3>& neighbour_elastic_contact_force) override;

void ComputeBallToRigidFaceContactForce(ParticleDataBuffer & data_buffer,
                                            array_1d<double, 3>& rElasticForce,
                                            array_1d<double, 3>& rContactForce,
                                            double& RollingResistance,
                                            array_1d<double, 3>& rigid_element_force,
                                            const ProcessInfo& r_process_info) override;

virtual void ComputeBallToBallStiffness(SphericParticle::ParticleDataBuffer & data_buffer,
                                                                            double& r_nodal_stiffness,
                                                                            double& r_nodal_rotational_stiffness);

virtual void ComputeBallToRigidFaceStiffness(SphericParticle::ParticleDataBuffer & data_buffer,
                                                                                double& r_nodal_stiffness,
                                                                                double& r_nodal_rotational_stiffness);

virtual double CalculateStiffnessNorm(const double& r_normal_stiffness, const double& r_tangential_stiffness);

private:

friend class Serializer;

virtual void save(Serializer& rSerializer) const override
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericParticle );

    // public members

    // protected members
}

virtual void load(Serializer& rSerializer) override
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericParticle );

    // public members

    // protected members
}

}; // Class SplitForwardEulerSphericParticle

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
            SplitForwardEulerSphericParticle& rThis){ return rIStream;}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
            const SplitForwardEulerSphericParticle& rThis)
{
rThis.PrintInfo(rOStream);
rOStream << std::endl;
rThis.PrintData(rOStream);

return rOStream;
}

}  // namespace Kratos.

#endif // KRATOS_SPLIT_FORWARD_EULER_SPHERIC_PARTICLE_H_INCLUDED  defined
