// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

// System includes

// External includes


// Project includes
#include "utilities/math_utils.h"

// Application includes
#include "custom_elements/small_displacement_explicit_split_scheme.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
SmallDisplacementExplicitSplitScheme::SmallDisplacementExplicitSplitScheme( IndexType NewId, GeometryType::Pointer pGeometry )
    : SmallDisplacement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

SmallDisplacementExplicitSplitScheme::SmallDisplacementExplicitSplitScheme( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : SmallDisplacement( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementExplicitSplitScheme::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<SmallDisplacementExplicitSplitScheme>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementExplicitSplitScheme::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<SmallDisplacementExplicitSplitScheme>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

SmallDisplacementExplicitSplitScheme::~SmallDisplacementExplicitSplitScheme()
{
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementExplicitSplitScheme::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    SmallDisplacementExplicitSplitScheme::Pointer p_new_elem = Kratos::make_intrusive<SmallDisplacementExplicitSplitScheme>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(BaseType::mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(BaseType::mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementExplicitSplitScheme::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<double>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // TODO: seguir

    auto& r_geom = this->GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType mat_size = number_of_nodes * dimension;

    // Compiting the nodal mass
    if (rDestinationVariable == NODAL_MASS ) {
        VectorType element_mass_vector(mat_size);
        this->CalculateLumpedMassVector(element_mass_vector);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * dimension;

            #pragma omp atomic
            r_geom[i].GetValue(NODAL_MASS) += element_mass_vector[index];
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementExplicitSplitScheme::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // TODO: seguir

    // Operation performed: rRightHandSideVector -= IntForce * IntegrationWeight
    noalias( rRightHandSideVector ) -= IntegrationWeight * prod( trans( rThisKinematicVariables.B ), rStressVector );

        auto& r_geometry = this->GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType dimension = r_geometry.WorkingSpaceDimension();
        const SizeType strain_size = GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Resizing as needed the LHS
        const SizeType mat_size = number_of_nodes * dimension;

        if ( CalculateStiffnessMatrixFlag ) { // Calculation of the matrix is required
            if ( rLeftHandSideMatrix.size1() != mat_size )
                rLeftHandSideMatrix.resize( mat_size, mat_size, false );

            noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
        }

        // Resizing as needed the RHS
        if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
            if ( rRightHandSideVector.size() != mat_size )
                rRightHandSideVector.resize( mat_size, false );

            rRightHandSideVector = ZeroVector( mat_size ); //resetting RHS
        }

        // Reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

        ConstitutiveLaw::Parameters Values(r_geometry,GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        if ( CalculateStiffnessMatrixFlag ) {
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        } else {
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        }

        // If strain has to be computed inside of the constitutive law with PK2
        Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

        // Some declarations
        array_1d<double, 3> body_force;
        double int_to_reference_weight;

        // Computing in all integrations points
        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            // Contribution to external forces
            noalias(body_force) = this->GetBodyForce(integration_points, point_number);

            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

            // Compute material reponse
            CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure());

            // Calculating weights for integration on the reference configuration
            int_to_reference_weight = GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0);

            if ( dimension == 2 && GetProperties().Has( THICKNESS ))
                int_to_reference_weight *= GetProperties()[THICKNESS];

            if ( CalculateStiffnessMatrixFlag ) { // Calculation of the matrix is required
                // Contributions to stiffness matrix calculated on the reference config
                this->CalculateAndAddKm( rLeftHandSideMatrix, this_kinematic_variables.B, this_constitutive_variables.D, int_to_reference_weight );
            }

            if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
                this->CalculateAndAddResidualVector(rRightHandSideVector, this_kinematic_variables, rCurrentProcessInfo, body_force, this_constitutive_variables.StressVector, int_to_reference_weight);
            }
        }



    auto& r_geom = this->GetGeometry();
    const auto& r_prop = this->GetProperties();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType element_size = dimension * number_of_nodes;

    Vector damping_residual_contribution = ZeroVector(element_size);

    // Calculate damping contribution to residual -->
    if (r_prop.Has(RAYLEIGH_ALPHA) || r_prop.Has(RAYLEIGH_BETA)) {
        Vector current_nodal_velocities = ZeroVector(element_size);
        this->GetFirstDerivativesVector(current_nodal_velocities);

        Matrix damping_matrix(element_size, element_size);
        this->CalculateDampingMatrixWithLumpedMass(damping_matrix, rCurrentProcessInfo);

        // Current residual contribution due to damping
        noalias(damping_residual_contribution) = prod(damping_matrix, current_nodal_velocities);
    }

    // Computing the force residual
    if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = dimension * i;

            array_1d<double, 3>& r_force_residual = r_geom[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

            for (IndexType j = 0; j < dimension; ++j) {
                #pragma omp atomic
                r_force_residual[j] += rRHSVector[index + j] - damping_residual_contribution[index + j];
            }
        }
    }

    KRATOS_CATCH("")
}


/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementExplicitSplitScheme::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SmallDisplacement );
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementExplicitSplitScheme::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SmallDisplacement );
}

} // Namespace Kratos


