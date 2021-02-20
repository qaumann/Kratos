// KRATOS
//
//  License:		 BSD License
//					 license: ../../license.txt
//
//  Main authors:
//
//

// System includes

// External includes


// Project includes
#include "utilities/geometry_utilities.h"
// Application includes
#include "custom_elements/acoustic_element.h"
#include "custom_elements/acoustic_pml_element.h"


namespace Kratos
{
AcousticPMLElement::AcousticPMLElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : AcousticElement( NewId, pGeometry )
{

}

/***********************************************************************************/
/***********************************************************************************/

AcousticPMLElement::AcousticPMLElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : AcousticElement( NewId, pGeometry, pProperties )
{

}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer AcousticPMLElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<AcousticPMLElement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer AcousticPMLElement::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<AcousticPMLElement>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

AcousticPMLElement::~AcousticPMLElement()
{

}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer AcousticPMLElement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    AcousticPMLElement::Pointer p_new_elem = Kratos::make_intrusive<AcousticPMLElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

int  AcousticPMLElement::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    int ier = Element::Check(rCurrentProcessInfo);

    return ier;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticPMLElement::Initialize()
{
    KRATOS_TRY
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticPMLElement::EquationIdVector(EquationIdVectorType& rResult,
                                        ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const SizeType num_nodes = GetGeometry().PointsNumber();

    if(rResult.size() != num_nodes)
        rResult.resize(num_nodes,false);

    for (SizeType i_node = 0; i_node < num_nodes; i_node++)
        rResult[i_node] = GetGeometry()[i_node].GetDof(PRESSURE).EquationId();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticPMLElement::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    const SizeType num_nodes = GetGeometry().PointsNumber();


    if(rElementalDofList.size() != num_nodes)
        rElementalDofList.resize(num_nodes);

    for (SizeType i_node = 0; i_node < num_nodes; i_node++)
        rElementalDofList[i_node] = GetGeometry()[i_node].pGetDof(PRESSURE);
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticPMLElement::ComplexJacobian( ComplexMatrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod, ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& geom = GetGeometry();
    const SizeType working_space_dimension = geom.WorkingSpaceDimension();
    const SizeType local_space_dimension = geom.LocalSpaceDimension();

    if(rResult.size1() != working_space_dimension || rResult.size2() != local_space_dimension)
        rResult.resize( working_space_dimension, local_space_dimension, false );


    const Matrix& r_shape_functions_gradient_in_integration_point = geom.ShapeFunctionsLocalGradients( ThisMethod )[ IntegrationPointIndex ];

    rResult.clear();

    const SizeType points_number = geom.PointsNumber();


    for (IndexType j = 0; j < points_number; ++j ) {
        const array_1d<double, 3>& r_coordinates = geom[j].Coordinates();

        // PML_FACTOR is calculated in python
        const double factor = geom[j].GetValue(PML_FACTOR);
        const double bulk = GetProperties()[BULK_MODULUS];
        const double density = GetProperties()[DENSITY];
        const double wavespeed = std::sqrt(bulk/density);

        // Absorption function
        const double wavenum = rCurrentProcessInfo[FREQUENCY]/wavespeed;
        const array_1d<double, 3> r_imag = geom[j].GetValue(PML_DIRECTION) * factor / wavenum;

        for(IndexType k = 0; k< working_space_dimension; ++k)
        {
            const std::complex<double> value(r_coordinates[k], -r_imag[k]);  // Add imaginary part to coordinates (complex coordinate stretching)
            for(IndexType m = 0; m < local_space_dimension; ++m) {
                rResult(k,m) += value * r_shape_functions_gradient_in_integration_point(j,m);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticPMLElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    SizeType number_of_nodes = GetGeometry().PointsNumber();
    VectorType rRightHandSideVector = ZeroVector( number_of_nodes );
    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, true, false);
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticPMLElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& geom = GetGeometry();
    IntegrationMethod ThisIntegrationMethod = geom.GetDefaultIntegrationMethod();
    const SizeType number_of_nodes = geom.PointsNumber();


    const GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints(ThisIntegrationMethod);
    IndexType NumGauss = integration_points.size();
    const Matrix& NContainer = geom.ShapeFunctionsValues(ThisIntegrationMethod);

    SizeType working_space_dimension = geom.WorkingSpaceDimension();
    SizeType local_space_dimension = geom.LocalSpaceDimension();

    if( rMassMatrix.size1() != number_of_nodes || rMassMatrix.size2() != number_of_nodes ) {
        rMassMatrix.resize(number_of_nodes, number_of_nodes, false);
    }

    noalias(rMassMatrix) = ZeroMatrix( number_of_nodes, number_of_nodes );
    ComplexMatrix ComplexMassMatrix = ComplexZeroMatrix( number_of_nodes, number_of_nodes );

    const double bulk = GetProperties()[BULK_MODULUS];

    for (IndexType g = 0; g < NumGauss; g++)
    {
        ComplexMatrix Jac (working_space_dimension, local_space_dimension);
        ComplexMatrix Jinv (working_space_dimension, local_space_dimension);
        ComplexJacobian(Jac, g, ThisIntegrationMethod, rCurrentProcessInfo);
        std::complex<double> DetJ;

        MathUtils<std::complex<double>>::InvertMatrix(Jac, Jinv, DetJ);
        for (IndexType i = 0; i < number_of_nodes; i++)
        {
            for (IndexType j = 0; j < number_of_nodes; j++)
            {
                double GaussWeight = integration_points[g].Weight();
                ComplexMassMatrix(i,j) += DetJ * NContainer(g, i) * NContainer(g, j) * GaussWeight / bulk;
            }
        }
    }

    if (rCurrentProcessInfo[BUILD_LEVEL] == 291) {
        for (IndexType i =0; i < number_of_nodes; i++) {
            for (IndexType j = 0; j < number_of_nodes; j++) {
                rMassMatrix(i, j) = ComplexMassMatrix(i,j).imag();
            }
        }
    } else {
        for (IndexType i =0; i < number_of_nodes; i++) {
            for (IndexType j = 0; j < number_of_nodes; j++) {
                rMassMatrix(i, j) = ComplexMassMatrix(i,j).real();
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticPMLElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{

    const GeometryType& geom = GetGeometry();
    const SizeType number_of_nodes = geom.PointsNumber();
     MatrixType rLeftHandSideVector = ZeroMatrix(number_of_nodes, number_of_nodes);

    // Resizing as needed the RHS

    if ( rRightHandSideVector.size() != number_of_nodes )
        rRightHandSideVector.resize( number_of_nodes, false );

    rRightHandSideVector = ZeroVector( number_of_nodes ); //resetting RHS


    CalculateAll(rLeftHandSideVector, rRightHandSideVector, rCurrentProcessInfo, true, true);
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticPMLElement::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag)
{
    const GeometryType& geom = GetGeometry();
    IntegrationMethod ThisIntegrationMethod = geom.GetDefaultIntegrationMethod();
    const SizeType number_of_nodes = geom.PointsNumber();
    const GeometryType::IntegrationPointsArrayType integration_points = geom.IntegrationPoints(ThisIntegrationMethod);
    const SizeType working_space_dimension = geom.WorkingSpaceDimension();
    const SizeType local_space_dimension = geom.LocalSpaceDimension();

    if (CalculateStiffnessMatrixFlag) {

        if( rLeftHandSideMatrix.size1() != number_of_nodes || rLeftHandSideMatrix.size2() != number_of_nodes ) {
            rLeftHandSideMatrix.resize(number_of_nodes, number_of_nodes, false);
        }

        noalias(rLeftHandSideMatrix) = ZeroMatrix( number_of_nodes, number_of_nodes );

        // for solving acoustic problem
        if (rCurrentProcessInfo[BUILD_LEVEL] < 100 ) {
            const double density = GetProperties()[DENSITY];
            ComplexMatrix ComplexLeftHandSideMatrix = ComplexZeroMatrix( number_of_nodes, number_of_nodes );
            ShapeFunctionDerivativesArrayType DN_DE = geom.ShapeFunctionsLocalGradients(ThisIntegrationMethod);
            std::complex<double> DetJ;

            for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number )
            {
                ComplexMatrix Jac (working_space_dimension, local_space_dimension);
                ComplexMatrix Jinv (working_space_dimension, local_space_dimension);
                ComplexJacobian(Jac, point_number, ThisIntegrationMethod, rCurrentProcessInfo);
                MathUtils<std::complex<double>>::InvertMatrix(Jac, Jinv, DetJ);
                double weight = integration_points[point_number].Weight();
                ComplexMatrix DN_DX = prod( DN_DE[point_number], Jinv);
                noalias( ComplexLeftHandSideMatrix ) += DetJ * weight/density * prod( DN_DX, trans(DN_DX));
            }

            if (rCurrentProcessInfo[BUILD_LEVEL] == 91) {
                for (IndexType i =0; i < number_of_nodes; i++) {
                    for (IndexType j = 0; j < number_of_nodes; j++) {
                        rLeftHandSideMatrix(i, j) = ComplexLeftHandSideMatrix(i,j).imag();
                    }
                }
            } else {
                for (IndexType i =0; i < number_of_nodes; i++) {
                    for (IndexType j = 0; j < number_of_nodes; j++) {
                        rLeftHandSideMatrix(i, j) = ComplexLeftHandSideMatrix(i,j).real();
                    }
                }
            }
        }

        // for solving problem of pml_direction_process
        if (rCurrentProcessInfo[BUILD_LEVEL] == 401) {
            noalias(rLeftHandSideMatrix) = ZeroMatrix( number_of_nodes, number_of_nodes );
            Vector DetJ;
            ShapeFunctionDerivativesArrayType DN_DX;
            geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, DetJ, ThisIntegrationMethod);
            for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
                double int_weight = integration_points[point_number].Weight() * DetJ(point_number);
                rLeftHandSideMatrix += int_weight * prod( DN_DX[point_number], trans(DN_DX[point_number]));
            }
        }
    }

    // RHS
    if(CalculateResidualVectorFlag)
    {
        // For acoustic problem
        if ( rRightHandSideVector.size() != number_of_nodes ) {
            rRightHandSideVector.resize( number_of_nodes, false );
        }
        rRightHandSideVector = ZeroVector( number_of_nodes ); //resetting RHS

        // Dirichlet condition for pml_direction_process
        if (rCurrentProcessInfo[BUILD_LEVEL] == 401) {
            VectorType p_potential = ZeroVector(number_of_nodes);
            for (SizeType i = 0; i < number_of_nodes; i++) {
                p_potential[i] = geom[i].GetValue(PML_PRESCRIBED_POTENTIAL);
            }
            noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, p_potential);
        }
    }


}

/***********************************************************************************/
/***********************************************************************************/

void AcousticPMLElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, true, true);
    KRATOS_CATCH( "" )

}

} // Namespace Kratos