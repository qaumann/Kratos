// KRATOS
//
//  License:		 BSD License
//					 license: ../../license.txt
//
//  Main authors:	Sebastian Schopper
//
//

// System includes

// External includes

// Project includes
#include "utilities/geometry_utilities.h"
// Application includes
#include "custom_elements/porous_element.h"

namespace Kratos
{
PorousElement::PorousElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

PorousElement::PorousElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : Element( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer PorousElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<PorousElement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer PorousElement::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<PorousElement>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

PorousElement::~PorousElement()
{
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer PorousElement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    PorousElement::Pointer p_new_elem = Kratos::make_intrusive<PorousElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

int PorousElement::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY
    const double numerical_limit = std::numeric_limits<double>::epsilon();

    int ier = Element::Check(rCurrentProcessInfo);

	// Check that all required variables have resonable values
    KRATOS_ERROR_IF(!GetProperties().Has(DENSITY_SOLID) ||
                GetProperties()[DENSITY_SOLID] <= numerical_limit)
        << "Please provide a reasonable value for \"DENSITY_SOLID\" for element #"
        << Id() << std::endl;

	//TODO: add sanity check

    return ier;

    KRATOS_CATCH( "" );
}

void PorousElement::Initialize()
{
    KRATOS_TRY
    // std::cout << "I am initializing a porous element\n";
    KRATOS_CATCH("")
}


/***********************************************************************************/
/***********************************************************************************/

void PorousElement::EquationIdVector(EquationIdVectorType& rResult,
                                        ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

	const SizeType number_of_nodes = GetGeometry().size();
	const SizeType dimension = GetGeometry().WorkingSpaceDimension();

	// 2D -> 3 dofs, 3D -> 4 dofs per node (pressure + displacement)
	const SizeType dofs_size = number_of_nodes * (dimension + 1);

    if(rResult.size() != dofs_size) {
        rResult.resize(dofs_size,false);
	}

	SizeType index;
	// order of the dofs: displacement, pressure
	if (dimension == 2) {
		for (SizeType i = 0; i < number_of_nodes; i++) {
			index = i * 2;
			rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
			rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
		}
		for (SizeType i = 0; i < number_of_nodes; i++)
		{
			index = 2 * number_of_nodes + i;
			rResult[index] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
		}
	}
	else {
		for (SizeType i = 0; i < number_of_nodes; i++) {
			index = i * 3;
			rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
			rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
			rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
		}
		for (SizeType i = 0; i < number_of_nodes; i++)
		{
			index = 3 * number_of_nodes + i;
			rResult[index] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
		}
	}

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void PorousElement::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo)
{
	const SizeType dimension = GetGeometry().WorkingSpaceDimension();

	rElementalDofList.resize(0);

	// order: displacements of all nodes first, then pressure dofs for all nodes at the end

	if (dimension == 2) {
		for (SizeType i = 0; i < GetGeometry().size(); i++)
		{

			rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
			rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
		}
		for (SizeType i = 0; i < GetGeometry().size(); i++)
		{
			rElementalDofList.push_back(GetGeometry()[i].pGetDof(PRESSURE));
		}
	}
	else {
		for (SizeType i = 0; i < GetGeometry().size(); i++)
		{

			rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
			rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
			rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
		}
		for (SizeType i = 0; i < GetGeometry().size(); i++)
		{
			rElementalDofList.push_back(GetGeometry()[i].pGetDof(PRESSURE));

		}
	}
}

/***********************************************************************************/
/***********************************************************************************/

// double AcousticElement::CalculateDerivativesOnReferenceConfiguration(
//     Matrix& rJ0,
//     Matrix& rInvJ0,
//     Matrix& rDN_DX,
//     const IndexType PointNumber,
//     IntegrationMethod ThisIntegrationMethod
//     ) const
// {
//     const GeometryType& r_geom = GetGeometry();
//     GeometryUtils::JacobianOnInitialConfiguration(
//         r_geom,
//         r_geom.IntegrationPoints(ThisIntegrationMethod)[PointNumber], rJ0);
//     double detJ0;
//     MathUtils<double>::InvertMatrix(rJ0, rInvJ0, detJ0);
//     const Matrix& rDN_De =
//         GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[PointNumber];
//     GeometryUtils::ShapeFunctionsGradients(rDN_De, rInvJ0, rDN_DX);
//     return detJ0;
// }

/***********************************************************************************/
/***********************************************************************************/

Matrix PorousElement::CalculateNMatrix(const Vector& rN)
{
	const auto& r_geometry = GetGeometry();
	const SizeType number_of_nodes = r_geometry.PointsNumber();
	const SizeType dimension = r_geometry.WorkingSpaceDimension();

	// initialize result
	Matrix N_mat;
	N_mat = ZeroMatrix(dimension, number_of_nodes * dimension);

	// loop over Matrix to fill the entries
	for (IndexType node = 0; node < number_of_nodes; node++)
	{
		for (IndexType dim = 0; dim < dimension; dim++)
		{
			N_mat(dim, node * dimension + dim) = rN(node);
		}
	}

	return N_mat;
}

/***********************************************************************************/
/***********************************************************************************/

Vector PorousElement::CalculateBuMatrix(
	const Matrix& rDN_DX)
{
	const auto& r_geometry = GetGeometry();
	const SizeType number_of_nodes = r_geometry.PointsNumber();
	const SizeType dimension = r_geometry.WorkingSpaceDimension();

	vector<double> Bu(number_of_nodes * dimension);

	for (IndexType node = 0; node < number_of_nodes; node++)
	{
		for (IndexType dim = 0; dim < dimension; dim++)
		{
			Bu(dimension * node + dim) = rDN_DX(node, dim);
		}
	}

	return Bu;
}

/***********************************************************************************/
/***********************************************************************************/

 Matrix PorousElement::CalculateBMatrix(
     const Matrix& rDN_DX)
 {
     const auto& r_geometry = GetGeometry();
     const SizeType number_of_nodes = r_geometry.PointsNumber();
     const SizeType dimension = r_geometry.WorkingSpaceDimension();

	 Matrix B;
	 B = ZeroMatrix(3 * dimension - 3, dimension * number_of_nodes);

     if(dimension == 2) {
         for ( IndexType i = 0; i < number_of_nodes; ++i ) {
             B(0, 2*i    ) = rDN_DX(i, 0);
             B(1, 2*i + 1) = rDN_DX(i, 1);
             B(2, 2*i    ) = rDN_DX(i, 1);
			 B(2, 2*i + 1) = rDN_DX(i, 0);
         }
     } else if(dimension == 3) {
         for ( IndexType i = 0; i < number_of_nodes; ++i ) {
             const IndexType initial_index = i*3;
             B(0, initial_index    ) = rDN_DX(i, 0);
             B(1, initial_index + 1) = rDN_DX(i, 1);
             B(2, initial_index + 2) = rDN_DX(i, 2);
             B(3, initial_index    ) = rDN_DX(i, 1);
             B(3, initial_index + 1) = rDN_DX(i, 0);
             B(4, initial_index + 1) = rDN_DX(i, 2);
             B(4, initial_index + 2) = rDN_DX(i, 1);
             B(5, initial_index    ) = rDN_DX(i, 2);
             B(5, initial_index + 2) = rDN_DX(i, 0);
         }
     }

	 return B;
 }

/***********************************************************************************/
/***********************************************************************************/

// void AcousticElement::CalculateKinematicVariables(
//     KinematicVariables& rThisKinematicVariables,
//     const IndexType PointNumber,
//     const GeometryType::IntegrationMethod& rIntegrationMethod
//     )
// {
//     const auto& r_geometry = GetGeometry();
//     const SizeType number_of_nodes = r_geometry.PointsNumber();
//     const SizeType dimension = r_geometry.WorkingSpaceDimension();

//     const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(rIntegrationMethod);
//     // Shape functions
//     rThisKinematicVariables.N = r_geometry.ShapeFunctionsValues(Vector &rResult, const CoordinatesArrayType& rCoordinates);
//     rThisKinematicVariables.detJ0 = CalculateDerivativesOnReferenceConfiguration(rThisKinematicVariables.J0, rThisKinematicVariables.InvJ0, rThisKinematicVariables.DN_DX, PointNumber, rIntegrationMethod);

//     KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0) << "WARNING:: ELEMENT ID: " << this->Id() << " INVERTED. DETJ0: " << rThisKinematicVariables.detJ0 << std::endl;

//     // Compute B

//     rThisKinematicVariables.B = rThisKinematicVariables.DN_DX

// }

/***********************************************************************************/
/***********************************************************************************/

void PorousElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
        // (Matrix& rJ0,
        // Matrix& rInvJ0,
        // Matrix& rDN_DX,
        // const IndexType PointNumber,
        // IntegrationMethod ThisIntegrationMethod
        // ) const

		// obtain parameters from geometry
		const GeometryType& geom = GetGeometry();
		IntegrationMethod ThisIntegrationMethod = geom.GetDefaultIntegrationMethod();
		const SizeType number_of_nodes = geom.PointsNumber();
		const SizeType dimension = geom.WorkingSpaceDimension();
		const GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints(ThisIntegrationMethod);

		// get material parameters
		const double lamda_solid = GetProperties()[LAMBDA_SOLID];
		const double mue_solid = GetProperties()[MUE_SOLID];
		// const double damping_solid = GetProperties()[DAMPING_SOLID];
		// const double density_solid = GetProperties()[DENSITY_SOLID];
		// const double density_fluid = GetProperties()[DENSITY_FLUID];
		// const double standard_pressure_fluid = GetProperties()[STANDARD_PRESSURE_FLUID];
		const double porosity = GetProperties()[POROSITY];
		// const double tortuosity = GetProperties()[TORTUOSITY];

		// Determine relevant densities
		// double APPARENT_MASS_DENSITY = porosity * density_fluid * (tortuosity - 1.0);

		// resize leftHandSideMatrix if the size is incorrect
		if (rLeftHandSideMatrix.size1() != number_of_nodes * (dimension + 1)  || rLeftHandSideMatrix.size2() != number_of_nodes * (dimension + 1))
        {
            rLeftHandSideMatrix.resize(number_of_nodes * (dimension + 1), number_of_nodes * (dimension + 1), false);
        }
		noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_nodes * (dimension + 1), number_of_nodes * (dimension + 1));

		// initialize result for shape function derivatives and determinant of the Jacobi Matrix
        ShapeFunctionDerivativesArrayType DN_DX;
        Vector DetJ;

		// get shape functions and their derivatives
		const Matrix& NContainer = geom.ShapeFunctionsValues(ThisIntegrationMethod);
        DN_DX = geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, DetJ, ThisIntegrationMethod);

		// compute material matrix
		Matrix D;

		if (dimension == 2)
		{
			D = ZeroMatrix(3, 3);
			D(0, 0) = lamda_solid + 2.0 * mue_solid;
			D(0, 1) = lamda_solid;
			D(1, 0) = lamda_solid;
			D(1, 1) = lamda_solid + 2.0 * mue_solid;
			D(2, 2) = mue_solid;
		}
		else if (dimension == 3)
		{
			D = ZeroMatrix(6, 6);
			D(0, 0) = lamda_solid + 2.0 * mue_solid;
			D(0, 1) = lamda_solid;
			D(0, 2) = lamda_solid;
			D(1, 0) = lamda_solid;
			D(1, 1) = lamda_solid + 2.0 * mue_solid;
			D(1, 2) = lamda_solid;
			D(2, 0) = lamda_solid;
			D(2, 1) = lamda_solid;
			D(2, 2) = lamda_solid + 2.0 * mue_solid;
			D(3, 3) = mue_solid;
			D(4, 4) = mue_solid;
			D(5, 5) = mue_solid;
		}


		// initialize results for submatrices

		// NPE 1 stiffness
		Matrix Ks_Matrix = ZeroMatrix(number_of_nodes * dimension, number_of_nodes * dimension);
		Matrix C1_Matrix = ZeroMatrix(number_of_nodes * dimension, number_of_nodes);
		Matrix C2_Matrix = ZeroMatrix(number_of_nodes * dimension, number_of_nodes);

		// NPE 2 stiffness
		Matrix C_Matrix = ZeroMatrix(number_of_nodes * dimension, number_of_nodes);

		// NPE 3 stiffness
		Matrix Kf_Matrix = ZeroMatrix(number_of_nodes, number_of_nodes);

        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number )
        {
			// compute coordinate transformation factor (including Gauss weight)
            double int_weight = integration_points[point_number].Weight() * DetJ(point_number);

			// get Be-Matrix, N_mat Matrix and Bu-Matrix
			Matrix Be = CalculateBMatrix(DN_DX[point_number]);
			Matrix N_mat = CalculateNMatrix(row(NContainer, point_number));
			Vector Bu = CalculateBuMatrix(DN_DX[point_number]);

			// compute submatrices

			// NPE 1 - stiffness
			noalias(Ks_Matrix) += int_weight * prod(trans(Be), Matrix(prod(D, Be)));
			noalias(C1_Matrix) += int_weight * porosity * prod(trans(N_mat), trans(DN_DX[point_number]));
			noalias(C2_Matrix) += int_weight * outer_prod(trans(Bu), row(NContainer, point_number));

			// NPE 2 - stiffness
			noalias(C_Matrix) += int_weight * prod(trans(N_mat), trans(DN_DX[point_number]));

			// NPE 3 stiffness
			noalias(Kf_Matrix) += int_weight * prod(DN_DX[point_number], trans(DN_DX[point_number]));

			// NPE 4 stiffness
			// empty
        }


		if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 81)
		{
			// K*1
			noalias(project(rLeftHandSideMatrix, range(0, number_of_nodes * dimension),
				range(number_of_nodes * dimension, number_of_nodes * (dimension + 1)))) = C_Matrix;
		}
		else if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 82)
		{
			// K*2
			noalias(project(rLeftHandSideMatrix, range(number_of_nodes * dimension, number_of_nodes * (dimension + 1)),
				range(number_of_nodes * dimension, number_of_nodes * (dimension + 1)))) = Kf_Matrix;

		}
		else
		{

			// K^ (real part)
			noalias(project(rLeftHandSideMatrix, range(0, number_of_nodes * dimension), range(0, number_of_nodes * dimension))) = Ks_Matrix;
			noalias(project(rLeftHandSideMatrix, range(0, number_of_nodes * dimension),
				range(number_of_nodes * dimension, number_of_nodes * (dimension + 1)))) = - C1_Matrix - C2_Matrix;


		}

		// KRATOS_WATCH(rCurrentProcessInfo[BUILD_LEVEL])
		// KRATOS_WATCH(rLeftHandSideMatrix)


        // GeometryUtils::JacobianOnInitialConfiguration(r_geom, r_geom.IntegrationPoints(ThisIntegrationMethod)[PointNumber], rJ0);
        // double detJ0;
        // MathUtils<double>::InvertMatrix(rJ0, rInvJ0, detJ0);
        // const Matrix& rDN_De =
        // GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[PointNumber];
        // GeometryUtils::ShapeFunctionsGradients(rDN_De, rInvJ0, rDN_DX);
        // return detJ0;
}

/***********************************************************************************/
/***********************************************************************************/

void PorousElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
	// obtain parameters from geometry
	const GeometryType& geom = GetGeometry();
	IntegrationMethod ThisIntegrationMethod = geom.GetDefaultIntegrationMethod();
	const SizeType number_of_nodes = geom.PointsNumber();
	const SizeType dimension = geom.WorkingSpaceDimension();
	const GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints(ThisIntegrationMethod);

	// get material parameters
	// const double lamda_solid = GetProperties()[LAMBDA_SOLID];
	// const double mue_solid = GetProperties()[MUE_SOLID];
	// const double damping_solid = GetProperties()[DAMPING_SOLID];
	const double density_solid = GetProperties()[DENSITY_SOLID];
	const double density_fluid = GetProperties()[DENSITY_FLUID];
	const double standard_pressure_fluid = GetProperties()[STANDARD_PRESSURE_FLUID];
	const double porosity = GetProperties()[POROSITY];
	const double tortuosity = GetProperties()[TORTUOSITY];

	// Determine relevant densities
	double APPARENT_MASS_DENSITY = porosity * density_fluid * (tortuosity - 1.0);

	// resize leftHandSideMatrix if the size is incorrect
	if (rMassMatrix.size1() != number_of_nodes * (dimension + 1) || rMassMatrix.size2() != number_of_nodes * (dimension + 1))
	{
		rMassMatrix.resize(number_of_nodes * (dimension + 1), number_of_nodes * (dimension + 1), false);
	}
	noalias(rMassMatrix) = ZeroMatrix(number_of_nodes * (dimension + 1), number_of_nodes * (dimension + 1));

	// initialize result for shape function derivatives and determinant of the Jacobi Matrix
	ShapeFunctionDerivativesArrayType DN_DX;
	Vector DetJ;

	// get shape functions and their derivatives
	const Matrix& NContainer = geom.ShapeFunctionsValues(ThisIntegrationMethod);
	DN_DX = geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, DetJ, ThisIntegrationMethod);


	// initialize results for submatrices

	// NPE 1 mass
	Matrix Ms1_Matrix = ZeroMatrix(number_of_nodes * dimension, number_of_nodes * dimension);
	Matrix Mf1_Matrix = ZeroMatrix(number_of_nodes, number_of_nodes);
	Matrix C1_Matrix = ZeroMatrix(number_of_nodes * dimension, number_of_nodes);
	Matrix C2_Matrix = ZeroMatrix(number_of_nodes * dimension, number_of_nodes);

	// NPE 2 mass
	Matrix C_Matrix = ZeroMatrix(number_of_nodes * dimension, number_of_nodes);

	// NPE 3 mass
	Matrix Ms3_Matrix = ZeroMatrix(number_of_nodes * dimension, number_of_nodes * dimension);

	// NPE 4
	Matrix Mf4_Matrix = ZeroMatrix(number_of_nodes, number_of_nodes);

	for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
	{
		// compute coordinate transformation factor (including Gauss weight)
		double int_weight = integration_points[point_number].Weight() * DetJ(point_number);

		// get Be-Matrix, N_mat Matrix and Bu-Matrix
		Matrix Be = CalculateBMatrix(DN_DX[point_number]);
		Matrix N_mat = CalculateNMatrix(row(NContainer, point_number));
		Vector Bu = CalculateBuMatrix(DN_DX[point_number]);

		// compute submatrices

		// NPE 1 - mass
		noalias(Ms1_Matrix) += int_weight * ((1.0 - porosity) * density_solid + APPARENT_MASS_DENSITY) * prod(trans(N_mat), N_mat);
		noalias(Mf1_Matrix) += int_weight * (pow(porosity, 2) / standard_pressure_fluid) *
			outer_prod(row(NContainer, point_number), row(NContainer, point_number));
		noalias(C1_Matrix) += int_weight * porosity * prod(trans(N_mat), trans(DN_DX[point_number]));
		noalias(C2_Matrix) += int_weight * outer_prod(trans(Bu), row(NContainer, point_number));

		// NPE 2 - mass
		noalias(C_Matrix) += int_weight * prod(trans(N_mat), trans(DN_DX[point_number]));

		// NPE 3 mass
		noalias(Ms3_Matrix) += int_weight * prod(trans(N_mat), N_mat);

		// NPE 4 mass
		noalias(Mf4_Matrix) += int_weight * outer_prod(row(NContainer, point_number), row(NContainer, point_number));
	}

	if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 281)
	{
		// M*1
		noalias(project(rMassMatrix, range(number_of_nodes * dimension, number_of_nodes * (dimension + 1)),
			range(0, number_of_nodes * dimension))) = -trans(C_Matrix);
	}
	else if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 282)
	{
		// M*2
		noalias(project(rMassMatrix, range(0, number_of_nodes * dimension), range(0, number_of_nodes * dimension))) = Ms3_Matrix;
	}
	else if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 283)
	{
		// M*3
		noalias(project(rMassMatrix, range(number_of_nodes * dimension, number_of_nodes * (dimension + 1)),
			range(number_of_nodes * dimension, number_of_nodes * (dimension + 1)))) = Mf4_Matrix;
	}
	else
	{
		// M^
		noalias(project(rMassMatrix, range(0, number_of_nodes * dimension), range(0, number_of_nodes * dimension))) = real(Ms1_Matrix);
		noalias(project(rMassMatrix, range(number_of_nodes * dimension, number_of_nodes * (dimension + 1)),
			range(0, number_of_nodes * dimension))) = trans(C1_Matrix) + trans(C2_Matrix);
		noalias(project(rMassMatrix, range(number_of_nodes * dimension, number_of_nodes * (dimension + 1)),
			range(number_of_nodes * dimension, number_of_nodes * (dimension + 1)))) = Mf1_Matrix;
	}

	// KRATOS_WATCH(rCurrentProcessInfo[BUILD_LEVEL])
	// KRATOS_WATCH(rMassMatrix)

}


/***********************************************************************************/
void PorousElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{

	const SizeType number_of_nodes = GetGeometry().size();
	const SizeType dimension = GetGeometry().WorkingSpaceDimension();
	const SizeType dofs_size = number_of_nodes * (dimension + 1);

    if( rDampingMatrix.size1() != dofs_size || rDampingMatrix.size2() != dofs_size ) {
        rDampingMatrix.resize(dofs_size, dofs_size, false);
    }
	noalias(rDampingMatrix) = ZeroMatrix( dofs_size, dofs_size );

}

/***********************************************************************************/

void PorousElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{

	const SizeType number_of_nodes = GetGeometry().size();
	const SizeType dimension = GetGeometry().WorkingSpaceDimension();
	const SizeType dofs_size = number_of_nodes * (dimension + 1);

    // Resizing as needed the RHS
	if ( rRightHandSideVector.size() != dofs_size ) {
		rRightHandSideVector.resize( dofs_size, false );
	}
	rRightHandSideVector = ZeroVector( dofs_size ); //resetting RHS
    // }


}


/***********************************************************************************/
/***********************************************************************************/

void PorousElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH( "" )

}


/***********************************************************************************/


void PorousElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
}

/***********************************************************************************/
/***********************************************************************************/

void PorousElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
}

} // Namespace Kratos


