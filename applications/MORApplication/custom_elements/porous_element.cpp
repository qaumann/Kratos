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
#include "includes/ublas_complex_interface.h"

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

    // // The vector containing the constitutive laws
    // p_new_elem->SetConstitutiveLawVector(BaseType::mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

int  PorousElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    int ier = Element::Check(rCurrentProcessInfo);

    return ier;

    KRATOS_CATCH( "" );
}

void PorousElement::Initialize()
{
    KRATOS_TRY
    std::cout << "i am initializing an porous element\n";
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

    if(rResult.size() != dofs_size)
        rResult.resize(dofs_size,false);

	SizeType index;
	// order of the dofs: displacement, pressure
	if (dimension == 2) {
		for (SizeType i = 0; i < number_of_nodes; i++) {
			index = i * 3;
			rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
			rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
			rResult[index + 2] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
		}
	}
	else {
		for (SizeType i = 0; i < number_of_nodes; i++) {
			index = i * 4;
			rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
			rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
			rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
			rResult[index + 3] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
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

	if (dimension == 2) {
		for (SizeType i = 0; i < GetGeometry().size(); i++)
		{

			rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
			rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
			rElementalDofList.push_back(GetGeometry()[i].pGetDof(PRESSURE));
		}
	}
	else {
		for (SizeType i = 0; i < GetGeometry().size(); i++)
		{

			rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
			rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
			rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
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


		/*const double LAMBDA_SOLID = GetProperties()[LAMBDA_SOLID];
		const double MUE_SOLID = GetProperties()[MUE_SOLID];
		const double DAMPING_SOLID = GetProperties()[DAMPING_SOLID];
		const double DENSITY_SOLID = GetProperties()[DENSITY_SOLID];
		const double DENSITY_FLUID = GetProperties()[DENSITY_FLUID];
		const double VISCOSITY_FLUID = GetProperties()[VISCOSITY_FLUID];
		const double STANDARD_PRESSURE_FLUID = GetProperties()[STANDARD_PRESSURE_FLUID];
		const double HEAT_CAPACITY_FLUID = GetProperties()[HEAT_CAPACITY_FLUID];
		const double PRANDTL_NUMBER_FLUID = GetProperties()[prandtl_number_fluid];
		const double POROSITY = GetProperties()[POROSITY];
		const double THERMAL_LENGTH = GetProperties()[THERMAL_LENGTH];
		const double OMEGA = GetProperties()[FREQUENCY];
		const double TORTUOSITY = GetProperties()[TORTUOSITY];
		const double FLOW_RESISTIVITY = GetProperties()[FLOW_RESISTIVITY];
		const double VISCOUS_LENGTH = GetProperties()[VISCOUS_LENGTH];
		*/

		// hard coded parameters to check the result
		const double LAMBDA_SOLID = 500.0;
		const double MUE_SOLID = 30.0;
		const double DAMPING_SOLID = 0.2;
		const double DENSITY_SOLID = 8.0;
		const double DENSITY_FLUID = 1.21;
		const double VISCOSITY_FLUID = 1.84e-5;
		const double STANDARD_PRESSURE_FLUID = 101e3;
		const double HEAT_CAPACITY_FLUID = 1.4;
		const double PRANDTL_NUMBER_FLUID = 0.71;
		const double POROSITY = 0.3;
		const double THERMAL_LENGTH = 1.5;
		const double OMEGA = 15.0;
		const double TORTUOSITY = 2.07;
		const double FLOW_RESISTIVITY = 4410;
		const double VISCOUS_LENGTH = 1.1e-2;
		//const int BUILD_LEVEL = 201;


		// Determine relevant densities
		std::complex<double> VISCOUS_DRAG = FLOW_RESISTIVITY * pow(POROSITY, 2.0) * pow(1.0 +
			std::complex<double>(0,4) * OMEGA * pow(TORTUOSITY, 2.0) * VISCOSITY_FLUID * DENSITY_FLUID / ( pow(FLOW_RESISTIVITY, 2.0) * pow(VISCOUS_LENGTH, 2) * pow(POROSITY, 2)), 0.5);
		double APPARENT_MASS_DENSITY = POROSITY * DENSITY_FLUID * (TORTUOSITY - 1.0);
		std::complex<double> EQUIVALENT_COUPLING_DENSITY = - APPARENT_MASS_DENSITY + std::complex<double>(0,1) * VISCOUS_DRAG / OMEGA;
		std::complex<double> EQUIVALENT_FLUID_DENSITY = POROSITY * DENSITY_FLUID - EQUIVALENT_COUPLING_DENSITY;
		std::complex<double> EQUIVALENT_SOLID_DENSITY = (1.0 - POROSITY) * DENSITY_SOLID - EQUIVALENT_COUPLING_DENSITY;
		std::complex<double> EQUIVALENT_DENSITY = EQUIVALENT_SOLID_DENSITY - pow(EQUIVALENT_COUPLING_DENSITY, 2) / EQUIVALENT_FLUID_DENSITY;

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

		std::complex<double> K = HEAT_CAPACITY_FLUID * STANDARD_PRESSURE_FLUID / pow(HEAT_CAPACITY_FLUID -
			(HEAT_CAPACITY_FLUID - 1) * ((1 + 8 * VISCOSITY_FLUID / (std::complex<double>(0,1) * OMEGA * PRANDTL_NUMBER_FLUID * pow(THERMAL_LENGTH, 2) * DENSITY_FLUID) *
			pow(1 + std::complex<double>(0,1) * OMEGA * PRANDTL_NUMBER_FLUID * pow(THERMAL_LENGTH, 2) * DENSITY_FLUID / (16 * VISCOSITY_FLUID), 0.5))), -1);
		std::complex<double> R = POROSITY * K;
		std::complex<double> Q = (1 - POROSITY) * K;
		std::complex<double> GAMMA = POROSITY * ((EQUIVALENT_COUPLING_DENSITY / EQUIVALENT_FLUID_DENSITY) - (Q / R));

		// compute material matrix

		// Determine Coefficients for the material matrix
		std::complex<double> LAMBDA = (1 + std::complex<double>(0,1) * DAMPING_SOLID) * LAMBDA_SOLID;
		std::complex<double> MUE = (1 + std::complex<double>(0,1) * DAMPING_SOLID) * MUE_SOLID;

		// real part
		ComplexMatrix D;

		if (dimension == 2)
		{
			D = ZeroMatrix(3, 3);
			D(0, 0) = LAMBDA + (2.0 + std::complex<double>(0,0)) * MUE;
			D(0, 1) = LAMBDA;
			D(1, 0) = LAMBDA;
			D(1, 1) = LAMBDA + (2.0 + std::complex<double>(0,0)) * MUE;
			D(2, 2) = MUE;
		}

		// initialize results for submatrices


		// NPE 1 stiffness
		ComplexMatrix Ks_Matrix = ZeroMatrix(number_of_nodes * dimension, number_of_nodes * dimension);
		Matrix C1_Matrix = ZeroMatrix(number_of_nodes * dimension, number_of_nodes);
		Matrix C2_Matrix = ZeroMatrix(number_of_nodes * dimension, number_of_nodes);

		// NPE 1 mass
		Matrix Ms1_Matrix = ZeroMatrix(number_of_nodes * dimension, number_of_nodes * dimension);
		Matrix Mf1_Matrix = ZeroMatrix(number_of_nodes, number_of_nodes);

		// NPE 2 stiffness and mass
		Matrix C_Matrix = ZeroMatrix(number_of_nodes * dimension, number_of_nodes);

		// NPE 3 stiffness
		Matrix Kf_Matrix = ZeroMatrix(number_of_nodes, number_of_nodes);

		// NPE 3 mass
		Matrix Ms3_Matrix = ZeroMatrix(number_of_nodes * dimension, number_of_nodes * dimension);

		// NPE 4
		Matrix Mf4_Matrix = ZeroMatrix(number_of_nodes, number_of_nodes);

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
			noalias(Ks_Matrix) += int_weight * prod(trans(Be), ComplexMatrix(prod(D, Be)));
			noalias(C1_Matrix) += int_weight * POROSITY * prod(trans(N_mat), trans(DN_DX[point_number]));
			noalias(C2_Matrix) += int_weight * outer_prod(trans(Bu), row(NContainer, point_number));

			// NPE 1 - mass
			noalias(Ms1_Matrix) += int_weight * ( (1.0 - POROSITY) * DENSITY_SOLID + APPARENT_MASS_DENSITY) * prod(trans(N_mat), N_mat);
			noalias(Mf1_Matrix) += int_weight * (pow(POROSITY, 2) / STANDARD_PRESSURE_FLUID) *
				outer_prod(row(NContainer, point_number), row(NContainer, point_number));
			// C1_Matrix
			// C2_Matrix

			// NPE 2 - stiffness
			noalias(C_Matrix) += int_weight * prod(trans(N_mat), trans(DN_DX[point_number]));

			// NPE 2 - mass
			// C_Matrix

			// NPE 3 stiffness
			noalias(Kf_Matrix) += int_weight * prod(DN_DX[point_number], trans(DN_DX[point_number]));

			// NPE 3 mass
			noalias(Ms3_Matrix) += int_weight * prod(trans(N_mat), N_mat);

			// NPE 4 stiffness
			// empty

			// NPE 4 mass
			noalias(Mf4_Matrix) += int_weight * outer_prod(row(NContainer, point_number), row(NContainer, point_number));
        }

		KRATOS_WATCH(Mf4_Matrix)
		KRATOS_WATCH(rCurrentProcessInfo[BUILD_LEVEL])

		if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 91)
		{
			// real part of NPE 1 stiffness
			noalias(project(rLeftHandSideMatrix, range(0, number_of_nodes * dimension), range(0, number_of_nodes * dimension))) = real(Ks_Matrix);
			noalias(project(rLeftHandSideMatrix, range(0, number_of_nodes * dimension),
				range(number_of_nodes * dimension, number_of_nodes * (dimension + 1)))) = -C1_Matrix;
			noalias(project(rLeftHandSideMatrix, range(0, number_of_nodes * dimension),
				range(number_of_nodes * dimension, number_of_nodes * (dimension + 1)))) = -C2_Matrix;
		}
		else if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 92)
		{
			// imaginary part of NPE 1 stiffness
			noalias(project(rLeftHandSideMatrix, range(0, number_of_nodes * dimension), range(0, number_of_nodes * dimension))) = imag(Ks_Matrix);
			// C1 and C2 do not have an imaginary part
		}
		else if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 93)
		{
			// real part of NPE 2 stiffness
			noalias(project(rLeftHandSideMatrix, range(0, number_of_nodes * dimension),
				range(number_of_nodes * dimension, number_of_nodes * (dimension + 1)))) = C_Matrix;
		}
		else if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 94)
		{
			// imaginary part of NPE 2 stiffness
			// there is no imaginary part of the C_Matrix
		}
		else if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 95)
		{
			// real part of NPE 3 stiffness
			noalias(project(rLeftHandSideMatrix, range(number_of_nodes * dimension, number_of_nodes * (dimension + 1)),
				range(number_of_nodes * dimension, number_of_nodes * (dimension + 1)))) = Kf_Matrix;
		}
		else if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 96)
		{
			// imaginary part of NPE 3 stiffness
			// there is no imaginary part of the Kf_Matrix
		}
		else if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 97)
		{
			// NPE 4 stiffness is empty
		}
		else if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 98)
		{
			// NPE 4 stiffness is empty
		}

		if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 291)
		{
			// real part of NPE 1 mass
			noalias(project(rLeftHandSideMatrix, range(0, number_of_nodes* dimension), range(0, number_of_nodes* dimension))) = real(Ms1_Matrix);
			noalias(project(rLeftHandSideMatrix, range(number_of_nodes* dimension, number_of_nodes* (dimension + 1)),
				range(0, number_of_nodes* dimension))) = trans(C1_Matrix);
			noalias(project(rLeftHandSideMatrix, range(number_of_nodes* dimension, number_of_nodes* (dimension + 1)),
				range(0, number_of_nodes* dimension))) = trans(C2_Matrix);
			noalias(project(rLeftHandSideMatrix, range(number_of_nodes* dimension, number_of_nodes* (dimension + 1)),
				range(number_of_nodes* dimension, number_of_nodes* (dimension + 1)))) = Mf1_Matrix;
		}
		else if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 292)
		{
			// imaginary part of NPE 1 mass
			noalias(project(rLeftHandSideMatrix, range(0, number_of_nodes* dimension), range(0, number_of_nodes* dimension))) = imag(Ms1_Matrix);
			// C1_Matrix and C2_Matrix do not have an imaginary part
		}
		else if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 293)
		{
			// real part of NPE 2 mass
			noalias(project(rLeftHandSideMatrix, range(number_of_nodes* dimension, number_of_nodes* (dimension + 1)),
				range(0, number_of_nodes* dimension))) = - trans(C_Matrix);
		}
		else if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 294)
		{
			// imaginary part of NPE 2 mass
			// there is no imaginary part of the C_Matrix
		}
		else if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 295)
		{
			// real part of NPE 3 mass
			noalias(project(rLeftHandSideMatrix, range(0, number_of_nodes* dimension), range(0, number_of_nodes* dimension))) = real(Ms3_Matrix);
		}
		else if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 296)
		{
			// imaginary part of NPE 3 mass
			// there is no imaginary part of the Ms3_Matrix
		}
		else if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 297)
		{
			// real part of NPE 4 mass
			noalias(project(rLeftHandSideMatrix, range(number_of_nodes* dimension, number_of_nodes* (dimension + 1)),
				range(number_of_nodes* dimension, number_of_nodes* (dimension + 1)))) = Mf4_Matrix;
		}
		else if (rCurrentProcessInfo.Has(BUILD_LEVEL) && rCurrentProcessInfo[BUILD_LEVEL] == 298)
		{
			// imaginary part of the NPE 4 mass
			// there is no imaginary part of the Mf4_Matrix
		}


		KRATOS_WATCH(rLeftHandSideMatrix)


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
	std::cout << "mass matrix\n";
	const SizeType number_of_nodes = GetGeometry().size();
	const SizeType dimension = GetGeometry().WorkingSpaceDimension();

	// 2D -> 3 dofs, 3D -> 4 dofs per node (pressure + displacement)
	const SizeType dofs_size = number_of_nodes * (dimension + 1);


    // const GeometryType& geom = GetGeometry();
    // const SizeType number_of_nodes = geom.PointsNumber();
    if( rMassMatrix.size1() != dofs_size || rMassMatrix.size2() != dofs_size )
    {
        rMassMatrix.resize(dofs_size, dofs_size, false);
        noalias(rMassMatrix) = ZeroMatrix( dofs_size, dofs_size );
    }
    // const GeometryType& geom = GetGeometry();
    // IntegrationMethod ThisIntegrationMethod = geom.GetDefaultIntegrationMethod();
    // const SizeType number_of_nodes = geom.PointsNumber();
    // const GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints(ThisIntegrationMethod);
    // IndexType NumGauss = integration_points.size();
    // const Matrix& NContainer = geom.ShapeFunctionsValues(ThisIntegrationMethod);
    // double Vol = geom.Volume();
    // if( rMassMatrix.size1() != number_of_nodes || rMassMatrix.size2() != number_of_nodes )
    // {
    //     rMassMatrix.resize(number_of_nodes, number_of_nodes, false);
    //     noalias(rMassMatrix) = ZeroMatrix( number_of_nodes, number_of_nodes );
    // }

    // // KRATOS_WATCH(BULK_MODULUS)

    // const double p = GetProperties()[DENSITY];
    // const double G = GetProperties()[BULK_MODULUS];

    // KRATOS_WATCH(NContainer)

    // for (IndexType i = 0; i < number_of_nodes; i++)
    // {
    //     for (IndexType j = 0; j < number_of_nodes; j++)
    //     {
    //         for (IndexType g = 0; g < NumGauss; g++)
    //             {
    //                 double DetJ = geom.DeterminantOfJacobian(g, ThisIntegrationMethod);
    //                 double GaussWeight = DetJ * integration_points[g].Weight();
    //                 rMassMatrix(i,j) += NContainer(g, i) * NContainer(g, j) * GaussWeight * (p/G);
    //             }
    //     }
    // }
    // KRATOS_WATCH(rMassMatrix)
}


/***********************************************************************************/
void PorousElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{

	const SizeType number_of_nodes = GetGeometry().size();
	const SizeType dimension = GetGeometry().WorkingSpaceDimension();
	const SizeType dofs_size = number_of_nodes * (dimension + 1);

    if( rDampingMatrix.size1() != dofs_size || rDampingMatrix.size2() != dofs_size )
    {
        rDampingMatrix.resize(dofs_size, dofs_size, false);
        noalias(rDampingMatrix) = ZeroMatrix( dofs_size, dofs_size );
    }

}

/***********************************************************************************/

void PorousElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{

	const SizeType number_of_nodes = GetGeometry().size();
	const SizeType dimension = GetGeometry().WorkingSpaceDimension();
	const SizeType dofs_size = number_of_nodes * (dimension + 1);

    // Resizing as needed the RHS
        if ( rRightHandSideVector.size() != dofs_size )
            rRightHandSideVector.resize( dofs_size, false );

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

    // auto& r_geometry = this->GetGeometry();
    // const SizeType number_of_nodes = r_geometry.size();
    // const SizeType dimension = r_geometry.WorkingSpaceDimension();
    // const SizeType strain_size = GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();

    // KinematicVariables this_kinematic_variables(number_of_nodes, dimension);

    // const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());



    KRATOS_CATCH( "" )

}


/***********************************************************************************/


// void AcousticElement::save( Serializer& rSerializer ) const
// {
//     KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseSolidElement );
// }

/***********************************************************************************/
/***********************************************************************************/

// void AcousticElement::load( Serializer& rSerializer )
// {
//     KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseSolidElement );
// }

} // Namespace Kratos


