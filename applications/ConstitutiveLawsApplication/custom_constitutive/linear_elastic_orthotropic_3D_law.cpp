// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Peter Wilson
//       Contact:    A.Winterstein [at] tum.de
//

// System includes
#include <iostream>

// External includes
// #include<cmath>

// Project includes
#include "includes/checks.h"
#include "custom_constitutive/linear_elastic_orthotropic_3D_law.h"

#include "structural_mechanics_application_variables.h"

namespace Kratos
{
	//******************************CONSTRUCTOR*********************************
	//**************************************************************************

	LinearElasticOrthotropic3DLaw::LinearElasticOrthotropic3DLaw()
		: ConstitutiveLaw()
	{
	}

	//******************************COPY CONSTRUCTOR****************************
	//**************************************************************************

	LinearElasticOrthotropic3DLaw::LinearElasticOrthotropic3DLaw
		(const LinearElasticOrthotropic3DLaw& rOther)
		: ConstitutiveLaw(rOther)
	{
	}

	//********************************CLONE*************************************
	//**************************************************************************

	ConstitutiveLaw::Pointer LinearElasticOrthotropic3DLaw::Clone() const
	{
		LinearElasticOrthotropic3DLaw::Pointer p_clone
			(new LinearElasticOrthotropic3DLaw(*this));
		return p_clone;
	}

	//*******************************DESTRUCTOR*********************************
	//**************************************************************************

	LinearElasticOrthotropic3DLaw::~LinearElasticOrthotropic3DLaw()
	{
	}

	//*****************************MATERIAL RESPONSES***************************
	//**************************************************************************

	void  LinearElasticOrthotropic3DLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
	{
        KRATOS_TRY;
        // // 1.- Lame constants
        // // const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
        // // const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];

        // //a.-Check if the constitutive parameters are passed correctly to the law calculation
		// //CheckParameters(rValues);

		// //b.- Get Values to compute the constitutive law:
		// Flags &Options = rValues.GetOptions();

		// const Properties& MaterialProperties = rValues.GetMaterialProperties();

		// Vector& StrainVector = rValues.GetStrainVector();
		// // Options.Set(COMPUTE_STRESS, false);

		// // KRATOS_WATCH(Options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN))
		// // KRATOS_WATCH(Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
		// // KRATOS_WATCH(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))

		// //-----------------------------//

		// if (Options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN))
		// {
		// 	//only needed
		// 	const Matrix& DeformationGradientF = rValues.GetDeformationGradientF();
		// 	KRATOS_WATCH(DeformationGradientF)

		// 	//4.-Right Cauchy Green
		// 	Matrix RightCauchyGreen = prod(trans(DeformationGradientF), DeformationGradientF);

		// 	//5.-Green-Lagrange Strain:

		// 	//E= 0.5*(FT*F-1)
		// 	this->CalculateGreenLagrangeStrain(RightCauchyGreen, StrainVector);
		// }
		// // std::cout << "end" << std::endl;

		// //7.-Calculate Total PK2 stress

		// if (Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
		// {
		// 	Vector& StressVector = rValues.GetStressVector();
		// 	if (Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
		// 	{
		// 		Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
		// 		this->CalculateLinearElasticMatrix(ConstitutiveMatrix, MaterialProperties);
		// 		this->CalculateStress(StrainVector, ConstitutiveMatrix, StressVector);
		// 	}
		// 	else {
		// 		Matrix ConstitutiveMatrix(StrainVector.size(), StrainVector.size());
		// 		noalias(ConstitutiveMatrix) = ZeroMatrix(StrainVector.size(), StrainVector.size());

		// 		this->CalculateLinearElasticMatrix(ConstitutiveMatrix, MaterialProperties);
		// 		this->CalculateStress(StrainVector, ConstitutiveMatrix, StressVector);
		// 	}
		// }
		// else if (Options.IsNot(ConstitutiveLaw::COMPUTE_STRESS) && Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
		// {
		// 	// std::cout << "hier?" << std::endl;
		// 	Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
		// 	// KRATOS_WATCH(ConstitutiveMatrix)
		// 	this->CalculateLinearElasticMatrix(ConstitutiveMatrix, MaterialProperties);
		// }

		// Get the constitutive law options
		Flags & r_constitutive_law_options = rValues.GetOptions();

		Vector& r_strain_vector = rValues.GetStrainVector();

		//NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
		if(r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
			CalculateGreenLagrangeStrainVector(rValues, r_strain_vector);
		}

		if( r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
			Vector& r_stress_vector = rValues.GetStressVector();
			CalculatePK2Stress(r_strain_vector, r_stress_vector, rValues);
		}

		if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR )) {
			Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
			CalculateLinearElasticMatrix(r_constitutive_matrix, rValues.GetMaterialProperties());
		}


        KRATOS_CATCH("");
    }

	/***********************************************************************************/
	/***********************************************************************************/

	void LinearElasticOrthotropic3DLaw::CalculatePK2Stress(
		const Vector& rStrainVector,
		Vector& rStressVector,
		ConstitutiveLaw::Parameters& rValues
		)
	{
		const Properties& r_material_properties = rValues.GetMaterialProperties();
		const Matrix C = r_material_properties[ELASTICITY_TENSOR];
		noalias(rStressVector) = prod(C, rStrainVector);
	}

	/***********************************************************************************/
	/***********************************************************************************/

	void LinearElasticOrthotropic3DLaw::CalculateGreenLagrangeStrainVector(
		ConstitutiveLaw::Parameters& rValues,
		Vector& rStrainVector
		)
	{
		const SizeType dim = this->WorkingSpaceDimension();
		const SizeType strain_size = GetStrainSize();

		// Get the deformation gradient tensor F
		const Matrix& rF = rValues.GetDeformationGradientF();
		KRATOS_DEBUG_ERROR_IF(rF.size1()!= dim || rF.size2() != dim) << "expected size of F " << dim << "x" << dim << ", got " << rF.size1() << "x" << rF.size2() << std::endl;

		// Calculate the Cauchy - Green strain tensor
		Matrix left_cauchy_green = prod(trans(rF), rF);

		// Calculate Green - Lagrange strain tensor
		// ConstitutiveLawUtilities<6>::CalculateGreenLagrangianStrain(left_cauchy_green, rStrainVector);

		// Doing resize in case is needed
		if (rStrainVector.size() != strain_size)
			rStrainVector.resize(strain_size, false);

		// Identity matrix
		Matrix identity_matrix(dim, dim);
		for (IndexType i = 0; i < dim; ++i) {
			for (IndexType j = 0; j < dim; ++j) {
				if (i == j) identity_matrix(i, j) = 1.0;
				else identity_matrix(i, j) = 0.0;
			}
		}

		// Calculate E matrix
		const BoundedMatrix<double, 3, 3> E_matrix = 0.5 * (left_cauchy_green - identity_matrix);

		// Green-Lagrangian Strain Calculation
		rStrainVector = MathUtils<double>::StrainTensorToVector(E_matrix, strain_size);
	}

    //************************************************************************************
	//************************************************************************************

	bool& LinearElasticOrthotropic3DLaw::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
	{
		// This Constitutive Law has been checked with Stenberg Stabilization
		if (rThisVariable == STENBERG_SHEAR_STABILIZATION_SUITABLE)
			rValue = true;

		return rValue;
	}

	//***********************COMPUTE TOTAL STRAIN*****************************************
	//************************************************************************************

	void LinearElasticOrthotropic3DLaw::CalculateGreenLagrangeStrain(const Matrix & rRightCauchyGreen,
		Vector& rStrainVector)
	{
		// TAKEN FROM LinearElasticPlasticPlaneStrain2DLaw
		//E= 0.5*(FT*F-1)
		rStrainVector[0] = 0.5 * (rRightCauchyGreen(0, 0) - 1.00);
		rStrainVector[1] = 0.5 * (rRightCauchyGreen(1, 1) - 1.00);
		rStrainVector[2] = rRightCauchyGreen(0, 1);
	}

	//***********************COMPUTE TOTAL STRESS PK2*************************************
	//************************************************************************************

	void LinearElasticOrthotropic3DLaw::CalculateStress(const Vector & rStrainVector,
		const Matrix & rConstitutiveMatrix,
		Vector& rStressVector)
	{
		//1.-2nd Piola Kirchhoff StressVector increment
		if (rStressVector.size() != rStrainVector.size())
			rStressVector.resize(rStrainVector.size(), false);

		noalias(rStressVector) = prod(rConstitutiveMatrix, rStrainVector);
	}

	//***********************COMPUTE LINEAR ELASTIC MATRIX**********************
	//**************************************************************************

	void LinearElasticOrthotropic3DLaw::CalculateLinearElasticMatrix(Matrix& rConstitutiveMatrix,
		const Properties& rMaterialProperties)
	{
		//double G13 = G12;	// currently handled through "shell_cross_section.cpp"
		//double G23 = G12;	// currently handled through "shell_cross_section.cpp"

        double youngs_modulus_x, youngs_modulus_y, youngs_modulus_z;
		double poisson_ratio_xy, poisson_ratio_xz, poisson_ratio_yz;
		double shear_modulus_xy, shear_modulus_xz, shear_modulus_yz;
        if (rMaterialProperties.Has(SHELL_ORTHOTROPIC_LAYERS))
        {
            // Using the Values directly from the ply-definition
            youngs_modulus_x = rMaterialProperties[SHELL_ORTHOTROPIC_LAYERS](0,1);
            youngs_modulus_y = youngs_modulus_x;
			youngs_modulus_z = rMaterialProperties[SHELL_ORTHOTROPIC_LAYERS](0,2);
            poisson_ratio_xy = rMaterialProperties[SHELL_ORTHOTROPIC_LAYERS](0,3);
            poisson_ratio_xz = 0.0;
            poisson_ratio_yz = 0.0;
            shear_modulus_xy = rMaterialProperties[SHELL_ORTHOTROPIC_LAYERS](0,4);
            shear_modulus_xz = rMaterialProperties[SHELL_ORTHOTROPIC_LAYERS](0,5);
            shear_modulus_yz = rMaterialProperties[SHELL_ORTHOTROPIC_LAYERS](0,6);
        }
        else
        {
            youngs_modulus_x = rMaterialProperties[YOUNG_MODULUS_X];
            youngs_modulus_y = rMaterialProperties[YOUNG_MODULUS_Y];
            youngs_modulus_z = rMaterialProperties[YOUNG_MODULUS_Z];
            poisson_ratio_xy = rMaterialProperties[POISSON_RATIO_XY];
            poisson_ratio_xz = rMaterialProperties[POISSON_RATIO_XZ];
            poisson_ratio_yz = rMaterialProperties[POISSON_RATIO_YZ];
            shear_modulus_xy = rMaterialProperties[SHEAR_MODULUS_XY];
            shear_modulus_xz = rMaterialProperties[SHEAR_MODULUS_XZ];
            shear_modulus_yz = rMaterialProperties[SHEAR_MODULUS_YZ];
        }

		// KRATOS_WATCH(youngs_modulus_x)
		// KRATOS_WATCH(youngs_modulus_y)
		// KRATOS_WATCH(youngs_modulus_z)
		// KRATOS_WATCH(poisson_ratio_xy)
		// KRATOS_WATCH(poisson_ratio_xz)
		// KRATOS_WATCH(poisson_ratio_yz)
		// KRATOS_WATCH(shear_modulus_xy)
		// KRATOS_WATCH(shear_modulus_xz)
		// KRATOS_WATCH(shear_modulus_yz)

		const double v12 = poisson_ratio_xy;
		const double v13 = poisson_ratio_xz;
		const double v23 = poisson_ratio_yz;

		const double v21 = v12*youngs_modulus_y / youngs_modulus_x;
		const double v31 = 0.35;
		const double v32 = 0.35;

		const double Y = 1.0/(1.0 - v12*v21 - v23*v32 - v31*v13 - 2.0*v21*v32*v13);

		const double D1111 = youngs_modulus_x * (1.0 - v23*v32) * Y;
		const double D2222 = youngs_modulus_y * (1.0 - v13*v31) * Y;
		const double D3333 = youngs_modulus_z * (1.0 - v12*v21) * Y;
		const double D1122 = youngs_modulus_x * (v21 + v31*v23) * Y;
		const double D1133 = youngs_modulus_x * (v31 + v21*v32) * Y;
		const double D2233 = youngs_modulus_z * (v32 + v12*v31) * Y;
		const double D1212 = shear_modulus_xy;
		const double D1313 = shear_modulus_xz;
		const double D2323 = shear_modulus_yz;

		rConstitutiveMatrix.clear();

		rConstitutiveMatrix(0,0) = D1111;
		rConstitutiveMatrix(1,1) = D2222;
		rConstitutiveMatrix(2,2) = D3333;
		rConstitutiveMatrix(3,3) = D1212;
		rConstitutiveMatrix(4,4) = D1313;
		rConstitutiveMatrix(5,5) = D2323;
		rConstitutiveMatrix(0,1) = D1122;
		rConstitutiveMatrix(1,0) = D1122;
		rConstitutiveMatrix(0,2) = D1133;
		rConstitutiveMatrix(2,0) = D1133;
		rConstitutiveMatrix(1,2) = D2233;
		rConstitutiveMatrix(2,1) = D2233;

		// const double Q11 = youngs_modulus_x / (1.0 - v12*v21);
		// const double Q12 = v12*youngs_modulus_y / (1.0 - v12*v21);
		// const double Q22 = youngs_modulus_y / (1.0 - v12*v21);
		// const double Q66 = shear_modulus_xy;
		// //double Q44 = G23;
		// //double Q55 = G13;

		// const double theta = 0.0;	// rotation currently handled through
		// // "shell_cross_section.cpp" variable iPlyAngle. Left in for clarity.

		// const double c = std::cos(theta);
		// const double c2 = c*c;
		// const double c4 = c2 * c2;
		// const double s = std::sin(theta);
		// const double s2 = s*s;
		// const double s4 = s2*s2;

		// rConstitutiveMatrix.clear();

		// rConstitutiveMatrix(0, 0) = Q11*c4 + 2.0*(Q12 + 2.0*Q66)*s2*c2 + Q22*s4;				// Q11_hat
		// rConstitutiveMatrix(0, 1) = (Q11 + Q22 - 4.0*Q66)*s2*c2 + Q12*(s4 + c4);				// Q12_hat
		// rConstitutiveMatrix(0, 2) = (Q11 - Q12 - 2.0*Q66)*s*c2*c + (Q12 - Q22 + 2.0*Q66)*s*s2*c;// Q16_hat

		// rConstitutiveMatrix(1, 0) = rConstitutiveMatrix(0, 1);
		// rConstitutiveMatrix(1, 1) = Q11*s4 + 2.0 * (Q12 + 2.0*Q66)*s2*c2 + Q22*c4;				// Q22_hat
		// rConstitutiveMatrix(1, 2) = (Q11 - Q12 - 2.0*Q66)*s2*s*c + (Q12 - Q22 + 2.0*Q66)*c2*c*s;// Q16_hat

		// rConstitutiveMatrix(2, 0) = rConstitutiveMatrix(0, 2);
		// rConstitutiveMatrix(2, 1) = rConstitutiveMatrix(1, 2);
		// rConstitutiveMatrix(2, 2) = (Q11 + Q22 - 2.0*Q12 - 2.0*Q66)*s2*c2 + Q66*(s4 + c4);		//Q66_hat
	}

	//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
	//************************************************************************************

	void LinearElasticOrthotropic3DLaw::GetLawFeatures(Features& rFeatures)
	{
		//Set the type of law
		rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
		rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
		rFeatures.mOptions.Set(ANISOTROPIC);


		//Set strain measure required by the consitutive law
		rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
		rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

		//Set the strain size
		rFeatures.mStrainSize = 6;

		//Set the spacedimension
		rFeatures.mSpaceDimension = 3;
	}


	//******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
	//************************************************************************************

	bool LinearElasticOrthotropic3DLaw::CheckParameters(ConstitutiveLaw::Parameters& rValues)
	{
		return rValues.CheckAllParameters();
	}

	int LinearElasticOrthotropic3DLaw::Check(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const ProcessInfo& rCurrentProcessInfo)
	{
        if(!rMaterialProperties.Has(SHELL_ORTHOTROPIC_LAYERS)) {
            KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS_X));
            KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS_Y));
            KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS_Z));

            KRATOS_CHECK(rMaterialProperties.Has(POISSON_RATIO_XY));
            KRATOS_CHECK(rMaterialProperties.Has(POISSON_RATIO_XZ));
            KRATOS_CHECK(rMaterialProperties.Has(POISSON_RATIO_YZ));

            KRATOS_CHECK(rMaterialProperties.Has(DENSITY));
        }

		return 0;
	}
} // Namespace Kratos