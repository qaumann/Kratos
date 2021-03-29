// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Sebastian Schopper
//

// System includes

// External includes

// Project includes
#include "includes/process_info.h"
#include "testing/testing.h"
#include "containers/model.h"

// Application includes

// Constitutive law
#include "custom_advanced_constitutive/linear_elastic_orthotropic_3D_law.h"
#include "includes/model_part.h"
#include "geometries/tetrahedra_3d_4.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
	namespace Testing
	{
		// We test the associated orthotropic Constitutive law...
		typedef Node<3> NodeType;

		/**
		* Check the correct calculation of the integrated stress with the CL's
		*/

		//KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawOrthotropicInternalVariables, KratosStructuralMechanicsFastSuite)
		//{
		//	//
		//	// Test: check correct behavior of internal and calculated variables
		//	//

		//	ViscousGeneralizedKelvin<ElasticIsotropic3D> cl = ViscousGeneralizedKelvin<ElasticIsotropic3D>();

		//	KRATOS_CHECK_IS_FALSE(cl.Has(INTEGRATED_STRESS_TENSOR));  // = False, in order to use CalculateValue())

		//	// This constitutive law does not use internal variables
		//	// TODO (marandra): check that this is compatible con API
		//	KRATOS_CHECK_IS_FALSE(cl.Has(INTERNAL_VARIABLES));  // = False
		//}

		
		KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawOrthotropic, KratosStructuralMechanicsFastSuite)
		{
			ConstitutiveLaw::Parameters cl_parameters;
			Properties material_properties;
			ProcessInfo process_info;
			Vector stress_vector, strain_vector;

			Model current_model;
			ModelPart& test_model_part = current_model.CreateModelPart("Main");

			NodeType::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
			NodeType::Pointer p_node_2 = test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
			NodeType::Pointer p_node_3 = test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
			NodeType::Pointer p_node_4 = test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);

			Tetrahedra3D4<NodeType> Geom = Tetrahedra3D4<NodeType>(p_node_1, p_node_2, p_node_3, p_node_4);

			stress_vector = ZeroVector(6);
			strain_vector = ZeroVector(6);
			strain_vector[0] = 0.0;
			strain_vector[1] = 0.0;
			strain_vector[2] = -8.0e-5;
			strain_vector[3] = 0.0;
			strain_vector[4] = 0.0;
			strain_vector[5] = -1.6941e-21;

			material_properties.SetValue(YOUNG_MODULUS_X, 1.751147194402183e5);
			material_properties.SetValue(YOUNG_MODULUS_Y, 1.751147194402183e5);
			material_properties.SetValue(YOUNG_MODULUS_Z, 8.312525690438484e8);
			material_properties.SetValue(POISSON_RATIO_XY, 0.9996);
			material_properties.SetValue(POISSON_RATIO_XZ, 7.373228557305470e-5);
			material_properties.SetValue(POISSON_RATIO_YZ, 7.373228557305470e-5);
			material_properties.SetValue(SHEAR_MODULUS_XY, 4.339799568735845e4);
			material_properties.SetValue(SHEAR_MODULUS_XZ, 1.543754771078400e8);
			material_properties.SetValue(SHEAR_MODULUS_YZ, 1.543754771078400e8);
			
			Flags cl_options;
			cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
			cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
			cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

			cl_parameters.SetElementGeometry(Geom);
			cl_parameters.SetMaterialProperties(material_properties);
			cl_parameters.SetStrainVector(strain_vector);
			cl_parameters.SetStressVector(stress_vector);
			cl_parameters.SetOptions(cl_options);
			Matrix const_matrix;
			cl_parameters.SetConstitutiveMatrix(const_matrix);

			// Create the CL
			LinearElasticOrthotropic3DLaw orthotropic_cl = LinearElasticOrthotropic3DLaw();

			std::vector<double> orthotropic_res;
			orthotropic_res = { -1.407402260511011e+04, -1.407402260511011e+04, -7.635202134708494e+04, 0.0, 0.0, -2.544261493143944e-13 };

			Vector test_orthotropic_stress;
			orthotropic_cl.CalculateMaterialResponsePK2(cl_parameters);
			test_orthotropic_stress = cl_parameters.GetStressVector();

			// Check the results
			KRATOS_CHECK_VECTOR_NEAR(test_orthotropic_stress, orthotropic_res, 1.0);
		} 
	} // namespace Testing
} // namespace Kratos
