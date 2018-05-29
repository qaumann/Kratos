// System includes

// External includes

// Project includes
//#include "includes/define.h"
#include "testing/testing.h"
#include "utilities/read_materials_utility.hpp"
#include "geometries/triangle_2d_3.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(ReadMaterialsUtility, KratosCoreFastSuite)
{
    Parameters parameters(R"(
        {
            "materials_filename": "/home/mraschi/Projects/Kratos/kratos/tests/utilities/test_read_materials.json"
        })");

//        {
//        "properties": [{
//            "model_part_name": "MATERIAL_HOMOGENEOUS",
//            "properties_id": 1,
//            "Material": {
//                "name": "homog",
//                "constitutive_law": {
//                    "name": "LinearElastic3DLaw"
//                },
//                "Variables": {
//                    "DENSITY": 7850.0,
//                    "YOUNG_MODULUS": 210000000000.0,
//                    "POISSON_RATIO": 0.3
//                    },
//                 "Tables": {
//                    "Table1": {
//        		    "input_variable": "TEMPERATURE",
//        		    "output_variable": "YOUNG_MODULUS",
//        		    "data": [ [0, 100], [1, 200], [2, 300] ]
//                    },
//        		        "Table2": {
//                        "input_variable": "TEMPERATURE",
//        	                "output_variable": "YOUNG_MODULUS",
//        	                "data": [ [3, 400], [4, 500], [5, 600] ]
//                        }
//        	            }
//                }
//            }, {
//            "model_part_name": "MATERIAL_MULTISCALE",
//            "properties_id": 2,
//            "Material": {
//                "name": "multiscale",
//                "constitutive_law": {
//                    "name": "LinearElastic3DLaw"
//                },
//                "Variables": {
//                    "DENSITY": 7.0,
//                    "YOUNG_MODULUS": 2100.0,
//                    "POISSON_RATIO": 0.3
//                    },
//                "Tables":  {}
//                }
//            }]
//        }

    ModelPart model_part("dummy");
    ReadMaterialsUtility(model_part, parameters);

    KRATOS_CHECK_EQUAL(model_part.NumberOfProperties(), 2);
    Properties& r_property = model_part.GetProperties(1, 0);
	KRATOS_CHECK(r_property.HasVariables());
    KRATOS_CHECK(r_property.HasTables());
    KRATOS_CHECK(r_property.Has(CONSTITUTIVE_LAW));
    KRATOS_CHECK(r_property.HasTable(TEMPERATURE, YOUNG_MODULUS));
    KRATOS_CHECK_EQUAL(r_property(DENSITY), 7850.0);
    KRATOS_CHECK_EQUAL(r_property(CONSTITUTIVE_LAW)->Info(), "ConstitutiveLaw");
}
}
} // namespace Kratos.
