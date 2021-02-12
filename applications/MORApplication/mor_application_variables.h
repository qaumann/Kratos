//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//

#if !defined(KRATOS_MOR_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_MOR_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"

namespace Kratos
{
typedef Node<3> NodeType;
typedef NodeType::Pointer NodeTypePointer;

KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, FREQUENCY )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, ADMITTANCE )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, ACOUSTIC_LOAD )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, SYSTEM_DAMPING_RATIO )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, int, BUILD_LEVEL )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, bool, USE_FREQUENCY_DEPENDENT_MATERIAL)

KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, SCALAR_OUTPUT )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( MOR_APPLICATION, COMPONENT_OUTPUT )

KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, Vector, EIGENVALUE_VECTOR)
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, Matrix, EIGENVECTOR_MATRIX)
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, Vector, REAL_EIGENVALUE_VECTOR)
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, Matrix, REAL_EIGENVECTOR_MATRIX)
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, Vector, IMAG_EIGENVALUE_VECTOR)
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, Matrix, IMAG_EIGENVECTOR_MATRIX)

KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, std::vector<NodeTypePointer>, MAPPING_NODES)
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, Matrix, MAPPING_FACTOR)

// complex dof variables
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( MOR_APPLICATION, REAL_DISPLACEMENT )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( MOR_APPLICATION, IMAG_DISPLACEMENT )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, REAL_PRESSURE )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, IMAG_PRESSURE )

// Biot theory variables
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, LAMBDA_SOLID )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, MUE_SOLID )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, DAMPING_SOLID )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, DENSITY_SOLID )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, DENSITY_FLUID )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, VISCOSITY_FLUID )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, STANDARD_PRESSURE_FLUID )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, HEAT_CAPACITY_FLUID )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, PRANDTL_NUMBER_FLUID )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, THERMAL_LENGTH )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, TORTUOSITY )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, FLOW_RESISTIVITY )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, VISCOUS_LENGTH )

// PML variables
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, int, PARAMETER_EXPONENT )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, PARAMETER_SPATIAL )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, PRESCRIBED_POTENTIAL )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, IMAG_DISTANCE )
KRATOS_DEFINE_APPLICATION_VARIABLE( MOR_APPLICATION, double, LOCAL_PML_WIDTH )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( MOR_APPLICATION, ABSORBTION_VECTOR )
}

#endif	/* KRATOS_MOR_APPLICATION_VARIABLES_H_INCLUDED */
