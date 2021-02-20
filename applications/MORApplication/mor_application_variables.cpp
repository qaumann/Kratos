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

#include "mor_application_variables.h"

namespace Kratos
{
KRATOS_CREATE_VARIABLE(double, FREQUENCY)
KRATOS_CREATE_VARIABLE(double, ADMITTANCE)
KRATOS_CREATE_VARIABLE(double, ACOUSTIC_LOAD)
KRATOS_CREATE_VARIABLE(double, SYSTEM_DAMPING_RATIO)
KRATOS_CREATE_VARIABLE(int, BUILD_LEVEL)
KRATOS_CREATE_VARIABLE(bool, USE_FREQUENCY_DEPENDENT_MATERIAL)

KRATOS_CREATE_VARIABLE( double, SCALAR_OUTPUT )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( COMPONENT_OUTPUT)

KRATOS_CREATE_VARIABLE(Vector, EIGENVALUE_VECTOR)
KRATOS_CREATE_VARIABLE(Matrix, EIGENVECTOR_MATRIX)
KRATOS_CREATE_VARIABLE(Vector, REAL_EIGENVALUE_VECTOR)
KRATOS_CREATE_VARIABLE(Matrix, REAL_EIGENVECTOR_MATRIX)
KRATOS_CREATE_VARIABLE(Vector, IMAG_EIGENVALUE_VECTOR)
KRATOS_CREATE_VARIABLE(Matrix, IMAG_EIGENVECTOR_MATRIX)

KRATOS_CREATE_VARIABLE(std::vector<NodeTypePointer>, MAPPING_NODES)
KRATOS_CREATE_VARIABLE(Matrix, MAPPING_FACTOR)

// complex dof variables
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( REAL_DISPLACEMENT )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( IMAG_DISPLACEMENT )
KRATOS_CREATE_VARIABLE( double, REAL_PRESSURE )
KRATOS_CREATE_VARIABLE( double, IMAG_PRESSURE )

// Biot theory variables
KRATOS_CREATE_VARIABLE( double, LAMBDA_SOLID )
KRATOS_CREATE_VARIABLE( double, MUE_SOLID )
KRATOS_CREATE_VARIABLE( double, DAMPING_SOLID )
KRATOS_CREATE_VARIABLE( double, DENSITY_SOLID )
KRATOS_CREATE_VARIABLE( double, DENSITY_FLUID )
KRATOS_CREATE_VARIABLE( double, VISCOSITY_FLUID )
KRATOS_CREATE_VARIABLE( double, STANDARD_PRESSURE_FLUID )
KRATOS_CREATE_VARIABLE( double, HEAT_CAPACITY_FLUID )
KRATOS_CREATE_VARIABLE( double, PRANDTL_NUMBER_FLUID )
KRATOS_CREATE_VARIABLE( double, THERMAL_LENGTH )
KRATOS_CREATE_VARIABLE( double, TORTUOSITY )
KRATOS_CREATE_VARIABLE( double, FLOW_RESISTIVITY )
KRATOS_CREATE_VARIABLE( double, VISCOUS_LENGTH )

// PML variables
KRATOS_CREATE_VARIABLE( double, PML_PRESCRIBED_POTENTIAL )
KRATOS_CREATE_VARIABLE( double, PML_FACTOR )
KRATOS_CREATE_VARIABLE( double, PML_LOCAL_WIDTH )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( PML_DIRECTION )

}
