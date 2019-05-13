/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
*/

#include "iga_application_variables.h"

namespace Kratos
{

KRATOS_CREATE_VARIABLE(int, BREP_ID)

KRATOS_CREATE_VARIABLE(double, NURBS_CONTROL_POINT_WEIGHT)

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(COORDINATES)
KRATOS_CREATE_VARIABLE(Vector, LOCAL_COORDINATES)
KRATOS_CREATE_VARIABLE(Vector, TANGENTS)
KRATOS_CREATE_VARIABLE(Vector, TANGENTS_SLAVE)

KRATOS_CREATE_VARIABLE(double, CROSS_AREA)
KRATOS_CREATE_VARIABLE(double, PRESTRESS_CAUCHY)

KRATOS_CREATE_VARIABLE(Vector, SHAPE_FUNCTION_VALUES)
KRATOS_CREATE_VARIABLE(Matrix, SHAPE_FUNCTION_LOCAL_DERIVATIVES)
KRATOS_CREATE_VARIABLE(Matrix, SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES)
KRATOS_CREATE_VARIABLE(Matrix, SHAPE_FUNCTION_LOCAL_THIRD_DERIVATIVES)

KRATOS_CREATE_VARIABLE(Vector, SHAPE_FUNCTION_VALUES_SLAVE)
KRATOS_CREATE_VARIABLE(Matrix, SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE)
KRATOS_CREATE_VARIABLE(Matrix, SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES_SLAVE)

KRATOS_CREATE_VARIABLE(double, RAYLEIGH_ALPHA)
KRATOS_CREATE_VARIABLE(double, RAYLEIGH_BETA)


//Condition Variables
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(POINT_LOAD)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(LINE_LOAD)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SURFACE_LOAD)

KRATOS_CREATE_VARIABLE(double, PENALTY_FACTOR)

// Stress recovery
KRATOS_CREATE_VARIABLE(double, STRESS_CAUCHY_11)
KRATOS_CREATE_VARIABLE(double, STRESS_CAUCHY_12)
KRATOS_CREATE_VARIABLE(double, STRESS_CAUCHY_22)
KRATOS_CREATE_VARIABLE(double, STRESS_CAUCHY_TOP_11)
KRATOS_CREATE_VARIABLE(double, STRESS_CAUCHY_TOP_12)
KRATOS_CREATE_VARIABLE(double, STRESS_CAUCHY_TOP_22)
KRATOS_CREATE_VARIABLE(double, STRESS_CAUCHY_BOTTOM_11)
KRATOS_CREATE_VARIABLE(double, STRESS_CAUCHY_BOTTOM_12)
KRATOS_CREATE_VARIABLE(double, STRESS_CAUCHY_BOTTOM_22)
KRATOS_CREATE_VARIABLE(double, INTERNAL_FORCE_11)
KRATOS_CREATE_VARIABLE(double, INTERNAL_FORCE_12)
KRATOS_CREATE_VARIABLE(double, INTERNAL_FORCE_22)
KRATOS_CREATE_VARIABLE(double, INTERNAL_MOMENT_11)
KRATOS_CREATE_VARIABLE(double, INTERNAL_MOMENT_12)
KRATOS_CREATE_VARIABLE(double, INTERNAL_MOMENT_22)
KRATOS_CREATE_VARIABLE(double, SHEAR_FORCE_1)
KRATOS_CREATE_VARIABLE(double, SHEAR_FORCE_2)

} // namespace Kratos
