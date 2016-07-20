//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_POROMECHANICS_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_POROMECHANICS_APPLICATION_VARIABLES_H_INCLUDED

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"
#include "solid_mechanics_application_variables.h"

namespace Kratos
{
//Define Variables

//Warning: Note that the application variables must not be defined if they already exist in "includes/variables.h", 
//         in "includes/cfd_variables.h" or in "solid_mechanics_application_variables.h"

KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, double, NEWMARK_COEFFICIENT_U )
KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, double, NEWMARK_COEFFICIENT_P )

KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, double, DT_WATER_PRESSURE )
KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, double, NORMAL_FLUID_FLUX )

KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, double, DENSITY_SOLID )
KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, double, BULK_MODULUS_SOLID )
KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, double, BULK_MODULUS_FLUID )
KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, double, PERMEABILITY_XX )
KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, double, PERMEABILITY_YY )
KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, double, PERMEABILITY_ZZ )
KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, double, PERMEABILITY_XY )
KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, double, PERMEABILITY_YZ )
KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, double, PERMEABILITY_ZX )

KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, double, MINIMUM_JOINT_WIDTH )
KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, double, TRANSVERSAL_PERMEABILITY )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( POROMECHANICS_APPLICATION, FLUID_FLUX_VECTOR )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( POROMECHANICS_APPLICATION, LOCAL_FLUID_FLUX_VECTOR )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( POROMECHANICS_APPLICATION, LOCAL_STRESS_VECTOR )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( POROMECHANICS_APPLICATION, LOCAL_RELATIVE_DISPLACEMENT_VECTOR )
KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, Matrix, PERMEABILITY_MATRIX )
KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, Matrix, LOCAL_PERMEABILITY_MATRIX )

KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, double, CRITICAL_DISPLACEMENT )
KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, double, RESIDUAL_STRESS )

KRATOS_DEFINE_APPLICATION_VARIABLE(POROMECHANICS_APPLICATION, bool, IS_CONVERGED)

KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, Matrix, TOTAL_STRESS_TENSOR )

KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, double, STATE_VARIABLE )
KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, double, ARC_LENGTH_LAMBDA )
KRATOS_DEFINE_APPLICATION_VARIABLE( POROMECHANICS_APPLICATION, double, ARC_LENGTH_RADIUS_FACTOR )
}

#endif	/* KRATOS_POROMECHANICS_APPLICATION_VARIABLES_H_INCLUDED */
