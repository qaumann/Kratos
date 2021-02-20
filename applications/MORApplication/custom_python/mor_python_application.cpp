//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Fullname
//                   Author2 Fullname
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "mor_application.h"
#include "mor_application_variables.h"

#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosMORApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosMORApplication,
        KratosMORApplication::Pointer,
        KratosApplication>(m, "KratosMORApplication")
        .def(py::init<>())
        ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);

    //registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FREQUENCY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ADMITTANCE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ACOUSTIC_LOAD)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BUILD_LEVEL)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, USE_FREQUENCY_DEPENDENT_MATERIAL)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_OUTPUT)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, COMPONENT_OUTPUT)

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, REAL_DISPLACEMENT)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, IMAG_DISPLACEMENT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, REAL_PRESSURE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IMAG_PRESSURE)

    // Biot theory variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LAMBDA_SOLID)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MUE_SOLID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DAMPING_SOLID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DENSITY_SOLID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DENSITY_FLUID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VISCOSITY_FLUID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, STANDARD_PRESSURE_FLUID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HEAT_CAPACITY_FLUID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRANDTL_NUMBER_FLUID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, THERMAL_LENGTH )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TORTUOSITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FLOW_RESISTIVITY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VISCOUS_LENGTH )

    // PML variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PML_PRESCRIBED_POTENTIAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PML_FACTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PML_LOCAL_WIDTH )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, PML_DIRECTION)
}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
