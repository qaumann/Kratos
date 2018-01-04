//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    msandre
//

// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define_python.h"
#include "includes/ale_variables.h"
#include "python/add_ale_variables_to_python.h"

namespace Kratos
{

namespace Python
{
    using namespace pybind11;

    void  AddALEVariablesToPython(pybind11::module& m)
    {
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,MESH_DISPLACEMENT)
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,MESH_ACCELERATION)
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,MESH_REACTION)
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,MESH_RHS)
    }
}  // namespace Python.
} // Namespace Kratos
