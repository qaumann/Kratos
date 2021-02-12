//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt

// System includes
#include "custom_python/add_custom_processes_to_python.h"

// External includes

// Project includes
// #include "includes/define.h"
// #include "includes/define_python.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

#include "processes/process.h"
#include "custom_processes/monolithic_mapping_process.hpp"
#include "custom_processes/acoustic_pml_direction_process.hpp"

namespace Kratos
{

namespace Python
{

void  AddCustomProcessesToPython(pybind11::module& m)
{

    namespace py = pybind11;
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;
    typedef LinearSolver<SparseSpaceType, DenseSpaceType > LinearSolverType;
    typedef SparseSpaceType::MatrixType SparseMatrixType;

    py::class_<MonolithicMappingProcess<SparseMatrixType>, typename MonolithicMappingProcess<SparseMatrixType>::Pointer, Process>(m, "MonolithicMappingProcess")
        .def(py::init<ModelPart&, ModelPart&, SparseMatrixType&>())
        ;

    typedef AcousticPMLDirectionProcess<SparseSpaceType, DenseSpaceType, LinearSolverType> AcousticPMLDirectionProcessType;
    py::class_<AcousticPMLDirectionProcessType, typename AcousticPMLDirectionProcessType::Pointer, Process>(m, "AcousticPMLDirectionProcess")
        .def(py::init<ModelPart&, ModelPart&, ModelPart&, typename LinearSolverType::Pointer>())
        ;

}

}  // namespace Python.

} // Namespace Kratos
