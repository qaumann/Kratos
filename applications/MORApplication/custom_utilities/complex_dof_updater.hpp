#if !defined(KRATOS_COMPLEX_DOF_UPDATER)
#define KRATOS_COMPLEX_DOF_UPDATER

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
namespace ComplexDofUpdater
{

    void AssignComplexDofValue(
        Variable<double>& var,
        Variable<double>& var_real,
        Variable<double>& var_imag,
        Node<3>& node,
        const std::complex<double>& Z,
        bool signed_abs)
    {
        if( signed_abs ) {
            const int sign = std::real(Z) >= 0 ? 1 : -1;
            node.FastGetSolutionStepValue(var) = sign * std::abs(Z);
        } else {
            node.FastGetSolutionStepValue(var) = std::abs(Z);
        }
        if( node.SolutionStepsDataHas(var_real) && node.SolutionStepsDataHas(var_imag) ) {
            node.FastGetSolutionStepValue(var_real) = std::real(Z);
            node.FastGetSolutionStepValue(var_imag) = std::imag(Z);
        } else {
            KRATOS_WARNING("ComplexDofUpdater") << "Complex variable for " << var.Name() <<
                " is not available. Using absolute value only."  << std::endl;
        }
    }

#ifdef KRATOS_MOR_ASSIGN_COMPLEX_DOF
#undef KRATOS_MOR_ASSIGN_COMPLEX_DOF
#endif
#define KRATOS_MOR_ASSIGN_COMPLEX_DOF(name, it_node, signed_abs) \
    if( it_node->HasDofFor(name) ) \
    { \
        const std::size_t eq_id = it_node->GetDof(name).EquationId(); \
        AssignComplexDofValue(name, REAL_##name, IMAG_##name, (*it_node), rX(eq_id), signed_abs); \
    }

    /** @brief Assign new values for the problem's degrees of freedom using the complex vector rX.
     *  @details value = std::abs(rX[dof.EquationId()])
     *      REAL_value = std::real(rX[dof.EquationId()])
     *      IMAG_value = std::imag(rX[dof.EquationId()])
     *  @param[in/out] rModelPart The model part.
     *  @param[in] rX The solution vector.
     */
    template<typename ComplexSystemVectorType>
    void AssignDofs(
        ModelPart& rModelPart,
        const ComplexSystemVectorType& rX)
    {
        const int num_nodes = static_cast<int>( rModelPart.Nodes().size() );

        #pragma omp parallel for
        for( int i=0; i<num_nodes; ++i )
        {
            auto it_node = std::begin(rModelPart.Nodes()) + i;

            KRATOS_MOR_ASSIGN_COMPLEX_DOF(DISPLACEMENT_X, it_node, true);
            KRATOS_MOR_ASSIGN_COMPLEX_DOF(DISPLACEMENT_Y, it_node, true);
            KRATOS_MOR_ASSIGN_COMPLEX_DOF(DISPLACEMENT_Z, it_node, true);
            KRATOS_MOR_ASSIGN_COMPLEX_DOF(PRESSURE, it_node, false);
        }
    }

} // namespace ComplexDofUpdater

} // namespace Kratos

#endif