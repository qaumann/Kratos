#if !defined(KRATOS_DIRICHLET_UTILITY_H_INCLUDED)
#define KRATOS_DIRICHLET_UTILITY_H_INCLUDED

// System includes

// External includes
// #include <Eigen/QR>

// Project includes
#include "includes/define.h"

namespace Kratos
{
namespace DirichletUtility
{

    /**
     * @brief Finds all fixed dofs
     * @details vector<bool> is not thread safe, so we use unsigned int here
     * @param rDofSet the problem's dof set
     * @param rFixedDofs array where 1 specifies a fixed dof´
     */
    template<typename DofsArrayType>
    void GetDirichletConstraints(DofsArrayType& rDofSet, std::vector<unsigned int>& rFixedDofs)
    {
        const size_t n_dofs = rDofSet.size();

        if( rFixedDofs.size() != n_dofs )
            rFixedDofs.resize( n_dofs, false);

        //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        #pragma omp parallel for firstprivate(n_dofs)
        for( int k = 0; k < static_cast<int>(n_dofs); k++ )
        {
            auto dof_iterator = std::begin(rDofSet) + k;
            rFixedDofs[k] = dof_iterator->IsFixed() ? 1 : 0;
        }
    }

    /**
     * @brief Applies the dirichlet conditions. This operation may be very heavy or completely
     * unexpensive depending on the implementation choosen and on how the System Matrix is built.
     * This should be part of the builder and solver.
     * @details For explanation of how it works for a particular implementation the user
     * should refer to the particular Builder And Solver choosen
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     * @param b The RHS vector
     */
    template<typename SparseSpaceType>
    void ApplyDirichletConditions(
        typename SparseSpaceType::MatrixType& A,
        typename SparseSpaceType::VectorType& b,
        const std::vector<unsigned int>& FixedDofSet,
        const typename SparseSpaceType::DataType Factor)
    {
        std::size_t system_size = A.size1();
        typename SparseSpaceType::DataType* Avalues = A.value_data().begin();
        std::size_t* Arow_indices = A.index1_data().begin();
        std::size_t* Acol_indices = A.index2_data().begin();

        //detect if there is a line of all zeros and set the diagonal to a 1 if this happens
        if( std::abs(Factor) != 0.0 )
        {
            #pragma omp parallel for firstprivate(system_size)
            for (int k = 0; k < static_cast<int>(system_size); ++k){
                std::size_t col_begin = Arow_indices[k];
                std::size_t col_end = Arow_indices[k+1];
                bool empty = true;
                for (std::size_t j = col_begin; j < col_end; ++j)
                {
                    if(Avalues[j] != 0.0)
                    {
                        empty = false;
                        break;
                    }
                }

                if(empty == true)
                {
                    A(k,k) = 1.0;
                    b[k] = 0.0;
                }
            }
        }

        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(system_size); ++k)
        {
            std::size_t col_begin = Arow_indices[k];
            std::size_t col_end = Arow_indices[k+1];
            const unsigned int is_fixed = FixedDofSet[k];
            if( is_fixed == 1 )
            {
                // zero out the whole row, except the diagonal
                for (std::size_t j = col_begin; j < col_end; ++j)
                    if (static_cast<int>(Acol_indices[j]) != k )
                        Avalues[j] = 0.0;
                    else
                        Avalues[j] = Factor;

                // zero out the RHS
                b[k] = 0.0;
            }
            else
            {
                // zero out the column which is associated with the zero'ed row
                for( std::size_t j = col_begin; j < col_end; ++j )
                    if( FixedDofSet[ Acol_indices[j] ] )
                        Avalues[j] = 0.0;
            }
        }
    }

    template<typename ComplexSparseSpaceType>
    void SetComplexDirichletConditions(
        ModelPart& rModelPart,
        typename ComplexSparseSpaceType::MatrixType& rA,
        typename ComplexSparseSpaceType::VectorType& rb,
        const std::vector<unsigned int>& FixedDofSet,
        const Variable<double>& rVariable,
        const Variable<double>& rValueVariableReal,
        const Variable<double>& rValueVariableImag
    )
    {
        typename ComplexSparseSpaceType::VectorType x(rb.size());
        ComplexSparseSpaceType::SetToZero(x);
        const int num_nodes = static_cast<int>( rModelPart.Nodes().size() );
        bool empty = true;

        // add contribution to rhs
        #pragma omp parallel for
        for( int i=0; i<num_nodes; ++i ) {
            auto it_node = std::begin(rModelPart.Nodes()) + i;
            if( it_node->HasDofFor(rVariable) ) {
                if( it_node->IsFixed(rVariable) ) {
                    const size_t pos = it_node->GetDof(rVariable).GetId() - 1;
                    std::complex<double> val(it_node->GetValue(rValueVariableReal), it_node->GetValue(rValueVariableImag));
                    x[pos] = val;
                    if( std::abs(val) > std::numeric_limits<double>::epsilon() ) {
                        empty = false;
                    }
                }
            }
        }

        if( !empty ) {
            axpy_prod( rA, -x, rb, false);
        }

        DirichletUtility::ApplyDirichletConditions<ComplexSparseSpaceType>(rA, rb, FixedDofSet, 1.0);

        if( !empty ) {
            #pragma omp parallel for
            for( int i=0; i<static_cast<int>(rb.size()); ++i ) {
                if( std::abs(x[i]) > std::numeric_limits<double>::epsilon() ) {
                    rb[i] = x[i];
                }
            }
        }
    }


} // namespace DirichletUtility

} // namespace Kratos

#endif