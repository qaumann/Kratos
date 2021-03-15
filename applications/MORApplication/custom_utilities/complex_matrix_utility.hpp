#if !defined(KATOS_COMPLEX_MATRIX_UTILITY_H_INCLUDED)
#define KATOS_COMPLEX_MATRIX_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
namespace ComplexMatrixUtility
{

    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    typedef PointerVectorSet<Element, IndexedObject> ElementsContainerType;


    /**
     * @brief Constructs the matrix structure of a complex sparse matrix
     * @details ..
     * @param pScheme the problem's dof set
     * @param A array where 1 specifies a fixed dof´
     * @param rModelPart
     * @param equationSize
     */
    template<typename TSchemeType, typename TComplexMatrixType>
    void ConstructMatrixStructure(
        typename TSchemeType::Pointer pScheme,
        TComplexMatrixType& A,
        ModelPart& rModelPart,
        std::size_t equationSize)
    {
        //filling with zero the matrix (creating the structure)
        Timer::Start("MatrixStructure");

        // Getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        // Getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

        std::vector< LockObject > lock_array(equationSize);

        std::vector<std::unordered_set<std::size_t> > indices(equationSize);
        #pragma omp parallel for firstprivate(equationSize)
        for (int iii = 0; iii < static_cast<int>(equationSize); iii++) {
            indices[iii].reserve(40);
        }

        Element::EquationIdVectorType ids(3, 0);

        #pragma omp parallel for firstprivate(nelements, ids)
        for (int iii=0; iii<nelements; iii++) {
            typename ElementsContainerType::iterator i_element = el_begin + iii;
            pScheme->EquationId(*i_element, ids, CurrentProcessInfo);
            for (std::size_t i = 0; i < ids.size(); i++) {
                lock_array[ids[i]].SetLock();
                auto& row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());
                lock_array[ids[i]].UnSetLock();
            }
        }

        #pragma omp parallel for firstprivate(nconditions, ids)
        for (int iii = 0; iii<nconditions; iii++) {
            typename ConditionsArrayType::iterator i_condition = cond_begin + iii;
            pScheme->EquationId(*i_condition, ids, CurrentProcessInfo);
            for (std::size_t i = 0; i < ids.size(); i++) {
                lock_array[ids[i]].SetLock();
                auto& row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());
                lock_array[ids[i]].UnSetLock();
            }
        }

        //destroy locks
        lock_array = std::vector< LockObject >();

        //count the row sizes
        unsigned int nnz = 0;
        for (unsigned int i = 0; i < indices.size(); i++) {
            nnz += indices[i].size();
        }

        A = TComplexMatrixType(indices.size(), indices.size(), nnz);

        std::complex<double>* Avalues = A.value_data().begin();
        std::size_t* Arow_indices = A.index1_data().begin();
        std::size_t* Acol_indices = A.index2_data().begin();

        //filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
        Arow_indices[0] = 0;
        for (int i = 0; i < static_cast<int>(A.size1()); i++) {
            Arow_indices[i+1] = Arow_indices[i] + indices[i].size();
        }

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(A.size1()); i++) {
            const unsigned int row_begin = Arow_indices[i];
            const unsigned int row_end = Arow_indices[i+1];
            unsigned int k = row_begin;
            for (auto it = indices[i].begin(); it != indices[i].end(); it++) {
                Acol_indices[k] = *it;
                Avalues[k] = std::complex<double>(0.0, 0.0);
                k++;
            }

            indices[i].clear(); //deallocating the memory

            std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);

        }

        A.set_filled(indices.size()+1, nnz);

        Timer::Stop("MatrixStructure");
    }

    /**
     * @brief Compute y += a*x (axpy) for different matrix types
     * @details A must be double
     */
    template<typename MatrixType, typename OtherMatrixType>
    void axpy(MatrixType& rX, OtherMatrixType& rY, const typename MatrixType::value_type A)
    {
        const std::size_t size = rX.value_data().size();
        KRATOS_ERROR_IF_NOT(size == rY.value_data().size()) << "wrong size!" << std::endl;

        if( A == 1.0 ) {
            #pragma omp parallel for// schedule(static)
            for( int i=0; i<static_cast<int>(size); ++i ) {
                rY.value_data()[i] += rX.value_data()[i];
            }
        } else if(A == -1.0 ) {
            #pragma omp parallel for// schedule(static)
            for( int i=0; i<static_cast<int>(size); ++i ) {
                rY.value_data()[i] -= rX.value_data()[i];
            }
        } else {
            #pragma omp parallel for// schedule(static)
            for( int i=0; i<static_cast<int>(size); ++i ) {
                rY.value_data()[i] += A * rX.value_data()[i];
            }
        }
    }

} // namespace ComplexMatrixUtility

} // namespace Kratos

#endif