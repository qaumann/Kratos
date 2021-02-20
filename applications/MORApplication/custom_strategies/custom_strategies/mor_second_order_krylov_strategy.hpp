//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Quirin Aumann
//

#if !defined(MOR_SECOND_ORDER_KRYLOV_STRATEGY)
#define MOR_SECOND_ORDER_KRYLOV_STRATEGY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/builtin_timer.h"
#include "custom_strategies/custom_strategies/mor_offline_second_order_strategy.hpp"
#include "custom_utilities/orthogonalization_utility.hpp"
#include "custom_utilities/complex_sort_utility.hpp"

//default builder and solver
#include "custom_strategies/custom_builder_and_solvers/system_matrix_builder_and_solver.hpp"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}

///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class MorSecondOrderKrylovStrategy
 * @ingroup KratosCore
 * @brief This is the linear MOR matrix output strategy
 * @details This strategy builds the K and M matrices and outputs them
 * @author Aditya Ghantasala
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver, //= LinearSolver<TSparseSpace,TDenseSpace>
          class TReducedSparseSpace,
          class TReducedDenseSpace
          >
class MorSecondOrderKrylovStrategy
    // : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
    : public MorOfflineSecondOrderStrategy< TSparseSpace, TDenseSpace, TLinearSolver, TReducedSparseSpace, TReducedDenseSpace >
{
    using complex = std::complex<double>;

  public:
    ///@name Type Definitions
    ///@{
    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(MorSecondOrderKrylovStrategy);

    typedef TUblasSparseSpace<complex> ComplexSparseSpaceType;

    // typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef MorOfflineSecondOrderStrategy< TSparseSpace, TDenseSpace, TLinearSolver, TReducedSparseSpace, TReducedDenseSpace > BaseType;

    // typedef SystemMatrixBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver > TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef TDenseSpace DenseSpaceType;

    typedef typename TDenseSpace::MatrixType TDenseMatrixType;

    typedef typename TDenseSpace::MatrixPointerType TDenseMatrixPointerType;

    typedef typename BaseType::TSchemeType TSchemeType;

    //typedef typename BaseType::DofSetType DofSetType;

    typedef TLinearSolver TLinearSolverType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename TReducedSparseSpace::MatrixType TReducedSparseMatrixType;

    typedef typename TReducedSparseSpace::VectorType TReducedSparseVectorType;

    typedef typename TReducedDenseSpace::MatrixType TReducedDenseMatrixType;

    typedef typename ComplexSparseSpaceType::MatrixType ComplexSparseMatrixType;

    typedef typename ComplexSparseSpaceType::MatrixPointerType ComplexSparseMatrixPointerType;

    typedef typename ComplexSparseSpaceType::VectorType ComplexSparseVectorType;

    typedef typename ComplexSparseSpaceType::VectorPointerType ComplexSparseVectorPointerType;


    ///@}
    ///@name Life Cycle

    ///@{

    /**
     * Default constructor for the symmetric damped case and expansion points in complex conjugate pairs
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pBuilderAndSolver The builder and solver
     * @param pLinearSolver The complex linear solver to compute the basis
     * @param samplingPoints The complex sampling points in pairs of complex conjugate values
     */
    MorSecondOrderKrylovStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename BaseType::TBuilderAndSolverType::Pointer pBuilderAndSolver,
        typename TLinearSolverType::Pointer pNewLinearSolver,
        ComplexVector samplingPoints)
        : BaseType(rModelPart, pScheme, pBuilderAndSolver, pNewLinearSolver, true),
            mSamplingPoints(samplingPoints)
    {
        KRATOS_TRY;

        this->mSystemSizeR = mSamplingPoints.size();

        // check, if the sampling points are ordered in complex conjugate pairs
        ComplexSortUtility::PairComplexConjugates(mSamplingPoints);

        // use only one value of each pair to build a real-valued reduction basis
        mOrders.resize(this->mSystemSizeR);
        std::size_t n = 1;
        std::generate(mOrders.begin(), mOrders.end(), [&n] () { n++; return n%2; });


        KRATOS_CATCH("");
    }

    /**
     * Default constructor for the non-symmetric damped case and expansion points in complex conjugate pairs
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pBuilderAndSolver The builder and solver
     * @param pLinearSolver The complex linear solver to compute the basis
     * @param samplingPoints The complex sampling points in pairs of complex conjugate values
     */
    MorSecondOrderKrylovStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename BaseType::TBuilderAndSolverType::Pointer pBuilderAndSolver,
        typename TLinearSolverType::Pointer pNewLinearSolver,
        typename TLinearSolverType::Pointer pNewAdjointLinearSolver,
        ComplexVector samplingPoints)
        : BaseType(rModelPart, pScheme, pBuilderAndSolver, pNewLinearSolver, pNewAdjointLinearSolver, true),
            mSamplingPoints(samplingPoints)
    {
        KRATOS_TRY;

        this->mSystemSizeR = mSamplingPoints.size();

        // check, if the sampling points are ordered in complex conjugate pairs
        ComplexSortUtility::PairComplexConjugates(mSamplingPoints);

        // use only one value of each pair to build a real-valued reduction basis
        mOrders.resize(this->mSystemSizeR);
        std::size_t n = 0;
        std::generate(mOrders.begin(), mOrders.end(), [&n] () { n++; return n%2; });


        KRATOS_CATCH("");
    }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~MorSecondOrderKrylovStrategy() override
    {
        this->Clear();
    }


    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    virtual void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        if (this->mSolutionStepIsInitialized == false)
        {
            BaseType::InitializeSolutionStep();

            const std::size_t reduced_system_size = this->mSystemSizeR;
            const unsigned int system_size = this->GetBuilderAndSolver()->GetEquationSystemSize();

            TReducedDenseSpace::Resize(this->GetKr(), reduced_system_size, reduced_system_size);
            TReducedDenseSpace::Resize(this->GetDr(), reduced_system_size, reduced_system_size);
            TReducedDenseSpace::Resize(this->GetMr(), reduced_system_size, reduced_system_size);
            TReducedDenseSpace::Resize(this->GetBasis(), system_size, reduced_system_size);
            if( !this->SystemIsSymmetric() ) {
               TReducedDenseSpace::Resize(this->GetBasisLeft(), system_size, reduced_system_size);
            }


            this->mSolutionStepIsInitialized = true;
        }

        KRATOS_CATCH("");

    }


    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT: **/

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        KRATOS_TRY;
        // std::cout << "hello! this is where the second order krylov MOR magic happens" << std::endl;
        typename TSchemeType::Pointer p_scheme = this->GetScheme();
        typename BaseType::TBuilderAndSolverType::Pointer p_builder_and_solver = this->GetBuilderAndSolver();
        const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

        TSystemMatrixType& r_K = this->GetK();
        TSystemMatrixType& r_M = this->GetM();
        TSystemMatrixType& r_D = this->GetC();
        TSystemVectorType& r_RHS = this->GetSystemVector();
        TSystemVectorType& r_ov = this->GetOutputVector();
        TReducedDenseMatrixType& r_basis = this->GetBasis();
        TReducedDenseMatrixType& r_basis_l = this->GetBasisLeft();

        const std::size_t system_size = p_builder_and_solver->GetEquationSystemSize();
        const std::size_t n_sampling_points = mSamplingPoints.size();

        // copy mass matrix, damping matrix, and rhs to the complex space
        ComplexSparseMatrixType r_M_tmp(r_M.size1(), r_M.size2());
        ComplexMatrixUtility::ConstructMatrixStructure<TSchemeType, ComplexSparseMatrixType>(p_scheme,
            r_M_tmp,
            BaseType::GetModelPart(),
            p_builder_and_solver->GetEquationSystemSize());
        ComplexMatrixUtility::axpy<TSystemMatrixType, ComplexSparseMatrixType>(r_M, r_M_tmp, 1.0);

        ComplexSparseMatrixType r_D_tmp(r_D.size1(), r_D.size2());
        ComplexMatrixUtility::ConstructMatrixStructure<TSchemeType, ComplexSparseMatrixType>(p_scheme,
            r_D_tmp,
            BaseType::GetModelPart(),
            p_builder_and_solver->GetEquationSystemSize());
        ComplexMatrixUtility::axpy<TSystemMatrixType, ComplexSparseMatrixType>(r_D, r_D_tmp, 1.0);

        ComplexSparseVectorType r_RHS_tmp = ComplexSparseVectorType(r_RHS);
        ComplexSparseVectorType r_ov_tmp = ComplexSparseVectorType(r_ov);

        // create solution vector
        ComplexSparseVectorPointerType tmp_dx = ComplexSparseSpaceType::CreateEmptyVectorPointer();
        ComplexSparseVectorType& r_tmp_dx = *tmp_dx;
        ComplexSparseSpaceType::Resize(r_tmp_dx, system_size);
        ComplexSparseSpaceType::Set(r_tmp_dx, 0.0);

        // create dynamic stiffness matrix
        ComplexSparseMatrixType r_kdyn(system_size, system_size);
        ComplexMatrixUtility::ConstructMatrixStructure<TSchemeType, ComplexSparseMatrixType>(p_scheme,
            r_kdyn,
            BaseType::GetModelPart(),
            p_builder_and_solver->GetEquationSystemSize());

        // loop over sampling points
        std::size_t index = 0;
        std::size_t index_left = 0;
        BuiltinTimer basis_construction_time;

        for( std::size_t i=0; i<n_sampling_points; ++i ) {
            if( mOrders[i] == 0 ) {
                continue;
            }

            // build dynamic stiffness matrix
            ComplexSparseSpaceType::SetToZero(r_kdyn);
            ComplexMatrixUtility::axpy<TSystemMatrixType, ComplexSparseMatrixType>(r_K, r_kdyn, 1.0);
            ComplexMatrixUtility::axpy<ComplexSparseMatrixType, ComplexSparseMatrixType>(r_D_tmp, r_kdyn, mSamplingPoints(i));
            ComplexMatrixUtility::axpy<ComplexSparseMatrixType, ComplexSparseMatrixType>(r_M_tmp, r_kdyn, std::pow( mSamplingPoints(i), 2 ));

            // compute the basis vectors
            this->mpLinearSolver->Solve( r_kdyn, r_tmp_dx, r_RHS_tmp );
            UpdateBasis<typename TReducedSparseSpace::DataType>(r_basis, r_tmp_dx, index);

            if( !this->SystemIsSymmetric() ) {
                this->mpAdjointLinearSolver->Solve(r_kdyn, r_tmp_dx, r_ov_tmp);
                UpdateBasis<typename TReducedSparseSpace::DataType>(r_basis_l, r_tmp_dx, index_left);
            }
        }

        KRATOS_INFO_IF("Basis Construction Time", BaseType::GetEchoLevel() > 0 && rank == 0)
            << basis_construction_time.ElapsedSeconds() << std::endl;

        // orthogonalization step
        BuiltinTimer basis_orthogonalization_time;

        OrthogonalizationUtility::OrthogonalizeQR<TReducedDenseSpace>(r_basis);
        if( !this->SystemIsSymmetric() ) {
            OrthogonalizationUtility::OrthogonalizeQR<TReducedDenseSpace>(r_basis_l);
        }

        KRATOS_INFO_IF("Basis Orthogonalization Time", BaseType::GetEchoLevel() > 0 && rank == 0)
            << basis_orthogonalization_time.ElapsedSeconds() << std::endl;

		return true;

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "MorSecondOrderKrylovStrategy";
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

  protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

  private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    bool mUseComplexFlag;
    ComplexVector mSamplingPoints;
    std::vector<std::size_t> mOrders;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * Basis update for a complex basis
     * @param r_dx the complex vector to be included in the basis
     * @param index the column index
     */
    template<typename TScalar, typename std::enable_if<std::is_same<std::complex<double>, TScalar>::value, int>::type = 0>
    void UpdateBasis(TReducedDenseMatrixType& rBasis, ComplexSparseVectorType& rDx, std::size_t& index)
    {
        // TReducedDenseMatrixType& r_basis = this->GetBasis();
        column( rBasis, index ) = rDx;
        index++;
    }

    /**
     * Basis update for a real basis. Real and imaginary part will be placen in two consecutive columns
     * @param r_dx the complex vector to be included in the basis
     * @param index the column index
     */
    template<typename TScalar, typename std::enable_if<std::is_same<double, TScalar>::value, int>::type = 0>
    void UpdateBasis(TReducedDenseMatrixType& rBasis, ComplexSparseVectorType& rDx, std::size_t& index)
    {
        // TReducedDenseMatrixType& r_basis = this->GetBasis();
        column( rBasis, index ) = real(rDx);
        column( rBasis, index+1 ) = imag(rDx);
        index = index + 2;
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /**
     * Copy constructor.
     */

    MorSecondOrderKrylovStrategy(const MorSecondOrderKrylovStrategy &Other){};

    ///@}

}; /* Class MorSecondOrderKrylovStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* MOR_SECOND_ORDER_KRYLOV_STRATEGY  defined */