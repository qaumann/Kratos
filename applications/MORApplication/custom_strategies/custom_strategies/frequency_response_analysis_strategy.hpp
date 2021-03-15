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

#if !defined(FREQUENCY_RESPONSE_ANALYSIS_STRATEGY)
#define FREQUENCY_RESPONSE_ANALYSIS_STRATEGY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "utilities/builtin_timer.h"

#include "custom_utilities/dirichlet_utility.hpp"
#include "custom_utilities/complex_dof_updater.hpp"
#include "custom_utilities/complex_matrix_utility.hpp"

namespace Kratos
{

namespace { // settings namespace
    template<class SpaceType>
    struct FrequencyDependentSettings
    {
        FrequencyDependentSettings(ModelPart& rModelPart, int buildLevel, bool matrixContribution, bool vectorContribution) :
            model_part(&rModelPart),
            build_level(buildLevel),
            has_matrix_contribution(matrixContribution),
            has_vector_contribution(vectorContribution) {};

        ModelPart* model_part;
        int build_level;
        std::complex<double> factor = std::complex<double>(0,0);
        typename SpaceType::MatrixPointerType matrix = SpaceType::CreateEmptyMatrixPointer();
        typename SpaceType::VectorPointerType vector = SpaceType::CreateEmptyVectorPointer();

        bool initialized = false;
        bool has_matrix_contribution;
        bool has_vector_contribution;
    };
}

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
 * @class FrequencyResponseAnalysisStrategy
 * @ingroup MORApplication
 * @brief This is the harmonic frequency response strategy
 * @details This strategy performs a harmonic analysis and outputs the absolute values of the solution to the dofs.
 *      Rayleigh damping is enabled by default. Modal damping (i.e. complex stiffness damping) can be enabled using
 *      the UseModalDamping flag. The property RAYLEIGH_BETA is then treated as the modal damping factor. RAYLEIGH_ALPHA
 *      has to be set to zero if modal damping is to be used!
 * @author Quirin Aumann
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class FrequencyResponseAnalysisStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
    using complex = std::complex<double>;

  public:
    ///@name Type Definitions
    ///@
    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(FrequencyResponseAnalysisStrategy);

    typedef TUblasSparseSpace<complex> ComplexSparseSpaceType;
    typedef TUblasDenseSpace<complex> ComplexDenseSpaceType;

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver > TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef TDenseSpace DenseSpaceType;

    typedef typename TDenseSpace::MatrixType TDenseMatrixType;

    typedef typename TDenseSpace::MatrixPointerType TDenseMatrixPointerType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename ComplexSparseSpaceType::MatrixType TSolutionMatrixType;

    typedef typename ComplexSparseSpaceType::MatrixPointerType TSolutionMatrixPointerType;

    typedef typename ComplexSparseSpaceType::VectorType TSolutionVectorType;

    typedef typename ComplexSparseSpaceType::VectorPointerType TSolutionVectorPointerType;

    typedef ComplexSparseSpaceType TSolutionSpace;

    typedef std::complex<double> ComplexType;

    typedef LinearSolver<ComplexSparseSpaceType, ComplexDenseSpaceType> ComplexLinearSolverType;

    ///@}
    ///@name Life Cycle

    ///@{

    /**
     * Default constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration schemed
     * @param MoveMeshFlag The flag that allows to move the mesh
     * @param UseModalDamping True if modal damping is used
     */
    FrequencyResponseAnalysisStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver,
        typename ComplexLinearSolverType::Pointer pNewcomplexLinearSolver,
        bool MoveMeshFlag = false,
        bool UseModalDamping = false)
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, MoveMeshFlag),
            mpScheme(pScheme),
            mpBuilderAndSolver(pBuilderAndSolver),
            mpComplexLinearSolver(pNewcomplexLinearSolver),
            mUseModalDamping(UseModalDamping)
    {
        KRATOS_TRY;

        // Set flags to start correcty the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;
        mReformDofSetAtEachStep = false;

        // Tells to the builder and solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(false);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        // be reshaped at each step or not
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are never rebuilt
        this->SetRebuildLevel(0);

        mpA = ComplexSparseSpaceType::CreateEmptyMatrixPointer();
        mpK = TSparseSpace::CreateEmptyMatrixPointer();
        mpKi = TSparseSpace::CreateEmptyMatrixPointer();
        mpC = TSparseSpace::CreateEmptyMatrixPointer();
        mpM = TSparseSpace::CreateEmptyMatrixPointer();
        mpMi = TSparseSpace::CreateEmptyMatrixPointer();
        mpRHS = ComplexSparseSpaceType::CreateEmptyVectorPointer();
        mpDx = ComplexSparseSpaceType::CreateEmptyVectorPointer();
        mpB = TSparseSpace::CreateEmptyVectorPointer();

        KRATOS_CATCH("");
    }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~FrequencyResponseAnalysisStrategy() override
    {
        Clear();
    }

    /**
     * @brief Set method for the time scheme
     * @param pScheme The pointer to the time scheme considered
     */
    void SetScheme(typename TSchemeType::Pointer pScheme)
    {
        mpScheme = pScheme;
    };

    /**
     * @brief Get method for the time scheme
     * @return mpScheme: The pointer to the time scheme considered
     */
    typename TSchemeType::Pointer GetScheme()
    {
        return mpScheme;
    };

    /**
     * @brief Set method for the builder and solver
     * @param pNewBuilderAndSolver The pointer to the builder and solver considered
     */
    void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver)
    {
        mpBuilderAndSolver = pNewBuilderAndSolver;
    };

    /**
     * @brief Get method for the builder and solver
     * @return mpBuilderAndSolver: The pointer to the builder and solver considered
     */
    typename TBuilderAndSolverType::Pointer GetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

    /**
     * @brief This method sets the flag mInitializeWasPerformed
     * @param InitializePerformedFlag The flag that tells if the initialize has been computed
     */
    void SetInitializePerformedFlag(bool InitializePerformedFlag = true)
    {
        mInitializeWasPerformed = InitializePerformedFlag;
    }

    /**
     * @brief This method gets the flag mInitializeWasPerformed
     * @return mInitializeWasPerformed: The flag that tells if the initialize has been computed
     */
    bool GetInitializePerformedFlag()
    {
        return mInitializeWasPerformed;
    }

    /**
     * @brief This method sets the flag mReformDofSetAtEachStep
     * @param Flag The flag that tells if each time step the system is rebuilt
     */
    void SetReformDofSetAtEachStepFlag(bool Flag)
    {
        mReformDofSetAtEachStep = Flag;
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
    }

    /**
     * @brief This method returns the flag mReformDofSetAtEachStep
     * @return The flag that tells if each time step the system is rebuilt
     */
    bool GetReformDofSetAtEachStepFlag()
    {
        return mReformDofSetAtEachStep;
    }

    /**
     * @brief It sets the level of echo for the solving strategy
     * @param Level The level to set
     * @details The different levels of echo are:
     * - 0: Mute... no echo at all
     * - 1: Printing time and basic informations
     * - 2: Printing linear solver data
     * - 3: Print of debug informations: Echo of stiffness matrix, Dx, b...
     */

    void SetEchoLevel(int Level) override
    {
        BaseType::mEchoLevel = Level;
        GetBuilderAndSolver()->SetEchoLevel(Level);
    }

    void SetFrequencyDependentMaterial(FrequencyDependentSettings<TSparseSpace>* settings)
    {
        mFrequencyDependentSettings.push_back(settings);
        mUseFrequencyDependentMaterials = true;
    }

    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize() override
    {
        KRATOS_TRY;

        if (mInitializeWasPerformed == false)
        {

            //pointers needed in the solution
            typename TSchemeType::Pointer p_scheme = GetScheme();

            //Initialize The Scheme - OPERATIONS TO BE DONE ONCE
            if (p_scheme->SchemeIsInitialized() == false)
                p_scheme->Initialize(BaseType::GetModelPart());

            //Initialize The Elements - OPERATIONS TO BE DONE ONCE
            if (p_scheme->ElementsAreInitialized() == false)
                p_scheme->InitializeElements(BaseType::GetModelPart());

            //Initialize The Conditions - OPERATIONS TO BE DONE ONCE
            if (p_scheme->ConditionsAreInitialized() == false)
                p_scheme->InitializeConditions(BaseType::GetModelPart());

            mInitializeWasPerformed = true;
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Clears the internal storage
     */
    void Clear() override
    {
        KRATOS_TRY;

        // if the preconditioner is saved between solves, it
        // should be cleared here.
        GetBuilderAndSolver()->GetLinearSystemSolver()->Clear();

        if (mpA != nullptr)
            TSolutionSpace::Clear(mpA);
        if (mpK != nullptr)
            SparseSpaceType::Clear(mpK);
        if (mpKi != nullptr)
            SparseSpaceType::Clear(mpKi);
        if (mpM != nullptr)
            SparseSpaceType::Clear(mpM);
        if (mpMi != nullptr)
            SparseSpaceType::Clear(mpMi);
        if (mpC != nullptr)
            SparseSpaceType::Clear(mpC);
        if (mpRHS != nullptr)
            TSolutionSpace::Clear(mpRHS);
        if (mpDx != nullptr)
            TSolutionSpace::Clear(mpDx);
        if (mpB != nullptr)
            SparseSpaceType::Clear(mpB);

        //setting to zero the internal flag to ensure that the dof sets are recalculated
        GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
        GetBuilderAndSolver()->Clear();
        GetScheme()->Clear();

        mInitializeWasPerformed = false;

        KRATOS_CATCH("");
    }


    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        if (mSolutionStepIsInitialized == false)
        {
            typename TSchemeType::Pointer p_scheme = GetScheme();
            typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();
            ModelPart& r_model_part = BaseType::GetModelPart();

            const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();
            //set up the system, operation performed just once unless it is required
            //to reform the dof set at each iteration
            TSystemMatrixType& r_K  = *mpK;
            TSystemMatrixType& r_Ki  = *mpKi;
            TSystemMatrixType& r_M = *mpM;
            TSystemMatrixType& r_Mi = *mpMi;
            TSystemMatrixType& r_C  = *mpC;
            TSolutionMatrixType& r_A  = *mpA;
            // TSolutionVectorType& r_RHS  = *mpRHS;
            TSystemVectorType& r_B = *mpB;

            if (p_builder_and_solver->GetDofSetIsInitializedFlag() == false ||
                mReformDofSetAtEachStep == true)
            {
                BuiltinTimer system_construction_time;
                //setting up the list of the DOFs to be solved
                BuiltinTimer setup_dofs_time;
                p_builder_and_solver->SetUpDofSet(p_scheme, BaseType::GetModelPart());
                KRATOS_INFO_IF("Setup Dofs Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << setup_dofs_time.ElapsedSeconds() << std::endl;

                //shaping correctly the system
                BuiltinTimer setup_system_time;
                p_builder_and_solver->SetUpSystem(BaseType::GetModelPart());
                KRATOS_INFO_IF("Setup System Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << setup_system_time.ElapsedSeconds() << std::endl;

                //setting up the Vectors involved to the correct size
                // TSystemVectorPointerType tmp_RHS = TSparseSpace::CreateEmptyVectorPointer();
                // TSystemVectorType& r_tmp_RHS = *tmp_RHS;

                BuiltinTimer system_matrix_resize_time;
                p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpK, mpB, mpB,
                                                                 BaseType::GetModelPart());
                p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpM, mpB, mpB,
                                                                 BaseType::GetModelPart());
                p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpC, mpB, mpB,
                                                                 BaseType::GetModelPart());
                if (mpRHS->size() != p_builder_and_solver->GetEquationSystemSize())
                    mpRHS->resize(p_builder_and_solver->GetEquationSystemSize(), false);
                if (mpDx->size() != p_builder_and_solver->GetEquationSystemSize())
                    mpDx->resize(p_builder_and_solver->GetEquationSystemSize(), false);

                KRATOS_INFO_IF("System Matrix Resize Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << system_matrix_resize_time.ElapsedSeconds() << std::endl;

                //set up system matrices
                std::vector<unsigned int> fixed_dofs;
                DirichletUtility::GetDirichletConstraints<typename TBuilderAndSolverType::DofsArrayType>(p_builder_and_solver->GetDofSet(), fixed_dofs);

                TSystemVectorPointerType p_tmp = TSparseSpace::CreateEmptyVectorPointer();
                TSystemVectorType tmp(r_K.size1(), 0.0);

                //set up the stiffness matrix and rhs
                r_model_part.GetProcessInfo()[BUILD_LEVEL] = 1;
                p_scheme->InitializeNonLinIteration(BaseType::GetModelPart(), r_K, tmp, tmp);
                p_builder_and_solver->Build(p_scheme, BaseType::GetModelPart(), r_K, r_B);

                //if required, set up the imaginary part of the stiffness and mass
                if( mUseModalDamping )
                {
                    p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpKi, p_tmp, p_tmp,
                                                                    BaseType::GetModelPart());
                    r_model_part.GetProcessInfo()[BUILD_LEVEL] = 111;
                    p_builder_and_solver->Build(p_scheme, BaseType::GetModelPart(), r_Ki, tmp);

                    p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpMi, p_tmp, p_tmp,
                                                                    BaseType::GetModelPart());
                    r_model_part.GetProcessInfo()[BUILD_LEVEL] = 121;
                    p_builder_and_solver->Build(p_scheme, BaseType::GetModelPart(), r_Mi, tmp);
                }

                //set up the damping matrix
                r_model_part.GetProcessInfo()[BUILD_LEVEL] = 101;
                p_builder_and_solver->Build(p_scheme, BaseType::GetModelPart(), r_C, tmp);

                //set up the mass matrix
                r_model_part.GetProcessInfo()[BUILD_LEVEL] = 201;
                p_builder_and_solver->Build(p_scheme, BaseType::GetModelPart(), r_M, tmp);

                //if required, set up porous Biot material matrices
                if( mUseFrequencyDependentMaterials ) {
                    for( auto& setting : mFrequencyDependentSettings ) {
                        InitializeFrequencyDependentMatrices(*setting, fixed_dofs);
                    }
                }

                p_scheme->FinalizeNonLinIteration(BaseType::GetModelPart(), r_K, tmp, tmp);

                //resize and itialize working matrix
                ComplexMatrixUtility::ConstructMatrixStructure<TSchemeType, TSolutionMatrixType>(p_scheme,
                    r_A,
                    BaseType::GetModelPart(),
                    p_builder_and_solver->GetEquationSystemSize());

            KRATOS_INFO_IF("System Construction Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << system_construction_time.ElapsedSeconds() << std::endl;
            }


            //initial operations ... things that are constant over the Solution Step
            //nothing to do here

            mSolutionStepIsInitialized = true;
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        typename TSchemeType::Pointer p_scheme = GetScheme();

        //Cleaning memory after the solution
        p_scheme->Clean();

        //reset flags for next step
        mSolutionStepIsInitialized = false;

        if (mReformDofSetAtEachStep == true) //deallocate the systemvectors
        {
            TSolutionSpace::Clear(mpA);
            SparseSpaceType::Clear(mpK);
            SparseSpaceType::Clear(mpKi);
            SparseSpaceType::Clear(mpM);
            SparseSpaceType::Clear(mpMi);
            SparseSpaceType::Clear(mpC);
            TSolutionSpace::Clear(mpRHS);
            TSolutionSpace::Clear(mpDx);
            SparseSpaceType::Clear(mpB);

            if( mUseFrequencyDependentMaterials ) {
                for( auto& setting : mFrequencyDependentSettings ) {
                    SparseSpaceType::Clear(setting->matrix);
                    SparseSpaceType::Clear(setting->vector);
                }
            }

            this->Clear();
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        KRATOS_TRY;

        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();
        const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

        TSolutionMatrixType& r_A = *mpA;
        TSystemMatrixType& r_K = *mpK;
        TSystemMatrixType& r_Ki = *mpKi;
        TSystemMatrixType& r_M = *mpM;
        TSystemMatrixType& r_Mi = *mpMi;
        TSystemMatrixType& r_C  = *mpC;
        TSolutionVectorType& r_RHS = *mpRHS;
        TSolutionVectorType& r_Dx = *mpDx;
        TSystemVectorType& r_B = *mpB;

        auto& r_process_info = BaseType::GetModelPart().GetProcessInfo();
        const double excitation_frequency = r_process_info[FREQUENCY];
        KRATOS_ERROR_IF(excitation_frequency <= 0) << "Invalid excitation frequency" << std::endl;
        const double excitation_frequency2 = std::pow(excitation_frequency, 2);

        TSolutionSpace::SetToZero(r_A);
        TSolutionSpace::SetToZero(r_RHS);

        BuiltinTimer build_time;
        ComplexMatrixUtility::axpy<TSystemMatrixType, TSolutionMatrixType>(r_C, r_A, excitation_frequency);

        if( mUseModalDamping ) {
            ComplexMatrixUtility::axpy<TSystemMatrixType, TSolutionMatrixType>(r_Ki, r_A, 1.0);
            ComplexMatrixUtility::axpy<TSystemMatrixType, TSolutionMatrixType>(r_Mi, r_A, -excitation_frequency2);
        }

        if( mUseFrequencyDependentMaterials ) {
            double factor;
            for( auto& setting : mFrequencyDependentSettings ) {
                // KRATOS_WATCH(setting->factor)
                factor = std::imag(setting->factor);
                if( std::abs(factor) < std::numeric_limits<double>::epsilon() ) {
                    continue;
                }
                if( setting->has_matrix_contribution ) {
                    TSystemMatrixType& r_Mat = *(setting->matrix);
                    ComplexMatrixUtility::axpy<TSystemMatrixType, TSolutionMatrixType>(r_Mat, r_A, factor);
                }
                if( setting->has_vector_contribution ) {
                    TSystemVectorType& r_Vec = *(setting->vector);
                    r_RHS += factor * r_Vec;
                }
            }
        }

        r_A *= complex(0,1);
        r_RHS *= complex(0,1);

        ComplexMatrixUtility::axpy<TSystemMatrixType, TSolutionMatrixType>(r_K, r_A, 1.0);
        ComplexMatrixUtility::axpy<TSystemMatrixType, TSolutionMatrixType>(r_M, r_A, -excitation_frequency2);
        r_RHS += r_B;

        if( mUseFrequencyDependentMaterials ) {
            double factor;
            for( auto& setting : mFrequencyDependentSettings ) {
                // KRATOS_WATCH(setting->factor)
                factor = std::real(setting->factor);
                if( std::abs(factor) < std::numeric_limits<double>::epsilon() ) {
                    continue;
                }
                if( setting->has_matrix_contribution ) {
                    TSystemMatrixType& r_Mat = *(setting->matrix);
                    ComplexMatrixUtility::axpy<TSystemMatrixType, TSolutionMatrixType>(r_Mat, r_A, factor);
                }
                if( setting->has_vector_contribution ) {
                    TSystemVectorType& r_Vec = *(setting->vector);
                    r_RHS += factor * r_Vec;
                }
            }
        }

        // add complex contribution to rhs for complex dirichlet condition
        ApplyComplexDirichletConditions();

        KRATOS_INFO_IF("Dynamic Stiffness Matrix Build Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << build_time.ElapsedSeconds() << std::endl;

        // if required, write dynamic stiffness and rhs
        if( BaseType::GetEchoLevel() > 4 ) {
            MatrixOutput<ComplexSparseSpaceType>(r_A, "A_" + std::to_string(excitation_frequency));
            VectorOutput<ComplexSparseSpaceType>(r_RHS, "RHS_" + std::to_string(excitation_frequency));
        }

        //Solve the system
        BuiltinTimer solve_time;
        mpComplexLinearSolver->Solve(r_A, r_Dx, r_RHS);
        KRATOS_INFO_IF("Linear Solve Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << solve_time.ElapsedSeconds() << std::endl;

        //Assign the computed values
        BuiltinTimer assign_dofs_time;
        ComplexDofUpdater::AssignDofs<TSolutionVectorType>(BaseType::GetModelPart(), r_Dx);
        KRATOS_INFO_IF("Assign Dofs Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << assign_dofs_time.ElapsedSeconds() << std::endl;
		return true;

        KRATOS_CATCH("");
    }

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        BaseType::Check();

        GetBuilderAndSolver()->Check(BaseType::GetModelPart());

        GetScheme()->Check(BaseType::GetModelPart());

        return 0;

        KRATOS_CATCH("")
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

    /**
     * @brief This method returns the LHS matrix
     * @return The LHS matrix
     */
    TSolutionMatrixType &GetAssembledSystemMatrix()
    {
        TSolutionMatrixType &mA = *mpA;

        return mA;
    }

    TSystemMatrixType &GetStiffnessMatrix()
    {
        TSystemMatrixType &mK = *mpK;

        return mK;
    }

    TSystemMatrixType &GetImaginaryStiffnessMatrix()
    {
        TSystemMatrixType &mKi = *mpKi;

        return mKi;
    }

    TSystemMatrixType &GetMassMatrix()
    {
        TSystemMatrixType &mM = *mpM;

        return mM;
    }

    TSystemMatrixType &GetDampingMatrix()
    {
        TSystemMatrixType &mS = *mpC;

        return mS;
    }

    /**
     * @brief This method returns the RHS vector
     * @return The RHS vector
     */
    TSolutionVectorType& GetAssembledSystemVector()
    {
        TSolutionVectorType& mb = *mpRHS;

        return mb;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief This method returns the components of the system of equations depending of the echo level
     */
    virtual void EchoInfo()
    {
        TSystemMatrixType& rK  = *mpK;
        TSystemMatrixType& rKi = *mpKi;
        TSystemMatrixType& rM = *mpM;
        TSystemMatrixType& rMi = *mpMi;
        TSolutionVectorType& rRHS  = *mpRHS;
        TSystemMatrixType& rC  = *mpC;

        if (this->GetEchoLevel() == 2) //if it is needed to print the debug info
        {
            KRATOS_INFO("RHS") << "RHS  = " << rRHS << std::endl;
        }
        else if (this->GetEchoLevel() == 3) //if it is needed to print the debug info
        {
            KRATOS_INFO("K") << "SystemMatrix = " << rK << std::endl;
            KRATOS_INFO("M")  << "Mass Matrix = " << mpM << std::endl;
            KRATOS_INFO("C")  << "Damping Matrix = " << mpC << std::endl;
            KRATOS_INFO("RHS") << "RHS  = " << rRHS << std::endl;
        }
        std::stringstream matrix_market_name;
        matrix_market_name << "K_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_name.str()).c_str(), rK, false);

        std::stringstream matrix_market_ki_name;
        matrix_market_ki_name << "Ki_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_ki_name.str()).c_str(), rKi, false);

        std::stringstream matrix_market_mass_name;
        matrix_market_mass_name << "M_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_mass_name.str()).c_str(), rM, false);

        std::stringstream matrix_market_mi_name;
        matrix_market_mi_name << "Mi_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_mi_name.str()).c_str(), rMi, false);

        std::stringstream matrix_market_damping_name;
        matrix_market_damping_name << "C_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << ".mm";
        ComplexSparseSpaceType::WriteMatrixMarketMatrix((char *)(matrix_market_damping_name.str()).c_str(), rC, false);

        std::stringstream matrix_market_vectname;
        matrix_market_vectname << "RHS_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << ".mm.rhs";
        TSparseSpace::WriteMatrixMarketVector((char *)(matrix_market_vectname.str()).c_str(), rRHS);
    }

    template <typename SparseSpaceType>
    void MatrixOutput(typename SparseSpaceType::MatrixType matrix, const std::string name)
    {
        std::stringstream matrix_market_name;
        matrix_market_name << name << ".mm";
        SparseSpaceType::WriteMatrixMarketMatrix((char *)(matrix_market_name.str()).c_str(), matrix, false);
    }

    template <typename SparseSpaceType>
    void VectorOutput(typename SparseSpaceType::VectorType vector, const std::string name)
    {
        std::stringstream matrix_market_name;
        matrix_market_name << name << ".mm.rhs";
        SparseSpaceType::WriteMatrixMarketVector((char *)(matrix_market_name.str()).c_str(), vector);
    }

    void MaterialSettingsInfo()
    {
        for( auto& setting : mFrequencyDependentSettings ) {
            auto rsetting = *setting;
            KRATOS_INFO("Frequency dependent material") << "BUILD_LEVEL=" << rsetting.build_level <<
                ", FACTOR=" << rsetting.factor << " , Initialized? " << rsetting.initialized << std::endl;
        }
    }

    std::vector<unsigned int> GetFixedDofs()
    {
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();
        std::vector<unsigned int> fixed_dofs;
        DirichletUtility::GetDirichletConstraints<typename TBuilderAndSolverType::DofsArrayType>(p_builder_and_solver->GetDofSet(), fixed_dofs);

        return fixed_dofs;
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

  private:
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

    void InitializeFrequencyDependentMatrices(FrequencyDependentSettings<TSparseSpace>& setting,
        const std::vector<unsigned int>& FixedDofSet)
    {
        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();

        KRATOS_ERROR_IF_NOT( setting.has_matrix_contribution || setting.has_vector_contribution ) <<
            "No contribution type has been set for frequency dependent contribution with build level " <<
            setting.build_level << std::endl;

        // build on provided model part
        ModelPart& r_model_part = *setting.model_part;
        TSystemMatrixType& r_matrix = *setting.matrix;
        TSystemVectorType& r_vector = *setting.vector;

        p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, setting.matrix, setting.vector, setting.vector,
            BaseType::GetModelPart());
        r_model_part.GetProcessInfo()[BUILD_LEVEL] = setting.build_level;

        p_builder_and_solver->Build(p_scheme, r_model_part, r_matrix, r_vector);

        // remove temporaries
        if( !setting.has_matrix_contribution ) {
            TSparseSpace::Clear(setting.matrix);
        }
        if( !setting.has_vector_contribution ) {
            TSparseSpace::Clear(setting.vector);
        }

        setting.initialized = true;
        std::cout << "initialized frequency dependent matrix\n";
    }

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

  protected:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{
    typename TLinearSolver::Pointer mpLinearSolver; /// The pointer to the linear solver considered
    typename TSchemeType::Pointer mpScheme; /// The pointer to the time scheme employed
    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver; /// The pointer to the builder and solver
    ComplexLinearSolverType::Pointer mpComplexLinearSolver;

    TSolutionVectorPointerType mpRHS; /// The RHS vector
    TSolutionVectorPointerType mpDx; /// The solution vector
    TSolutionMatrixPointerType mpA; /// The system matrix (dynamic stiffness matrix)
    TSystemMatrixPointerType mpK; /// The stiffness matrix (real part)
    TSystemMatrixPointerType mpKi; /// The stiffness matrix (imaginary part)
    TSystemMatrixPointerType mpM; /// The mass matrix (real part)
    TSystemMatrixPointerType mpMi; /// The mass matrix (imaginary part)
    TSystemMatrixPointerType mpC; /// The damping matrix
    TSystemVectorPointerType mpB; // The constant part of the RHS vector

    bool mReformDofSetAtEachStep;

    bool mSolutionStepIsInitialized; /// Flag to set as initialized the solution step

    bool mInitializeWasPerformed; /// Flag to set as initialized the strategy

    bool mUseModalDamping; // Flag to set if modal damping is used

    bool mUseFrequencyDependentMaterials = false; // Flag to set if Biot material is used

    std::list<FrequencyDependentSettings<TSparseSpace>*> mFrequencyDependentSettings;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void ApplyComplexDirichletConditions()
    {
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();
        std::vector<unsigned int> fixed_dofs;
        DirichletUtility::GetDirichletConstraints<typename TBuilderAndSolverType::DofsArrayType>(p_builder_and_solver->GetDofSet(), fixed_dofs);

        // this should not be hard coded
        DirichletUtility::SetComplexDirichletConditions<TSolutionSpace>(
            BaseType::GetModelPart(),
            *mpA,
            *mpRHS,
            fixed_dofs,
            PRESSURE,
            REAL_PRESSURE,
            IMAG_PRESSURE
        );

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

    FrequencyResponseAnalysisStrategy(const FrequencyResponseAnalysisStrategy &Other){};

    ///@}

}; /* Class MorOfflineStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* MOR_OFFLINE_STRATEGY  defined */