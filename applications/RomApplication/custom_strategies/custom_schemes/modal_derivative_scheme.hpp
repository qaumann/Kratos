//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Altug Emiroglu, https://github.com/emiroglu
//


#if !defined(KRATOS_MODAL_DERIVATIVE_SCHEME )
#define  KRATOS_MODAL_DERIVATIVE_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/dof.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "solving_strategies/schemes/scheme.h"

/* Application includes */
#include "rom_application_variables.h"
#include "custom_utilities/rom_finite_difference_utility.h"

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
 * @class ModalDerivativeScheme
 * @ingroup KratosCore
 * @brief This class provides the implementation of the basic tasks that are needed by the solution strategy.
 * @details It is intended to be the place for tailoring the solution strategies to problem specific tasks.
 * @tparam TSparseSpace The sparse space considered
 * @tparam TDenseSpace The dense space considered
 * @author Altug Emiroglu
 */
template<class TSparseSpace,
         class TDenseSpace //= DenseSpace<double>
         >
class ModalDerivativeScheme : public Scheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ModalDerivativeScheme
    KRATOS_CLASS_POINTER_DEFINITION(ModalDerivativeScheme);

    // Base type definition
    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;

    /// Data type definition
    typedef typename BaseType::TDataType TDataType;
    /// Matrix type definition
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    /// Vector type definition
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    /// Local system matrix type definition
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    /// Local system vector type definition
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    /// DoF type definition
    typedef Dof<double> TDofType;
    /// DoF array type definition
    typedef typename BaseType::DofsArrayType DofsArrayType;
    /// DoF iterator type definition
    typedef typename BaseType::DofIterator DofIterator;
    /// DoF constant iterator type definition
    typedef typename BaseType::DofConstantIterator DofConstantIterator;

    /// Elements containers definition
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    /// Conditions containers definition
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default Constructor
     * @details Initiliazes the flags
     */
    /// Constructor.
    ModalDerivativeScheme(Variable<double> DerivativeParameter, double FiniteDifferenceStepSize, bool FiniteDifferenceTypeFlag=false)
    : 
    Scheme<TSparseSpace,TDenseSpace>(), 
    mDerivativeParameter(DerivativeParameter),
    mFiniteDifferenceStepSize(FiniteDifferenceStepSize), 
    mFiniteDifferenceTypeFlag(FiniteDifferenceTypeFlag)
    {

        if (mDerivativeParameter == DENSITY)
            mDerivativeMatrixType = true;
        else if (mDerivativeParameter == MODAL_COORDINATE || mDerivativeParameter == YOUNG_MODULUS || mDerivativeParameter == POISSON_RATIO)
            mDerivativeMatrixType = false;
        else
            KRATOS_ERROR;
    }

    /// Destructor.
    ~ModalDerivativeScheme() override 
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief unction to be called when it is needed to initialize an iteration. It is designed to be called at the beginning of each non linear iteration
     * @note Take care: the elemental function with the same name is NOT called here.
     * @warning Must be defined in derived classes
     * @details The function is called in the builder for memory efficiency
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    void InitializeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Definition of the first element iterator
        const auto it_elem_begin = rModelPart.ElementsBegin();

        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); ++i) {
            auto it_elem = it_elem_begin + i;
            it_elem->InitializeNonLinearIteration(r_current_process_info);
        }

        // Definition of the first condition iterator
        const auto it_cond_begin = rModelPart.ConditionsBegin();

        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); ++i) {
            auto it_cond = it_cond_begin + i;
            it_cond->InitializeNonLinearIteration(r_current_process_info);
        }

        // Definition of the first constraint iterator
        const auto it_const_begin = rModelPart.MasterSlaveConstraintsBegin();

        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.MasterSlaveConstraints().size()); ++i) {
            auto it_const = it_const_begin + i;
            it_const->InitializeNonLinearIteration(r_current_process_info);
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief It initializes a non-linear iteration (for an individual condition)
     * @warning Must be defined in derived classes
     * @param pCurrentElement The element to compute
     * @param rCurrentProcessInfo The current process info instance
     */
    void InitializeNonLinearIteration(
        Element::Pointer pCurrentElement,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        (pCurrentElement)->InitializeNonLinearIteration(rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    /**
     * @brief It initializes a non-linear iteration (for an individual condition)
     * @warning Must be defined in derived classes
     * @param rCurrentCondition The condition to compute
     * @param rCurrentProcessInfo The current process info instance
     */
    void InitializeNonLinearIteration(
        Condition::Pointer pCurrentCondition,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        (pCurrentCondition)->InitializeNonLinearIteration(rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    // Element contributions
    void CalculateSystemContributions(
        Element& rElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        this->CalculateLHSContribution(rElement, rLHS_Contribution, rEquationId, rCurrentProcessInfo);

        this->CalculateRHSContribution(rElement, rRHS_Contribution, rEquationId, rCurrentProcessInfo);

        rElement.EquationIdVector(rEquationId,rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateLHSContribution(
        Element& rElement,
        LocalSystemMatrixType& rLHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        if (rCurrentProcessInfo[BUILD_LEVEL] == 1 ) // Mass matrix
            rElement.CalculateMassMatrix(rLHS_Contribution, rCurrentProcessInfo);
        else if (rCurrentProcessInfo[BUILD_LEVEL] == 2) // Stiffness matrix
            rElement.CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
        else 
            KRATOS_ERROR <<"Invalid BUILD_LEVEL: " << rCurrentProcessInfo[BUILD_LEVEL] << std::endl;
        
        rElement.EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateRHSContribution(
        Element& rElement,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY
        
        // Derivative of basis_i
        std::size_t basis_i = rCurrentProcessInfo[BASIS_I];
        
        // Create PhiElemental
        LocalSystemVectorType PhiElemental;
        std::size_t num_element_nodal_dofs = 0;
        for (auto& node_i : rElement.GetGeometry())
            num_element_nodal_dofs += node_i.GetDofs().size();
        PhiElemental.resize(num_element_nodal_dofs);
        
        // Get PhiElemental
        std::size_t dof_ctr = 0;
        for (auto& node_i : rElement.GetGeometry()) {
            auto& node_i_dofs = node_i.GetDofs();
            
            const Matrix *pPhiNodal = &(node_i.GetValue(ROM_BASIS));
            for (std::size_t dof_idx = 0; dof_idx < node_i_dofs.size(); dof_idx++){
                PhiElemental[dof_ctr + dof_idx] = (*pPhiNodal)(dof_idx, basis_i);
            }
            dof_ctr += node_i_dofs.size();
        }
        
        // Build RHS contribution

        // Initialize RHS contribution
        rRHS_Contribution.clear();
        std::vector<Dof<double>::Pointer> rElementalDofList;
        rElement.GetDofList(rElementalDofList, rCurrentProcessInfo);
        std::size_t elementalDofSize = rElementalDofList.size();
        rRHS_Contribution.resize(elementalDofSize, 0);
        
        // Compute element LHS derivative
        Matrix element_matrix_derivative;
        element_matrix_derivative.resize(elementalDofSize,elementalDofSize,false);
        if (mDerivativeParameter == MODAL_COORDINATE)
        {   // Modal parameter            

            // Derivative wrt basis_j
            std::size_t basis_j = rCurrentProcessInfo[BASIS_J];

            if (mFiniteDifferenceTypeFlag) // Central Difference
                this->CentralDifferencingWithBasis(rElement, element_matrix_derivative, basis_j, rCurrentProcessInfo);
            else // Forward difference
                this->ForwardDifferencingWithBasis(rElement, element_matrix_derivative, basis_j, rCurrentProcessInfo);

        }
        else if (mDerivativeParameter == DENSITY || mDerivativeParameter == YOUNG_MODULUS || mDerivativeParameter == POISSON_RATIO)
        {   // Mass parameter
            
            if (mFiniteDifferenceTypeFlag) // Central Difference
                this->CentralDifferencingWithMaterialParameter(rElement, element_matrix_derivative, mDerivativeParameter, rCurrentProcessInfo);
            else // Forward difference
                this->ForwardDifferencingWithMaterialParameter(rElement, element_matrix_derivative, mDerivativeParameter, rCurrentProcessInfo);

        }

        // Compute RHS contribution
        if (element_matrix_derivative.size1() != PhiElemental.size()){
            // TODO: this is a workaround. use a map as in ROM analysis
            // retrieve only displacement dofs from PhiElemental
            unsigned int disp_shifter = 3;
            Vector tmpPhiElemental;
            tmpPhiElemental.resize(PhiElemental.size() / 2);
            for (std::size_t iNode = 0; iNode < rElement.GetGeometry().size(); iNode++)
            {
                for (std::size_t iXYZ = 0; iXYZ < 3; iXYZ++)
                {
                    tmpPhiElemental[iNode * 3 + iXYZ] = PhiElemental[iNode * 6 + disp_shifter + iXYZ];
                }
            }
            rRHS_Contribution -= prod(element_matrix_derivative, tmpPhiElemental);
        }
        else
        {
            rRHS_Contribution -= prod(element_matrix_derivative, PhiElemental);
        }

        rElement.EquationIdVector(EquationId,rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
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

    // This function locks element nodes for finite differencing
    void LockElementNodes(Element& rElement)
    {
        // Loop over element nodes
        for (auto& node_i : rElement.GetGeometry()) 
            node_i.SetLock();
    }

    // This function locks element nodes for finite differencing
    void UnlockElementNodes(Element& rElement)
    {
        // Loop over element nodes
        for (auto& node_i : rElement.GetGeometry()) 
            node_i.UnSetLock();
    }

    void ForwardDifferencingWithMaterialParameter(
        Element& rElement,
        Matrix& rElementMatrixDerivative,
        const Variable<double>& MaterialParameter,
        const ProcessInfo& rCurrentProcessInfo)
    {

        if ( rElement.GetProperties().Has(MaterialParameter) )
        {

            Matrix element_matrix_initial;
            Matrix element_matrix_p_perturbed;

            // Compute initial matrix
            if (mDerivativeMatrixType)
                rElement.CalculateMassMatrix(element_matrix_initial, rCurrentProcessInfo);
            else
                rElement.CalculateLeftHandSide(element_matrix_initial, rCurrentProcessInfo);
            
            // Save property pointer
            Properties::Pointer p_global_properties = rElement.pGetProperties();
            const double initial_property_value = p_global_properties->GetValue(MaterialParameter);

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
            rElement.SetProperties(p_local_property);
            rElement.Initialize(rCurrentProcessInfo);
            
            // Positive perturbation
            // perturb the variable
            p_local_property->SetValue(MaterialParameter, (initial_property_value + mFiniteDifferenceStepSize));
            rElement.Initialize(rCurrentProcessInfo);
            // Compute element matrix after perturbation
            if (mDerivativeMatrixType)
                rElement.CalculateMassMatrix(element_matrix_p_perturbed, rCurrentProcessInfo);
            else
                rElement.CalculateLeftHandSide(element_matrix_p_perturbed, rCurrentProcessInfo);

            // Reset value
            p_local_property->SetValue(MaterialParameter, initial_property_value);
            rElement.Initialize(rCurrentProcessInfo);
            
            // Give element original properties back
            rElement.SetProperties(p_global_properties);
            rElement.Initialize(rCurrentProcessInfo);

            noalias(rElementMatrixDerivative) = (element_matrix_p_perturbed - element_matrix_initial) / (mFiniteDifferenceStepSize);
        }
    }

    void CentralDifferencingWithMaterialParameter(
        Element& rElement,
        Matrix& rElementMatrixDerivative,
        const Variable<double>& MaterialParameter,
        const ProcessInfo& rCurrentProcessInfo)
    {

        if ( rElement.GetProperties().Has(MaterialParameter) )
        {

            Matrix element_matrix_p_perturbed;
            Matrix element_matrix_n_perturbed;

            // Save property pointer
            Properties::Pointer p_global_properties = rElement.pGetProperties();
            const double initial_property_value = rElement.GetProperties()[MaterialParameter];

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
            rElement.SetProperties(p_local_property);
            rElement.Initialize(rCurrentProcessInfo);

            // Positive perturbation
            // perturb the variable
            p_local_property->SetValue(MaterialParameter, (initial_property_value + mFiniteDifferenceStepSize));
            rElement.Initialize(rCurrentProcessInfo);
            // Compute element matrix after perturbation
            if (mDerivativeMatrixType)
                rElement.CalculateMassMatrix(element_matrix_p_perturbed, rCurrentProcessInfo);
            else
                rElement.CalculateLeftHandSide(element_matrix_p_perturbed, rCurrentProcessInfo);

            // Negative perturbation
            // perturb the variable
            p_local_property->SetValue(MaterialParameter, (initial_property_value - mFiniteDifferenceStepSize));
            rElement.Initialize(rCurrentProcessInfo);
            // Compute element matrix after perturbation
            if (mDerivativeMatrixType)
                rElement.CalculateMassMatrix(element_matrix_n_perturbed, rCurrentProcessInfo);
            else
                rElement.CalculateLeftHandSide(element_matrix_n_perturbed, rCurrentProcessInfo);

            // Reset value
            p_local_property->SetValue(MaterialParameter, initial_property_value);
            rElement.Initialize(rCurrentProcessInfo);

            noalias(rElementMatrixDerivative) = (element_matrix_p_perturbed - element_matrix_n_perturbed) / (2.0*mFiniteDifferenceStepSize);

            // Give element original properties back
            p_local_property->SetValue(MaterialParameter, initial_property_value);
            rElement.SetProperties(p_global_properties);
        }
    }

    void ForwardDifferencingWithBasis(
         Element&rElement,
         Matrix& rElementMatrixDerivative,
         const std::size_t basis_j,
         const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        // Lock element nodes for OMP parallelism
        this->LockElementNodes(rElement);

        std::size_t matrix_size = rElementMatrixDerivative.size1();

        // Positive perturbation
        Matrix element_matrix_p_perturbed;
        element_matrix_p_perturbed.resize(matrix_size, matrix_size, false);
        this->PerturbElementWithBasis(rElement, 1.0, basis_j, rCurrentProcessInfo);
        rElement.CalculateLeftHandSide(element_matrix_p_perturbed, rCurrentProcessInfo);

        // Reset perturbation
        this->PerturbElementWithBasis(rElement, -1.0, basis_j, rCurrentProcessInfo);

        // Neutral state
        Matrix element_matrix;
        element_matrix.resize(matrix_size, matrix_size, false);
        rElement.CalculateLeftHandSide(element_matrix, rCurrentProcessInfo);

        // Compute ElementMatrixDerivative
        noalias(rElementMatrixDerivative) = (element_matrix_p_perturbed - element_matrix) / mFiniteDifferenceStepSize;

        // Unlock element nodes
        this->UnlockElementNodes(rElement);

        KRATOS_CATCH("");
    }

    void CentralDifferencingWithBasis(
         Element&rElement,
         Matrix& rElementMatrixDerivative,
         const std::size_t basis_j,
         const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        // Lock element nodes for OMP parallelism
        this->LockElementNodes(rElement);

        std::size_t matrix_size = rElementMatrixDerivative.size1();

        // Positive perturbation
        Matrix element_matrix_p_perturbed;
        element_matrix_p_perturbed.resize(matrix_size, matrix_size, false);
        this->PerturbElementWithBasis(rElement, 1.0, basis_j, rCurrentProcessInfo);
        rElement.CalculateLeftHandSide(element_matrix_p_perturbed, rCurrentProcessInfo);

        // Negative perturbation
        Matrix element_matrix_n_perturbed;
        element_matrix_n_perturbed.resize(matrix_size, matrix_size, false);
        this->PerturbElementWithBasis(rElement, -2.0, basis_j, rCurrentProcessInfo);
        rElement.CalculateLeftHandSide(element_matrix_n_perturbed, rCurrentProcessInfo);

        // Reset perturbation
        this->PerturbElementWithBasis(rElement, 1.0, basis_j, rCurrentProcessInfo);

        // Compute LHS derivative
        noalias(rElementMatrixDerivative) = (element_matrix_p_perturbed - element_matrix_n_perturbed) / (2.0 * mFiniteDifferenceStepSize);

        this->UnlockElementNodes(rElement);

        KRATOS_CATCH("");
    }

    // This function perturbs the element with the vector with given eigenvector index basis_j
    void PerturbElementWithBasis(
        Element& rElement,
        const double Step,
        const std::size_t basis_j,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        // Initialize is necessary for updating the section properties of shell elements
        rElement.InitializeNonLinearIteration(rCurrentProcessInfo);

        // Loop over element nodes
        for (auto& node_i : rElement.GetGeometry()) {
            auto& node_i_dofs = node_i.GetDofs();
            // Loop over nodal DOFs
            auto it_dof_i = node_i_dofs.begin();
            const int node_i_num_dofs = node_i_dofs.size();
            for (std::size_t dof_idx = 0; dof_idx < node_i_dofs.size(); dof_idx++){

                // Compute and assign the perturbation
                const double dof_perturbation = Step*mFiniteDifferenceStepSize*node_i.GetValue(ROM_BASIS)(dof_idx, basis_j);

                // Some elements need solution step value while others need the current coordinate.
                // Thus perturb all
                // Update solution step value
                (*it_dof_i)->GetSolutionStepValue() += dof_perturbation;

                // Update current nodal coordinates
                if (node_i_num_dofs > 3 && dof_idx > 2) // disp dofs when rots included
                    node_i.Coordinates()[dof_idx-3] += dof_perturbation;
                else if (node_i_num_dofs < 4 ) // disp dofs when rots excluded
                    node_i.Coordinates()[dof_idx] += dof_perturbation;
                
                // Increment the dof iterator
                ++it_dof_i;
            }
        }
        
        // Finalize is necessary for updating the section properties of shell elements
        rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
    }

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

    bool mSchemeIsInitialized;

    bool mDerivativeMatrixType;

    Variable<double> mDerivativeParameter;

    double mFiniteDifferenceStepSize;

    bool mFiniteDifferenceTypeFlag;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class ModalDerivativeScheme

} // namespace Kratos.

#endif /* KRATOS_MODAL_DERIVATIVE_SCHEME  defined */

