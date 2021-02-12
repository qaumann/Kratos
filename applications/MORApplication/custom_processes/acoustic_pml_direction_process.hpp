//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt

#if !defined(KRATOS_ACOUSTIC_PML_DIRECTION_PROCESS_H_INCLUDED)
#define KRATOS_ACOUSTIC_PML_DIRECTION_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "processes/compute_nodal_gradient_process.h"
#include "includes/model_part.h"
#include "includes/define.h"
#include "mor_application_variables.h"

#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "utilities/variable_utils.h"



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

/**
 * @class AcousticPMLDirectionProcess
 *
 * @ingroup MORApplication
 *
 * @brief This method provides PML damping vector field.
 * @details Based on solving a convection problem from the inner interface of the PML to the outer boundary, the PML damping direction is calculated for every node in the PML.
*/
template< class TSparseSpace, class TDenseSpace, class TLinearSolver >
class AcousticPMLDirectionProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{
    typedef std::size_t SizeType;

    typedef Node < 3 > NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;

    /// The definitionof the geometry
    typedef Geometry<NodeType> GeometryType;

    typedef std::vector<NodeTypePointer> NodeVector;

    typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;

    typedef SolvingStrategy< TSparseSpace, TDenseSpace, TLinearSolver > SolvingStrategyType;

    typedef typename TDenseSpace::MatrixType TDenseMatrixType;

    typedef typename TDenseSpace::MatrixPointerType TDenseMatrixPointerType;


    /// Pointer definition of MonolithicMappingProcess
    KRATOS_CLASS_POINTER_DEFINITION(AcousticPMLDirectionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AcousticPMLDirectionProcess(
        ModelPart& rModelPartPML,
        ModelPart& rModelPartInterface,
        ModelPart& rModelPartBoundary,
        typename TLinearSolver::Pointer plinear_solver
        ) : mrModelPartPML(rModelPartPML),
        mrModelPartInterface(rModelPartInterface),
        mrModelPartBoundary(rModelPartBoundary)
    {
        KRATOS_TRY

        // Check that there is at least one element and node in the model
        const auto n_nodes = mrModelPartPML.NumberOfNodes();
        const auto n_elems = mrModelPartPML.NumberOfElements();

        KRATOS_ERROR_IF(n_nodes == 0) << "The model has no nodes." << std::endl;
        KRATOS_ERROR_IF(n_elems == 0) << "The model has no elements." << std::endl;

        // Set build level

        mrModelPartPML.GetProcessInfo()[BUILD_LEVEL] = 401;

        // Generate a linear strategy
        typename SchemeType::Pointer pScheme = Kratos::make_shared< ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace > >();
        typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;

        bool CalculateReactions = false;
        bool ReformDofAtEachIteration = false;
        bool CalculateNormDxFlag = false;

        BuilderSolverTypePointer pBuilderSolver = Kratos::make_shared< ResidualBasedBlockBuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver > >(plinear_solver);
        mpSolvingStrategy = Kratos::make_unique< ResidualBasedLinearStrategy<TSparseSpace,TDenseSpace,TLinearSolver > >(
            mrModelPartPML,
            pScheme,
            // plinear_solver,
            pBuilderSolver,
            CalculateReactions,
            ReformDofAtEachIteration,
            CalculateNormDxFlag);

        mpSolvingStrategy->SetEchoLevel(0);

        //TODO: check flag DO_EXPENSIVE_CHECKS
        mpSolvingStrategy->Check();

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~AcousticPMLDirectionProcess() override = default;


    ///@}
    ///@name Operators
    ///@{

    void operator()(){
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
        KRATOS_TRY

        ModelPart::NodesContainerType& rInterfaceNodes = mrModelPartInterface.Nodes();
        ModelPart::NodesContainerType& rBoundaryNodes = mrModelPartBoundary.Nodes();
        VariableUtils().SetFlag(BOUNDARY, true, rBoundaryNodes);



        #pragma omp parallel for
        for ( ModelPart::NodesContainerType::iterator it_node = rInterfaceNodes.begin();
                it_node != rInterfaceNodes.end() ; ++it_node)
        {
            it_node->SetValue(PRESCRIBED_POTENTIAL, -1);        // the value used for the rhs
            it_node->GetSolutionStepValue(PRESSURE, 0) = -1;    // to have a consistent result
            it_node->Fix(PRESSURE);                             // dof must be fixed
        }

        #pragma omp parallel for
        for ( ModelPart::NodesContainerType::iterator it_node = rBoundaryNodes.begin();
                it_node != rBoundaryNodes.end() ; ++it_node)
        {
            it_node->SetValue(PRESCRIBED_POTENTIAL, 1);         // the value used for the rhs
            it_node->GetSolutionStepValue(PRESSURE, 0) = 1;     // to have a consistent result
            it_node->Fix(PRESSURE);                             // dof must be fixed

        }


        mpSolvingStrategy->Solve();

        typedef ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable> GradientType;
        GradientType process = GradientType(mrModelPartPML, PRESSURE, ABSORBTION_VECTOR, NODAL_AREA, false);
        process.Execute();

        // normalize gradients
        #pragma omp parallel for
        for( auto& it_node : mrModelPartPML.Nodes() ) {
            array_1d<double, 3>& grad = it_node.GetValue(ABSORBTION_VECTOR);
            const double norm = norm_2( grad );
            if( std::abs(norm) > std::numeric_limits<double>::epsilon() ) {
                std::for_each( grad.begin(), grad.end(), [norm] (double &a) {a /= norm;});
            }
        }

        // use mean of adjacent nodes for zero gradients
        #pragma omp parallel for
        for( auto& it_elem : mrModelPartPML.Elements() ) {
            for( auto& it_node : it_elem.GetGeometry() ) {
                array_1d<double, 3> grad = it_node.GetValue(ABSORBTION_VECTOR);
                double norm = norm_2( grad );
                if( std::abs(norm) < std::numeric_limits<double>::epsilon() ) {
                    for( auto& it_elem_node : it_elem.GetGeometry() ) {
                        grad += it_elem_node.GetValue(ABSORBTION_VECTOR);
                    }
                    norm = norm_2( grad );
                    std::for_each( grad.begin(), grad.end(), [norm] (double &a) {a /= norm;});
                    it_node.SetValue(ABSORBTION_VECTOR, grad);
                }
            }
        }

        // unfix all dofs and reset PRESSURE to be able to start computation with a clean model
        #pragma omp parallel for
        for( auto& it_node : mrModelPartPML.Nodes() ) {
            it_node.Free(PRESSURE);
            it_node.SetValue(PRESCRIBED_POTENTIAL, it_node.GetSolutionStepValue(PRESSURE, 0));
            it_node.GetSolutionStepValue(PRESSURE, 0) = 0;
        }


        // compute local pml width at boundary nodes
        for( auto& boundary_node : mrModelPartBoundary.Nodes() ) {

                double min_width = 10e7;
                double x_b = boundary_node.X();
                double y_b = boundary_node.Y();
                double z_b = boundary_node.Z();

                for( auto& interface_node : mrModelPartInterface.Nodes() ) {
                    double x_i = interface_node.X();
                    double y_i = interface_node.Y();
                    double z_i = interface_node.Z();
                    double width = std::sqrt(std::pow(x_b - x_i, 2) + std::pow(y_b - y_i, 2) + std::pow(z_b - z_i, 2));
                    if (width < min_width) {
                        min_width = width;
                    }

                }
            boundary_node.SetValue(LOCAL_PML_WIDTH, min_width);
        }

        // compute distance of node to interface, set local pml width of nodes
        for( auto& pml_node : mrModelPartPML.Nodes() ) {
            double x = pml_node.X();
            double y = pml_node.Y();
            double z = pml_node.Z();
            double min_dist_i = 10e7;
            double min_dist_b = 10e7;
            double width;

            for( auto& boundary_node : mrModelPartBoundary.Nodes() ) {
                double x_b = boundary_node.X();
                double y_b = boundary_node.Y();
                double z_b = boundary_node.Z();
                double dist = std::sqrt(std::pow(x - x_b, 2) + std::pow(y - y_b, 2) + std::pow(z - z_b, 2));
                if (dist < min_dist_b) {
                    min_dist_b = dist;
                    width = boundary_node.GetValue(LOCAL_PML_WIDTH);
                }
            }

            for( auto& interface_node : mrModelPartInterface.Nodes() ) {

                double x_i = interface_node.X();
                double y_i = interface_node.Y();
                double z_i = interface_node.Z();

                double dist = std::sqrt(std::pow(x - x_i, 2) + std::pow(y - y_i, 2) + std::pow(z - z_i, 2));
                if (dist < min_dist_i) {
                    min_dist_i = dist;
                }


            }

            pml_node.SetValue(IMAG_DISTANCE, min_dist_i);
            if( !pml_node.Is(BOUNDARY) ){
                pml_node.SetValue(LOCAL_PML_WIDTH, width);
            }


        }
        KRATOS_CATCH("")
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "AcousticAcousticPMLDirectionProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AcousticAcousticPMLDirectionProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

private:
    ///@}
    ///@name Member Variables
    ///@{
    ModelPart& mrModelPartPML;
    ModelPart& mrModelPartInterface;
    ModelPart& mrModelPartBoundary;
    typename SolvingStrategyType::UniquePointer mpSolvingStrategy;


};

}
#endif /* KRATOS_ACOUSTIC_PML_DIRECTION_PROCESS_H_INCLUDED defined */
