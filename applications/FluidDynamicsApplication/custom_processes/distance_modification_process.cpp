//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "utilities/openmp_utils.h"
#include "processes/find_nodal_h_process.h"

// Application includes
#include "distance_modification_process.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

/* Public functions *******************************************************/
DistanceModificationProcess::DistanceModificationProcess(
    ModelPart& rModelPart,
    const double FactorCoeff, //TODO: Remove it (here for legacy reasons)
    const double DistanceThreshold,
    const bool CheckAtEachStep,
    const bool NegElemDeactivation,
    const bool RecoverOriginalDistance)
    : Process(), mrModelPart(rModelPart) {

    mDistanceThreshold = DistanceThreshold;
    mCheckAtEachStep = CheckAtEachStep;
    mNegElemDeactivation = NegElemDeactivation;
    mRecoverOriginalDistance = RecoverOriginalDistance;
}

DistanceModificationProcess::DistanceModificationProcess(
    ModelPart& rModelPart,
    Parameters& rParameters)
    : Process(), mrModelPart(rModelPart) {

    Parameters default_parameters( R"(
    {
        "model_part_name"                        : "default_model_part_name",
        "distance_factor"                        : 2.0, 
        "distance_threshold"                     : 0.001,
        "continuous_distance"                    : true,
        "check_at_each_time_step"                : false,
        "avoid_almost_empty_elements"            : true,
        "deactivate_full_negative_elements"      : true,
        "recover_original_distance_at_each_step" : false
    }  )" );

    rParameters.ValidateAndAssignDefaults(default_parameters);

    mDistanceThreshold = rParameters["distance_threshold"].GetDouble();
    mContinuousDistance = rParameters["continuous_distance"].GetBool();
    mCheckAtEachStep = rParameters["check_at_each_time_step"].GetBool();
    mAvoidAlmostEmptyElements = rParameters["avoid_almost_empty_elements"].GetBool();
    mNegElemDeactivation = rParameters["deactivate_full_negative_elements"].GetBool();
    mRecoverOriginalDistance = rParameters["recover_original_distance_at_each_step"].GetBool();
}

void DistanceModificationProcess::ExecuteInitialize() {

    KRATOS_TRY;

    // Continuous distance field required variables check 
    if (mContinuousDistance){
        const auto& r_node = *mrModelPart.NodesBegin();
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NODAL_H, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);
    }

    KRATOS_CATCH("");
}

void DistanceModificationProcess::ExecuteBeforeSolutionLoop() {

    KRATOS_TRY;

    // Modify the nodal distance values to avoid bad intersections
    if (mContinuousDistance) {
        // Compute NODAL_H (used for computing the distance tolerance)
        FindNodalHProcess nodal_h_calculator(mrModelPart);
        nodal_h_calculator.Execute();
        // Modify the continuous distance field
        this->ModifyDistance();
    } else {
        // Modify the discontinuous distance field
        this->ModifyDiscontinuousDistance();
    }

    // If proceeds (depending on the formulation), perform the deactivation
    // Deactivates the full negative elements and sets the inner values to 0
    if (mNegElemDeactivation) {
        this->DeactivateFullNegativeElements();
    }

    KRATOS_CATCH("");
}

void DistanceModificationProcess::ExecuteInitializeSolutionStep() {

    if(mCheckAtEachStep == true) {
        DistanceModificationProcess::ExecuteBeforeSolutionLoop();
    }
}

void DistanceModificationProcess::ExecuteFinalizeSolutionStep() {

    if(mRecoverOriginalDistance == true) {
        if (mContinuousDistance){
            this->RecoverOriginalDistance();
        } else {
            this->RecoverOriginalDiscontinuousDistance();
        }
    }
}

/* Protected functions ****************************************************/

void DistanceModificationProcess::ModifyDistance() {

    ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();
    ModelPart::ElementsContainerType& rElements = mrModelPart.Elements();

    // Distance modification
    // Case in where the original distance does not need to be preserved (e.g. CFD)
    if (mRecoverOriginalDistance == false) {
        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(rNodes.size()); ++k) {
            ModelPart::NodesContainerType::iterator itNode = rNodes.begin() + k;
            const double h = itNode->FastGetSolutionStepValue(NODAL_H);
            const double tol_d = mDistanceThreshold*h;
            double& d = itNode->FastGetSolutionStepValue(DISTANCE);

            // Check if the distance values are close to zero
            // If proceeds, set the tolerance as distance value
            if(std::abs(d) < tol_d){
                if (d <= 0.0){
                    d = -tol_d;
                } else {
                    // If selected, avoid almost empty elements
                    if (mAvoidAlmostEmptyElements){
                        d = -tol_d;
                    } else {
                        d = tol_d;
                    }
                }
            }
        }
    }
    // Case in where the original distance needs to be kept to track the interface (e.g. FSI)
    else {
        const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
        std::vector<std::vector<unsigned int>> AuxModifiedDistancesIDs(NumThreads);
        std::vector<std::vector<double>> AuxModifiedDistancesValues(NumThreads);

        #pragma omp parallel shared(AuxModifiedDistancesIDs, AuxModifiedDistancesValues)
        {
            const int ThreadId = OpenMPUtils::ThisThread();             // Get the thread id
            std::vector<unsigned int>   LocalModifiedDistancesIDs;      // Local modified distances nodes id vector
            std::vector<double>      LocalModifiedDistancesValues;      // Local modified distances original values vector

            #pragma omp for
            for (int k = 0; k < static_cast<int>(rNodes.size()); ++k) {
                ModelPart::NodesContainerType::iterator itNode = rNodes.begin() + k;
                const double h = itNode->FastGetSolutionStepValue(NODAL_H);
                const double tol_d = mDistanceThreshold*h;
                double &d = itNode->FastGetSolutionStepValue(DISTANCE);

                // Check if the distance values are close to zero
                // If proceeds, set the tolerance as distance value
                if(std::abs(d) < tol_d){

                    // Store the original distance to be recovered at the end of the step
                    LocalModifiedDistancesIDs.push_back(d);
                    LocalModifiedDistancesValues.push_back(itNode->Id());

                    if (d <= 0.0){
                        d = -tol_d;
                    } else {
                        // If selected, avoid almost empty elements
                        if (mAvoidAlmostEmptyElements){
                            d = -tol_d;
                        } else {
                            d = tol_d;
                        }
                    }
                }
            }

            AuxModifiedDistancesIDs[ThreadId] = LocalModifiedDistancesIDs;
            AuxModifiedDistancesValues[ThreadId] = LocalModifiedDistancesValues;
        }

        mModifiedDistancesIDs = AuxModifiedDistancesIDs;
        mModifiedDistancesValues = AuxModifiedDistancesValues;
    }

    // Syncronize data between partitions (the modified distance has always a lower value)
    mrModelPart.GetCommunicator().SynchronizeCurrentDataToMin(DISTANCE);
}

void DistanceModificationProcess::ModifyDiscontinuousDistance(){

    auto elems_begin = mrModelPart.ElementsBegin();
    const std::size_t n_elems = mrModelPart.NumberOfElements();

    // Distance modification
    if (mRecoverOriginalDistance == false) {
        // Case in where the original distance does not need to be preserved (e.g. CFD)
        for (unsigned int i_elem = 0; i_elem < n_elems; ++i_elem){
            auto it_elem = elems_begin + i_elem;

            // Compute the distance tolerance
            const double tol_d = this->ComputeDiscontinuousDistanceElementTolerance(it_elem);

            // Check if the distance values are close to zero
            Vector &r_elem_dist = it_elem->GetValue(ELEMENTAL_DISTANCES);
            for (unsigned int i_node = 0; i_node < r_elem_dist.size(); ++i_node){
                if (std::abs(r_elem_dist(i_node)) < tol_d){
                    r_elem_dist(i_node) = -tol_d;
                }
            }
        }
    } else {
        // Case in where the original distance needs to be kept to track the interface (e.g. FSI)
        const unsigned int n_threads = OpenMPUtils::GetNumThreads();
        std::vector<std::vector<unsigned int>> mod_dist_elems_ids(n_threads);
        std::vector<std::vector<Vector>> orig_dist_values(n_threads);

        #pragma omp parallel shared(mod_dist_elems_ids, orig_dist_values)
        {
            const int i_thread = OpenMPUtils::ThisThread();         // Get the thread id
            std::vector<unsigned int> local_mod_dist_elems_ids;     // Local modified distances elemental id vector
            std::vector<Vector> local_orig_dist_values;             // Local modified distances original values vector

            #pragma omp for
            for (unsigned int i_elem = 0; i_elem < n_elems; ++i_elem){
                auto it_elem = elems_begin + i_elem;

                // Compute the distance tolerance
                const double tol_d = this->ComputeDiscontinuousDistanceElementTolerance(it_elem);

                bool is_saved = false;
                Vector &r_elem_dist = it_elem->GetValue(ELEMENTAL_DISTANCES);
                for (unsigned int i_node = 0; i_node < r_elem_dist.size(); ++i_node){
                    if (std::abs(r_elem_dist(i_node)) < tol_d){
                        r_elem_dist(i_node) = -tol_d;
                        if (!is_saved){
                            local_mod_dist_elems_ids.push_back(it_elem->Id());
                            local_orig_dist_values.push_back(r_elem_dist);
                        }
                    }
                }

            }

            mod_dist_elems_ids[i_thread] = local_mod_dist_elems_ids;
            orig_dist_values[i_thread] = local_orig_dist_values;
        }

        mModifiedDistancesIDs = mod_dist_elems_ids;
        mModifiedElementalDistancesValues = orig_dist_values;
    }
}

double DistanceModificationProcess::ComputeDiscontinuousDistanceElementTolerance(ModelPart::ElementIterator itElem){
    auto &r_geometry = itElem->GetGeometry();
    const std::size_t work_dim = r_geometry.WorkingSpaceDimension();
    const double elem_size = r_geometry.Area(); // Area (2D) or volume (3D)
    const double tol_d = work_dim == 2 ? mDistanceThreshold*std::sqrt(elem_size) : mDistanceThreshold*std::cbrt(elem_size);
    return tol_d;
}

void DistanceModificationProcess::RecoverOriginalDistance() {
    #pragma omp parallel
    {
        const int ThreadId = OpenMPUtils::ThisThread();
        const std::vector<unsigned int> LocalModifiedDistancesIDs = mModifiedDistancesIDs[ThreadId];
        const std::vector<double> LocalModifiedDistancesValues = mModifiedDistancesValues[ThreadId];

        for(unsigned int i=0; i<LocalModifiedDistancesIDs.size(); ++i) {
            const unsigned int nodeId = LocalModifiedDistancesIDs[i];
            mrModelPart.GetNode(nodeId).FastGetSolutionStepValue(DISTANCE) = LocalModifiedDistancesValues[i];
        }
    }

    // Syncronize data between partitions (the modified distance has always a lower value)
    mrModelPart.GetCommunicator().SynchronizeCurrentDataToMin(DISTANCE);

    // Empty the modified distance vectors
    mModifiedDistancesIDs.resize(0);
    mModifiedDistancesValues.resize(0);
    mModifiedDistancesIDs.shrink_to_fit();
    mModifiedDistancesValues.shrink_to_fit();
}

void DistanceModificationProcess::RecoverOriginalDiscontinuousDistance() {
    #pragma omp parallel
    {
        const int i_thread = OpenMPUtils::ThisThread();
        const std::vector<unsigned int> local_mod_dist_elems_ids = mModifiedDistancesIDs[i_thread];
        const std::vector<Vector> local_orig_dist_values = mModifiedElementalDistancesValues[i_thread];

        for(unsigned int i = 0; i < local_mod_dist_elems_ids.size(); ++i) {
            const unsigned int i_elem = local_mod_dist_elems_ids[i];
            mrModelPart.GetElement(i_elem).SetValue(ELEMENTAL_DISTANCES,local_orig_dist_values[i]);
        }
    }

    // Empty the modified distance vectors
    mModifiedDistancesIDs.resize(0);
    mModifiedElementalDistancesValues.resize(0);
    mModifiedDistancesIDs.shrink_to_fit();
    mModifiedElementalDistancesValues.shrink_to_fit();
}

void DistanceModificationProcess::DeactivateFullNegativeElements() {

    ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();
    ModelPart::ElementsContainerType& rElements = mrModelPart.Elements();

    // Initialize the EMBEDDED_IS_ACTIVE variable flag to 0
    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>(rNodes.size()); ++i_node){
        ModelPart::NodesContainerType::iterator it_node = rNodes.begin() + i_node;
        it_node->SetValue(EMBEDDED_IS_ACTIVE, 0);
    }

    // Deactivate those elements whose negative distance nodes summation is equal to their number of nodes
    #pragma omp parallel for
    for (int k = 0; k < static_cast<int>(rElements.size()); ++k){
        unsigned int n_neg = 0;
        ModelPart::ElementsContainerType::iterator itElement = rElements.begin() + k;
        auto& rGeometry = itElement->GetGeometry();

        // Check the distance function sign at the element nodes
        for (unsigned int i_node=0; i_node<rGeometry.size(); i_node++){
            if (rGeometry[i_node].FastGetSolutionStepValue(DISTANCE) < 0.0){
                n_neg++;
            }
        }

        (n_neg == rGeometry.size()) ? itElement->Set(ACTIVE, false) : itElement->Set(ACTIVE, true);

        // If the element is ACTIVE, all its nodes are active as well
        if (itElement->Is(ACTIVE)){
            for (unsigned int i_node = 0; i_node < rGeometry.size(); ++i_node){
                int& activation_index = rGeometry[i_node].GetValue(EMBEDDED_IS_ACTIVE);
                #pragma omp atomic
                activation_index += 1;
            }
        }
    }

    // Synchronize the EMBEDDED_IS_ACTIVE variable flag
    mrModelPart.GetCommunicator().AssembleNonHistoricalData(EMBEDDED_IS_ACTIVE);

    // Set to zero and fix the DOFs in the remaining inactive nodes
    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>(rNodes.size()); ++i_node){
        ModelPart::NodesContainerType::iterator it_node = rNodes.begin() + i_node;
        if (it_node->GetValue(EMBEDDED_IS_ACTIVE) == 0){
            // Fix the nodal DOFs
            it_node->Fix(PRESSURE);
            it_node->Fix(VELOCITY_X);
            it_node->Fix(VELOCITY_Y);
            it_node->Fix(VELOCITY_Z);
            // Set to zero the nodal DOFs
            it_node->FastGetSolutionStepValue(PRESSURE) = 0.0;
            it_node->FastGetSolutionStepValue(VELOCITY) = ZeroVector(3);
        }
    }
}

/* Private functions ****************************************************/

};  // namespace Kratos.
