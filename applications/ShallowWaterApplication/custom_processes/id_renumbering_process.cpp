//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


// System includes


// External includes


// Project includes
#include "id_renumbering_process.h"


namespace Kratos
{

IdRenumberingProcess::IdRenumberingProcess(Model& rThisModel) : mrModel(rThisModel)
{
    mModelPartNames = GetRootModelPartNames();
}

IdRenumberingProcess::IdRenumberingProcess(Model& rThisModel, const StringVectorType& rModelPartNames)
 : mrModel(rThisModel)
{
    mModelPartNames = rModelPartNames;
}

void IdRenumberingProcess::RenumberNodes()
{
    // The new absolute unique_id
    IndexType unique_id = 0;
    mOriginNodesIds.clear();

    // Loop the model parts
    for (const auto& name : mModelPartNames)
    {
        ModelPart& model_part = mrModel.GetModelPart(name);
        for (int i = 0; i < static_cast<int>(model_part.NumberOfNodes()); ++i)
        {
            auto it_node = model_part.NodesBegin() + i;
            mOriginNodesIds[++unique_id] = it_node->Id();
            it_node->SetId(unique_id);
        }
    }
}

void IdRenumberingProcess::RenumberElements()
{
    // The new absolute unique_id
    IndexType unique_id = 0;
    mOriginElementsIds.clear();

    // Loop the model parts
    for (const auto& name : mModelPartNames)
    {
        ModelPart& model_part = mrModel.GetModelPart(name);
        for (int i = 0; i < static_cast<int>(model_part.NumberOfElements()); ++i)
        {
            auto it_elem = model_part.ElementsBegin() + i;
            mOriginElementsIds[++unique_id] = it_elem->Id();
            it_elem->SetId(unique_id);
        }
    }
}

void IdRenumberingProcess::RenumberConditions()
{
    // The new absolute unique_id
    IndexType unique_id = 0;
    mOriginConditionsIds.clear();

    // Loop the model parts
    for (const auto& name : mModelPartNames)
    {
        ModelPart& model_part = mrModel.GetModelPart(name);
        for (int i = 0; i < static_cast<int>(model_part.NumberOfConditions()); ++i)
        {
            auto it_cond = model_part.ConditionsBegin() + i;
            mOriginConditionsIds[++unique_id] = it_cond->Id();
            it_cond->SetId(unique_id);
        }
    }
}

void IdRenumberingProcess::RestoreNodes()
{
    if(mOriginNodesIds.size() != 0)
    {
        // Loop the model parts
        for (const auto& name : mModelPartNames)
        {
            ModelPart& model_part = mrModel.GetModelPart(name);
            for (int i = 0; i < static_cast<int>(model_part.NumberOfNodes()); ++i)
            {
                auto it_node = model_part.NodesBegin() + i;
                it_node->SetId(mOriginNodesIds[it_node->Id()]);
            }
        }
    }
    // Avoiding the possibility to restore the ids twice
    mOriginNodesIds.clear();
}

void IdRenumberingProcess::RestoreElements()
{
    if (mOriginElementsIds.size() != 0)
    {
        // Loop the model parts
        for (const auto& name : mModelPartNames)
        {
            ModelPart& model_part = mrModel.GetModelPart(name);
            for (int i = 0; i < static_cast<int>(model_part.NumberOfElements()); ++i)
            {
                auto it_elem = model_part.ElementsBegin() + i;
                it_elem->SetId(mOriginElementsIds[it_elem->Id()]);
            }
        }
    }
    // Avoiding the possibility to restore the ids twice
    mOriginElementsIds.clear();
}

void IdRenumberingProcess::RestoreConditions()
{
    if (mOriginConditionsIds.size() != 0)
    {
        // Loop the model parts
        for (const auto& name : mModelPartNames)
        {
            ModelPart& model_part = mrModel.GetModelPart(name);
            for (int i = 0; i < static_cast<int>(model_part.NumberOfConditions()); ++i)
            {
                auto it_cond = model_part.ConditionsBegin() + i;
                it_cond->SetId(mOriginConditionsIds[it_cond->Id()]);
            }
        }
    }
    // Avoiding the possibility to restore the ids twice
    mOriginConditionsIds.clear();
}

StringVectorType IdRenumberingProcess::GetRootModelPartNames() const
{
    StringVectorType root_names;
    StringVectorType all_names = mrModel.GetModelPartNames();
    for (const auto name : all_names)
    {
        // Get the root name of the model part.
        std::string root_name;
        std::string::size_type pos = name.find('.');
        if (pos != std::string::npos) {
            root_name = name.substr(0, pos);
        }
        else {
            root_name = name;
        }

        // Add the root name if it is not present in the vector.
        if (std::find(root_names.begin(), root_names.end(), root_name) == root_names.end())
        {
            root_names.push_back(root_name);
        }
    }
    return root_names;
}

}  // namespace Kratos.
