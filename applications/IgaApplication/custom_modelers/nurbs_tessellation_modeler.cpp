//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


// Project includes
#include "nurbs_tessellation_modeler.h"


namespace Kratos
{
    ///@name Stages
    ///@{

    void NurbsTessellationModeler::SetupModelPart()
    {
        KRATOS_ERROR_IF_NOT(mParameters.Has("cad_model_part_name")) << "Missing \"cad_model_part\" section" << std::endl;
        ModelPart& cad_model_part = mpModel->GetModelPart(mParameters["cad_model_part_name"].GetString());

        KRATOS_ERROR_IF_NOT(mParameters.Has("analysis_model_part_name")) << "Missing \"analysis_model_part_name\" section" << std::endl;
        ModelPart& analysis_model_part = mpModel->GetModelPart(mParameters["analysis_model_part_name"].GetString());

        const std::string& rDataFileName = mParameters.Has("physics_file_name")
            ? mParameters["physics_file_name"].GetString()
            : "physics.iga.json";

        KRATOS_WATCH(analysis_model_part)

        auto a = analysis_model_part.pGetGeometry(1);
        auto b = *a.get();
        KRATOS_WATCH(a->GetGeometryData());




        // KRATOS_WATCH(analysis_model_part.GetGeometry(1).Faces())
        // auto b_rep_model_vector = static_cast<GeometriesArrayType>(*a);


        // KRATOS_WATCH(b_rep_model_vector["breps"])




    }



    ///@}

}
