//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//


#if !defined(KRATOS_IGA_MODELER_H_INCLUDED )
#define  KRATOS_IGA_MODELER_H_INCLUDED


// System includes

// External includes

// Project includes
#include "modeler.h"
#include "includes/properties.h"
#include "includes/define.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(KRATOS_CORE) IgaModeler
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(IgaModeler);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::GeometriesArrayType GeometriesArrayType;

    typedef typename Properties::Pointer PropertiesPointerType;

    typedef typename ModelPart::ElementsContainerType ElementsContainerType;
    typedef typename ModelPart::ConditionsContainerType ConditionsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IgaModeler(const Parameters ModelerParameters = Parameters())
        : Modeler(ModelerParameters)
        , mEchoLevel(
            ModelerParameters.Has("echo_level")
            ? ModelerParameters["echo_level"].GetInt()
            : 0)
    {
    }

    /// Destructor.
    virtual ~IgaModeler() = default;

    /// Creates the Modeler Pointer
    Modeler::Pointer Create(const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<IgaModeler>(ModelParameters);
    }

    ///@}
    ///@name Stages
    ///@{

    void GenerateModelPart(
        Model& rModel) const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "IgaModeler";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream & rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream & rOStream) const override
    {
    }

    ///@}

private:
    ///@name Iga functionalities
    ///@{

    SizeType mEchoLevel;

    ///@}
    ///@name Iga functionalities
    ///@{

    /// Creates integration point geometries and applies elements and conditions
    void CreateIntegrationDomain(
        ModelPart& rCadModelPart,
        ModelPart& rModelPart,
        const Parameters rParameters) const;

    void CreateIntegrationDomainPerUnit(
        ModelPart& rCadModelPart,
        ModelPart& rModelPart,
        const Parameters rParameters) const;

    /// Creates list of rQuadraturePointGeometryList
    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rQuadraturePointGeometryList,
        ModelPart& rModelPart,
        const Parameters rParameters) const;

    ///@}
    ///@name CAD functionalities
    ///@{

    /// Gets list of geometries from rModelPart
    void GetCadGeometryList(
        GeometriesArrayType& rGeometryList,
        ModelPart& rModelPart,
        const Parameters rParameters) const;

    ///@}
    ///@name Generate Elements and Conditions
    ///@{

    /// Creates elements from geometries
    void CreateElements(
        typename GeometriesArrayType::iterator rGeometriesBegin,
        typename GeometriesArrayType::iterator rGeometriesEnd,
        ModelPart& rDestinationModelPart,
        std::string& rElementName,
        SizeType& rIdCounter,
        PropertiesPointerType pProperties) const;

    /// Creates conditions from geometries
    void CreateConditions(
        typename GeometriesArrayType::iterator rGeometriesBegin,
        typename GeometriesArrayType::iterator rGeometriesEnd,
        ModelPart& rDestinationModelPart,
        std::string& rConditionName,
        SizeType& rIdCounter,
        PropertiesPointerType pProperties) const;

    ///@}
    ///@name Utility
    ///@{

    Parameters ReadParamatersFile(
        const std::string& rDataFileName) const;

    ///@}

}; // Class CadModeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    IgaModeler& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const IgaModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_IGA_MODELER_H_INCLUDED  defined