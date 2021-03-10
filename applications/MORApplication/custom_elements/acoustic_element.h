// KRATOS
//
//  License:		BSD License
//					license: ../../license.txt
//
//  Main authors:   Elias Schwaiger
//                  Quirin Aumann
//
//


#if !defined(KRATOS_ACOUSTIC_ELEMENT_H_INCLUDED )
#define  KRATOS_ACOUSTIC_ELEMENT_H_INCLUDED

// System includes


// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "includes/serializer.h"

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
 * @class AcousticElement
 * @ingroup MORApplication
 * @brief Acoustic element
 * @details
 * @author
 * @author
 */

class AcousticElement
    : public Element
{
public:
    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// Type for shape function values container
    typedef Kratos::Vector ShapeFunctionsType;

    /// Type for a matrix containing the shape function gradients
    typedef Kratos::Matrix ShapeFunctionDerivativesType;

    /// Type for an array of shape function gradient matrices
    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    /// The base element type
    typedef Element BaseType;

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    /// Counted pointer of AcousticElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AcousticElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AcousticElement(IndexType NewId, GeometryType::Pointer pGeometry);
    AcousticElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    // Copy constructor
    AcousticElement(AcousticElement const& rOther)
        :BaseType(rOther)
    {};

    /// Destructor.
    ~AcousticElement() override;

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param pGeom The pointer to the geometry of the element
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param ThisNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief It creates a new element pointer and clones the previous element data
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone (
        IndexType NewId,
        NodesArrayType const& rThisNodes
        ) const override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The current process info instance
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Acoustic Element #" << Id() << "\n";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Acoustic Element #" << Id() << "\n";
    }


    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;

    // void Initialize(const ProcessInfo& rCurrentProcessInfo) override;


protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    IntegrationMethod mThisIntegrationMethod;

    ///@}
    ///@name Protected Operators
    ///@{

    AcousticElement() : Element()
    {
    }

    ///@}
    ///@name Protected Operations
    ///@{

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;


    void SetIntegrationMethod(const IntegrationMethod& ThisIntegrationMethod)
    {
         mThisIntegrationMethod = ThisIntegrationMethod;
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

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  AcSizeType
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    /// Assignment operator.
    //AcousticElement& operator=(const AcousticElement& rOther);
    /// Copy constructor.
    //AcousticElement(const AcousticElement& rOther);
    ///@}

}; // Class AcousticElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_ACOUSTIC_ELEMENT_H_INCLUDED  defined
