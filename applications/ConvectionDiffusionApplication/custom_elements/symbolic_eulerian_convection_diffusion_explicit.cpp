//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Riccardo Tosi
//

#include "symbolic_eulerian_convection_diffusion_explicit.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

SymbolicEulerianConvectionDiffusionExplicit::SymbolicEulerianConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

SymbolicEulerianConvectionDiffusionExplicit::SymbolicEulerianConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

SymbolicEulerianConvectionDiffusionExplicit::~SymbolicEulerianConvectionDiffusionExplicit() {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

Element::Pointer SymbolicEulerianConvectionDiffusionExplicit::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicEulerianConvectionDiffusionExplicit>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


Element::Pointer SymbolicEulerianConvectionDiffusionExplicit::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicEulerianConvectionDiffusionExplicit>(NewId, pGeom, pProperties);
}

void SymbolicEulerianConvectionDiffusionExplicit::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const unsigned int LocalSize = r_geometry.size();
    const unsigned int Dimension = r_geometry.WorkingSpaceDimension();

    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    // Element variables
    ElementVariables rVariables;
    this->InitializeEulerianElement(rVariables,rCurrentProcessInfo);

    // Reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Define local variables
    Element::GeometryType::JacobiansType J0;
    Matrix DN(LocalSize,Dimension);
    Matrix InvJ0(Dimension,Dimension);
    Vector temp(LocalSize);
    double DetJ0;
    MatrixType lhs;
    VectorType rhs;

    // Compute Jacobian
    r_geometry.Jacobian(J0,this->GetIntegrationMethod());

    // Iterate over integration points to evaluate local contribution
    for (unsigned int g = 0; g < integration_points.size(); g++) {

        const double k = rVariables.diffusivity;
        const Vector f = rVariables.forcing;
        const Vector phi = rVariables.unknown;

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(J0[g],InvJ0,DetJ0);
        // Calculate the cartesian derivatives on integration point "g"
        noalias(DN) = prod(DN_De[g],InvJ0);
        // Caluclate N on the gauss point "g"
        auto N = row(N_gausspoint,g);
        // Compute weight
        const double IntToReferenceWeight = integration_points[g].Weight() * DetJ0;

        if (Dimension == 2){

            const double clhs0 =             k*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1));
const double clhs1 =             k*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1));
const double clhs2 =             k*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1));
            lhs(0,0)=k*(pow(DN(0,0), 2) + pow(DN(0,1), 2));
            lhs(0,1)=clhs0;
            lhs(0,2)=clhs1;
            lhs(1,0)=clhs0;
            lhs(1,1)=k*(pow(DN(1,0), 2) + pow(DN(1,1), 2));
            lhs(1,2)=clhs2;
            lhs(2,0)=clhs1;
            lhs(2,1)=clhs2;
            lhs(2,2)=k*(pow(DN(2,0), 2) + pow(DN(2,1), 2));


            const double crhs0 =             N[0]*f[0] + N[1]*f[1] + N[2]*f[2];
const double crhs1 =             DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2];
const double crhs2 =             DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2];
            rhs[0]=k*(-DN(0,0)*crhs1 - DN(0,1)*crhs2 + N[0]*crhs0);
            rhs[1]=k*(-DN(1,0)*crhs1 - DN(1,1)*crhs2 + N[1]*crhs0);
            rhs[2]=k*(-DN(2,0)*crhs1 - DN(2,1)*crhs2 + N[2]*crhs0);

        }
        else if (Dimension == 3){

            const double clhs0 =             k*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1) + DN(0,2)*DN(1,2));
const double clhs1 =             k*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1) + DN(0,2)*DN(2,2));
const double clhs2 =             k*(DN(0,0)*DN(3,0) + DN(0,1)*DN(3,1) + DN(0,2)*DN(3,2));
const double clhs3 =             k*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1) + DN(1,2)*DN(2,2));
const double clhs4 =             k*(DN(1,0)*DN(3,0) + DN(1,1)*DN(3,1) + DN(1,2)*DN(3,2));
const double clhs5 =             k*(DN(2,0)*DN(3,0) + DN(2,1)*DN(3,1) + DN(2,2)*DN(3,2));
            lhs(0,0)=k*(pow(DN(0,0), 2) + pow(DN(0,1), 2) + pow(DN(0,2), 2));
            lhs(0,1)=clhs0;
            lhs(0,2)=clhs1;
            lhs(0,3)=clhs2;
            lhs(1,0)=clhs0;
            lhs(1,1)=k*(pow(DN(1,0), 2) + pow(DN(1,1), 2) + pow(DN(1,2), 2));
            lhs(1,2)=clhs3;
            lhs(1,3)=clhs4;
            lhs(2,0)=clhs1;
            lhs(2,1)=clhs3;
            lhs(2,2)=k*(pow(DN(2,0), 2) + pow(DN(2,1), 2) + pow(DN(2,2), 2));
            lhs(2,3)=clhs5;
            lhs(3,0)=clhs2;
            lhs(3,1)=clhs4;
            lhs(3,2)=clhs5;
            lhs(3,3)=k*(pow(DN(3,0), 2) + pow(DN(3,1), 2) + pow(DN(3,2), 2));


            const double crhs0 =             N[0]*f[0] + N[1]*f[1] + N[2]*f[2] + N[3]*f[3];
const double crhs1 =             DN(0,0)*phi[0] + DN(1,0)*phi[1] + DN(2,0)*phi[2] + DN(3,0)*phi[3];
const double crhs2 =             DN(0,1)*phi[0] + DN(1,1)*phi[1] + DN(2,1)*phi[2] + DN(3,1)*phi[3];
const double crhs3 =             DN(0,2)*phi[0] + DN(1,2)*phi[1] + DN(2,2)*phi[2] + DN(3,2)*phi[3];
            rhs[0]=k*(-DN(0,0)*crhs1 - DN(0,1)*crhs2 - DN(0,2)*crhs3 + N[0]*crhs0);
            rhs[1]=k*(-DN(1,0)*crhs1 - DN(1,1)*crhs2 - DN(1,2)*crhs3 + N[1]*crhs0);
            rhs[2]=k*(-DN(2,0)*crhs1 - DN(2,1)*crhs2 - DN(2,2)*crhs3 + N[2]*crhs0);
            rhs[3]=k*(-DN(3,0)*crhs1 - DN(3,1)*crhs2 - DN(3,2)*crhs3 + N[3]*crhs0);

        }

        noalias(rLeftHandSideMatrix) += lhs * IntToReferenceWeight;
        noalias(rRightHandSideVector) += rhs * IntToReferenceWeight;

    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

int SymbolicEulerianConvectionDiffusionExplicit::Check(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;
    int out = Element::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    return 0;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public I/O


///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected operations

void SymbolicEulerianConvectionDiffusionExplicit::InitializeEulerianElement(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const ProcessInfo& r_process_info = rCurrentProcessInfo;
        ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
        auto& r_settings = *p_settings;

        const auto& r_geometry = GetGeometry();
        const unsigned int LocalSize = r_geometry.size();

        rVariables.diffusivity = r_settings.GetDiffusionVariable();

        for(unsigned int node_element = 0; node_element<LocalSize; node_element++)
    {
        rVariables.unknown[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetUnknownVariable());
        rVariables.forcing[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_settings.GetVolumeSourceVariable());
    }

        KRATOS_CATCH( "" )
    }

Element::IntegrationMethod SymbolicEulerianConvectionDiffusionExplicit::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private serialization

}