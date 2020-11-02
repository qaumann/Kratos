//#include "includes/kratos_flags.h"

#include "custom_elements/embedded_fluid_element_discontinuous_edge.h"
#include "custom_elements/qs_vms.h"
#include "custom_elements/symbolic_navier_stokes.h"

//#include "utilities/element_size_calculator.h"
//#include "custom_utilities/embedded_discontinuous_data.h"
#include "custom_utilities/symbolic_navier_stokes_data.h"
#include "custom_utilities/time_integrated_qsvms_data.h"

#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions.h"

namespace Kratos {

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TBaseElement >
EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::EmbeddedFluidElementDiscontinuousEdge(IndexType NewId)
    : EmbeddedFluidElementDiscontinuous<TBaseElement>(NewId)
{}

template< class TBaseElement >
EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::EmbeddedFluidElementDiscontinuousEdge(IndexType NewId, const NodesArrayType& ThisNodes)
    : EmbeddedFluidElementDiscontinuous<TBaseElement>(NewId,ThisNodes)
{}


template< class TBaseElement >
EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::EmbeddedFluidElementDiscontinuousEdge(IndexType NewId, Geometry<NodeType>::Pointer pGeometry)
    : EmbeddedFluidElementDiscontinuous<TBaseElement>(NewId,pGeometry)
{}

template< class TBaseElement >
EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::EmbeddedFluidElementDiscontinuousEdge(IndexType NewId, Geometry<NodeType>::Pointer pGeometry, Properties::Pointer pProperties)
    : EmbeddedFluidElementDiscontinuous<TBaseElement>(NewId,pGeometry,pProperties)
{}


template< class TBaseElement >
EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::~EmbeddedFluidElementDiscontinuousEdge()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations TODO


///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <class TBaseElement>
int EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    int out = EmbeddedDiscontinuousElementData::Check(*this, rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Something is wrong with the elemental data of Element "
        << this->Info() << std::endl;

    return TBaseElement::Check(rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <class TBaseElement>
std::string EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::Info() const
{
    std::stringstream buffer;
    buffer << "EmbeddedFluidElementDiscontinuousEdge #" << this->Id();
    return buffer.str();
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "EmbeddedFluidElementDiscontinuousEdge" << Dim << "D" << NumNodes << "N"
             << std::endl
             << "on top of ";
    TBaseElement::PrintInfo(rOStream);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions TODO
///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
// Operations


///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions TODO
///////////////////////////////////////////////////////////////////////////////////////////////////

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::CalculateDragForce(
    EmbeddedDiscontinuousElementData& rData,
    array_1d<double,3>& rDragForce) const
{
    // Initialize the embedded element data
    const unsigned int number_of_positive_gauss_points = rData.PositiveSideWeights.size();
    const unsigned int number_of_negative_gauss_points = rData.NegativeSideWeights.size();
    const size_t volume_gauss_points = number_of_positive_gauss_points + number_of_negative_gauss_points;

    if (rData.IsCut()){
        // Integrate positive interface side drag
        const unsigned int n_int_pos_gauss = rData.PositiveInterfaceWeights.size();
        for (unsigned int g = 0; g < n_int_pos_gauss; ++g) {
            // Update the Gauss pt. data and the constitutive law
            this->UpdateIntegrationPointData(
                rData,
                g + volume_gauss_points,
                rData.PositiveInterfaceWeights[g],
                row(rData.PositiveInterfaceN, g),
                rData.PositiveInterfaceDNDX[g]);

            // Get the interface Gauss pt. unit noromal
            const auto &aux_unit_normal = rData.PositiveInterfaceUnitNormals[g];

            // Compute Gauss pt. pressure
            const double p_gauss = inner_prod(rData.N, rData.Pressure);

            // Get the normal projection matrix in Voigt notation
            BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
            FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

            // Add the shear and pressure drag contributions
            const array_1d<double, Dim> shear_proj = rData.Weight * prod(voigt_normal_proj_matrix, rData.ShearStress);
            for (unsigned int i = 0; i < Dim ; ++i){
                rDragForce(i) -= shear_proj(i);
            }
            rDragForce += rData.Weight * p_gauss * aux_unit_normal;
        }

        // Integrate negative interface side drag
        const unsigned int n_int_neg_gauss = rData.NegativeInterfaceWeights.size();
        for (unsigned int g = 0; g < n_int_neg_gauss; ++g) {
            // Update the Gauss pt. data and the constitutive law
            this->UpdateIntegrationPointData(
                rData,
                g + volume_gauss_points + n_int_pos_gauss,
                rData.NegativeInterfaceWeights[g],
                row(rData.NegativeInterfaceN, g),
                rData.NegativeInterfaceDNDX[g]);

            // Get the interface Gauss pt. unit noromal
            const auto &aux_unit_normal = rData.NegativeInterfaceUnitNormals[g];

            // Compute Gauss pt. pressure
            const double p_gauss = inner_prod(rData.N, rData.Pressure);

            // Get the normal projection matrix in Voigt notation
            BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
            FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

            // Add the shear and pressure drag contributions
            const array_1d<double, Dim> shear_proj = rData.Weight * prod(voigt_normal_proj_matrix, rData.ShearStress);
            for (unsigned int i = 0; i < Dim ; ++i){
                rDragForce(i) -= shear_proj(i);
            }
            rDragForce += rData.Weight * p_gauss * aux_unit_normal;
        }
    }
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::CalculateDragForceCenter(
    EmbeddedDiscontinuousElementData& rData,
    array_1d<double,3>& rDragForceLocation) const
{
    const auto &r_geometry = this->GetGeometry();
    array_1d<double,3> tot_drag = ZeroVector(3);
    const unsigned int number_of_positive_gauss_points = rData.PositiveSideWeights.size();
    const unsigned int number_of_negative_gauss_points = rData.NegativeSideWeights.size();
    const size_t volume_gauss_points = number_of_positive_gauss_points + number_of_negative_gauss_points;

    if (rData.IsCut()){
        // Get the positive interface continuous shape functions
        // We use these ones to interpolate the position of the intersection Gauss pt.
        // Note that we take advantage of the fact that the positive and negative interface Gauss pt. coincide
        Vector pos_int_continuous_weights;
        Matrix pos_int_continuous_N;
        typename EmbeddedDiscontinuousElementData::ShapeFunctionsGradientsType pos_int_continuous_DN_DX;
        auto p_continuous_sh_func_calculator = EmbeddedDiscontinuousInternals::GetContinuousShapeFunctionCalculator<Dim, NumNodes>(*this, rData.ElementalDistances);
        p_continuous_sh_func_calculator->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
            pos_int_continuous_N,
            pos_int_continuous_DN_DX,
            pos_int_continuous_weights,
            GeometryData::GI_GAUSS_2);

        // Integrate positive interface side drag
        const unsigned int n_int_pos_gauss = rData.PositiveInterfaceWeights.size();
        for (unsigned int g = 0; g < n_int_pos_gauss; ++g) {
            // Obtain the Gauss pt. coordinates using the standard shape functions
            array_1d<double,3> g_coords = ZeroVector(3);
            const auto g_shape_functions = row(pos_int_continuous_N, g);
            for (unsigned int i_node = 0; i_node < NumNodes; ++i_node) {
                g_coords += g_shape_functions[i_node] * r_geometry[i_node].Coordinates();
            }

            // Update the Gauss pt. data and the constitutive law
            this->UpdateIntegrationPointData(
                rData,
                g + volume_gauss_points,
                rData.PositiveInterfaceWeights[g],
                row(rData.PositiveInterfaceN, g),
                rData.PositiveInterfaceDNDX[g]);

            // Get the interface Gauss pt. unit noromal
            const auto &aux_unit_normal = rData.PositiveInterfaceUnitNormals[g];

            // Compute Gauss pt. pressure
            const double p_gauss = inner_prod(rData.N, rData.Pressure);

            // Get the normal projection matrix in Voigt notation
            BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
            FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

            // Add the shear and pressure drag contributions
            const array_1d<double, 3> p_proj = rData.Weight * p_gauss * aux_unit_normal;
            const array_1d<double, Dim> shear_proj = rData.Weight * prod(voigt_normal_proj_matrix, rData.ShearStress);
            for (unsigned int i = 0; i < Dim ; ++i){
                tot_drag(i) -= shear_proj(i);
                rDragForceLocation(i) += g_coords(i) * p_proj(i);
                rDragForceLocation(i) -= g_coords(i) * shear_proj(i);
            }
            tot_drag += p_proj;
        }

        // Integrate negative interface side drag
        const unsigned int n_int_neg_gauss = rData.NegativeInterfaceWeights.size();
        for (unsigned int g = 0; g < n_int_neg_gauss; ++g) {
            // Obtain the Gauss pt. coordinates using the standard shape functions
            array_1d<double,3> g_coords = ZeroVector(3);
            const auto g_shape_functions = row(pos_int_continuous_N, g);
            for (unsigned int i_node = 0; i_node < NumNodes; ++i_node) {
                g_coords += g_shape_functions[i_node] * r_geometry[i_node].Coordinates();
            }

            // Update the Gauss pt. data and the constitutive law
            this->UpdateIntegrationPointData(
                rData,
                g + volume_gauss_points + n_int_pos_gauss,
                rData.NegativeInterfaceWeights[g],
                row(rData.NegativeInterfaceN, g),
                rData.NegativeInterfaceDNDX[g]);

            // Get the interface Gauss pt. unit noromal
            const auto &aux_unit_normal = rData.NegativeInterfaceUnitNormals[g];

            // Compute Gauss pt. pressure
            const double p_gauss = inner_prod(rData.N, rData.Pressure);

            // Get the normal projection matrix in Voigt notation
            BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
            FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

            // Add the shear and pressure drag contributions
            const array_1d<double, 3> p_proj = rData.Weight * p_gauss * aux_unit_normal;
            const array_1d<double, Dim> shear_proj = rData.Weight * prod(voigt_normal_proj_matrix, rData.ShearStress);
            for (unsigned int i = 0; i < Dim ; ++i){
                tot_drag(i) -= shear_proj(i);
                rDragForceLocation(i) += g_coords(i) * p_proj(i);
                rDragForceLocation(i) -= g_coords(i) * shear_proj(i);
            }
            tot_drag += p_proj;
        }

        // Divide the obtained result by the total drag
        rDragForceLocation(0) /= tot_drag(0);
        rDragForceLocation(1) /= tot_drag(1);
        if (Dim == 3) {
            rDragForceLocation(2) /= tot_drag(2);
        }
    }
}

// serializer

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TBaseElement);
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuousEdge<TBaseElement>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TBaseElement);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Helper functions for template specialization
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace EmbeddedDiscontinuousEdgeInternals {

template <>
ModifiedShapeFunctions::Pointer GetShapeFunctionCalculator<2, 3>(const Element& rElement, const Vector& rElementalDistances)
{
    return ModifiedShapeFunctions::Pointer(new Triangle2D3AusasModifiedShapeFunctions(rElement.pGetGeometry(), rElementalDistances));
}

template <>
ModifiedShapeFunctions::Pointer GetShapeFunctionCalculator<3, 4>(const Element& rElement, const Vector& rElementalDistances)
{
    return ModifiedShapeFunctions::Pointer(new Tetrahedra3D4AusasModifiedShapeFunctions(rElement.pGetGeometry(), rElementalDistances));
}

template <>
ModifiedShapeFunctions::Pointer GetContinuousShapeFunctionCalculator<2, 3>(
    const Element& rElement,
    const Vector& rElementalDistances)
{
    return ModifiedShapeFunctions::Pointer(new Triangle2D3ModifiedShapeFunctions(rElement.pGetGeometry(), rElementalDistances));
}

template <>
ModifiedShapeFunctions::Pointer GetContinuousShapeFunctionCalculator<3, 4>(
    const Element& rElement,
    const Vector& rElementalDistances)
{
    return ModifiedShapeFunctions::Pointer(new Tetrahedra3D4ModifiedShapeFunctions(rElement.pGetGeometry(), rElementalDistances));
}

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class EmbeddedFluidElementDiscontinuousEdge< QSVMS< TimeIntegratedQSVMSData<2,3> > >;
template class EmbeddedFluidElementDiscontinuousEdge< QSVMS< TimeIntegratedQSVMSData<3,4> > >;

template class EmbeddedFluidElementDiscontinuousEdge< SymbolicNavierStokes< SymbolicNavierStokesData<2,3> > >;
template class EmbeddedFluidElementDiscontinuousEdge< SymbolicNavierStokes< SymbolicNavierStokesData<3,4> > >;

///////////////////////////////////////////////////////////////////////////////////////////////////

}