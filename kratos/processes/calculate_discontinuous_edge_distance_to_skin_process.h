//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Ruben Zorrilla
//

#if !defined(KRATOS_CALCULATE_DISCONTINUOUS_EDGE_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED )
#define  KRATOS_CALCULATE_DISCONTINUOUS_EDGE_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "geometries/plane_3d.h"
#include "includes/checks.h"
#include "processes/process.h"
#include "processes/find_intersected_geometrical_objects_process.h"
#include "processes/calculate_discontinuous_distance_to_skin_process.h"

namespace Kratos
{
///@addtogroup Kratos Core
///@{

///@name Kratos Classes
///@{

/// ((This only calculates the distance. Calculating the inside outside should be done by a derived class of this.))
/** This process takes a volume model part (with tetrahedra mesh) and a skin model part (with triangle mesh) and
     and calcualtes the distance to the skin for all the elements and nodes of the volume model part.
*/
template<std::size_t TDim = 3>
class KRATOS_API(KRATOS_CORE) CalculateDiscontinuousEdgeDistanceToSkinProcess : public CalculateDiscontinuousDistanceToSkinProcess<TDim>
{

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CalculateDiscontinuousEdgeDistanceToSkinProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateDiscontinuousEdgeDistanceToSkinProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor to be used.
    CalculateDiscontinuousEdgeDistanceToSkinProcess(
        ModelPart& rVolumePart,
        ModelPart& rSkinPart);

    /// Destructor.
    ~CalculateDiscontinuousEdgeDistanceToSkinProcess() override;

    ///@}
    ///@name Deleted
    ///@{

    /// Default constructor.
    CalculateDiscontinuousEdgeDistanceToSkinProcess() = delete;

    /// Copy constructor.
    CalculateDiscontinuousEdgeDistanceToSkinProcess(CalculateDiscontinuousEdgeDistanceToSkinProcess const& rOther) = delete;

    /// Assignment operator.
    CalculateDiscontinuousEdgeDistanceToSkinProcess& operator=(CalculateDiscontinuousEdgeDistanceToSkinProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Initializes discontinuous distance computation process
     * This method initializes the TO_SPLIT flag, the DISTANCE and
     * ELEMENTAL_DISTANCES variables as well as the EMBEDDED_VELOCITY
     *
     * NEW - initialize ELEMENTAL_EDGE_DISTANCES aswell
     */
    virtual void Initialize() override;

    /**
     * @brief Computes the elemental distance values
     * Given an intersecting objects vector, this method computes the elemental distance field
     * @param rIntersectedObjects array containing pointers to the intersecting geometries
     *
     * NEW ...
     */
    virtual void CalculateDistances(std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects) override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    ///@}
protected:
    ///@name Protected Operations
    ///@{

    ///@}
private:
    ///@name Member Variables
    ///@{

    ModelPart& mrSkinPart;
    ModelPart& mrVolumePart;
    //std::size_t mNumEdges = (TDim -1) * 3; //TODO: declare here?

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Computes the edge distances of one element (EMBEDDED_EDGE_DISTANCES)
     * This method computes ratios for intersections along the edges of a given element
     * @param rElement1 reference to the element of interest
     * @param rIntersectedObjects reference to the array containing the element of interest intersecting geometries
     *
     * NEW
     */
    void CalculateEdgeDistances(
        Element& rElement1,
        PointerVector<GeometricalObject>& rIntersectedObjects);

    /**
     * @brief Computes the edges intersections in one element
     * Provided a list of elemental intersecting geometries, this
     * method computes the edge intersections for a given element
     * @param rElement1 reference to the element of interest
     * @param rIntersectedObjects reference to the array containing the element of interest intersecting geometries
     * @param rCutEdgesVector array that classifies the edges depending on their cut / uncut status
     * @param rIntersectionRatiosArray array containing the edges intersection ratios
     * @return unsigned int number of cut edges
     *
     * NEW
     */
    unsigned int ComputeEdgeIntersectionRatios(
        Element& rElement1,
        const PointerVector<GeometricalObject>& rIntersectedObjects,
        std::vector<double> &rIntersectionRatios);

    /**
     * @brief Computes the intersection of a single edge
     * This method computes the intersection of a given edge with the candidate
     * intersecting geometry. This operation is performed accordingly to the working
     * space dimension using the intersection utilities implemented in intersection_utilities.h
     * @param rIntObjGeometry candidate intersecting geometry
     * @param rEdgePoint1 edge origin point
     * @param rEdgePoint2 edge end point
     * @param rIntersectionPoint intersection point
     * @return int type of intersection id (see intersection_utilities.h)
     */
    int ComputeEdgeIntersection(
        const Element::GeometryType& rIntObjGeometry,
        const Element::NodeType& rEdgePoint1,
        const Element::NodeType& rEdgePoint2,
        Point& rIntersectionPoint);

    /**
     * @brief Computes the discontinuous distance for nodes of one element (ELEMENTAL_DISTANCES)
     * This method computes the discontinuous distance field for a given element
     * @param rElement1 reference to the element of interest
     * @param rIntersectedObjects reference to the array containing the element of interest intersecting geometries
     */
    void CalculateNodeDistances(
        Element& rElement1,
        const PointerVector<GeometricalObject>& rIntersectedObjects);

    /**
     * @brief
     * @param rElement1 reference to the element of interest
     * @param rEdgeDistances
     * @param rIntersectionPointsArray
     *
     * NEW
     */
    void ComputeIntersectPtsFromRatios(
        Element& rElement1,
        const Vector &rEdgeDistances,
        std::vector<array_1d <double,3> > &rIntersectionPointsArray);

    /** TODO: make base class method protected?
     * @brief Computes the element intersection unit normal
     * This method computes the element intersection unit normal vector using the distance function gradient.
     * @param rGeometry reference to the geometry of the element of interest
     * @param rElementalDistances array containing the ELEMENTAL_DISTANCES values
     * @param rNormal obtained unit normal vector
     */
    void ComputeIntersectionNormal(
        Element::GeometryType& rGeometry,
        const Vector& rElementalDistances,
        array_1d<double,3> &rNormal);

    /** TODO: make base class method protected?
     * @brief Computes the intersection plane approximation
     * For complex intersection patterns, this method takes a list containing
     * all the intersecting points and computes the plane that minimizes the
     * distance from all these points in a least squares sense. The approximated
     * plane is defined in terms of an origin point and its normal vector.
     * @param rElement1 reference to the element of interest
     * @param rPointsCoord list containing the coordinates of al the intersecting points
     * @param rPlaneBasePointCoords base point defining the approximated plane
     * @param rPlaneNormal normal vector defining the approximated plane
     */
    void ComputePlaneApproximation(
        const Element& rElement1,
        const std::vector< array_1d<double,3> >& rPointsCoord,
        array_1d<double,3>& rPlaneBasePointCoords,
        array_1d<double,3>& rPlaneNormal);

    /** TODO: make base class method protected?
     * @brief Checks (and corrects if needed) the intersection normal orientation
     * This method checks the orientation of the previously computed intersection normal.
     * To do that, the normal vector to each one of the intersecting geometries is
     * computed and its directo is compared against the current one. If the negative
     * votes win, the current normal vector orientation is switched.
     * @param rGeometry element of interest geometry
     * @param rIntersectedObjects reference to the array containing the element of interest intersecting geometries
     * @param rElementalDistances array containing the ELEMENTAL_DISTANCES values
     */
    void CorrectDistanceOrientation(
        Element::GeometryType& rGeometry,
        const PointerVector<GeometricalObject>& rIntersectedObjects,
        Vector& rElementalDistances);

    /** TODO: make base class method protected?
     * @brief Computes the normal vector to an intersecting object geometry
     * This method computes the normal vector to an intersecting object geometry.
     * @param rGeometry reference to the geometry of the intersecting object
     * @param rIntObjNormal reference to the intersecting object normal vector
     */
    void inline ComputeIntersectionNormalFromGeometry(
        const Element::GeometryType &rGeometry,
        array_1d<double,3> &rIntObjNormal);

    /** TODO: make base class method protected?
     * @brief Computes the value of any embedded variable
     * For a given array variable in the skin mesh, this method calculates the value
     * of such variable in the embedded mesh. This is done in each element of the volume
     * mesh by computing the average value of all the edges intersections. This value
     * is averaged again according to the number of intersected edges.
     * @tparam TVarType variable type
     * @param rVariable origin variable in the skin mesh
     * @param rEmbeddedVariable elemental variable in the volume mesh to be computed
     */
    template<class TVarType>
	void CalculateEmbeddedVariableFromSkinSpecialization(
		const Variable<TVarType> &rVariable,
		const Variable<TVarType> &rEmbeddedVariable)
	{
		const auto &r_int_obj_vect= this->GetIntersections();
		const int n_elems = mrVolumePart.NumberOfElements();

		// Check requested variables
		KRATOS_ERROR_IF(rEmbeddedVariable.Key() == 0)
			<< rEmbeddedVariable << " key is 0. Check that the variable is correctly registered." << std::endl;

		KRATOS_ERROR_IF((mrSkinPart.NodesBegin())->SolutionStepsDataHas(rVariable) == false)
			<< "Skin model part solution step data missing variable: " << rVariable << std::endl;

		// Initialize embedded variable value
		#pragma omp parallel for
		for (int i_elem = 0; i_elem < n_elems; ++i_elem) {
			auto it_elem = mrVolumePart.ElementsBegin() + i_elem;
			it_elem->SetValue(rEmbeddedVariable, rEmbeddedVariable.Zero());
		}

		// Compute the embedded variable value for each element
		#pragma omp parallel for schedule(dynamic)
		for (int i_elem = 0; i_elem < n_elems; ++i_elem) {
			// Check if the current element has intersecting entities
			if (r_int_obj_vect[i_elem].size() != 0) {
				// Initialize the element values
				unsigned int n_int_edges = 0;
				auto it_elem = mrVolumePart.ElementsBegin() + i_elem;
				auto &r_geom = it_elem->GetGeometry();
				const auto edges = r_geom.GenerateEdges();

				// Loop the element of interest edges
				for (unsigned int i_edge = 0; i_edge < r_geom.EdgesNumber(); ++i_edge) {
					// Initialize edge values
					unsigned int n_int_obj = 0;
					TVarType i_edge_val = rEmbeddedVariable.Zero();

					// Check the edge intersection against all the candidates
					for (auto &r_int_obj : r_int_obj_vect[i_elem]) {
						Point intersection_point;
						const int is_intersected = this->ComputeEdgeIntersection(
							r_int_obj.GetGeometry(),
							edges[i_edge][0],
							edges[i_edge][1],
							intersection_point);

						// Compute the variable value in the intersection point
                        if (is_intersected == 1) {
							n_int_obj++;
							array_1d<double,3> local_coords;
                            r_int_obj.GetGeometry().PointLocalCoordinates(local_coords, intersection_point);
                            Vector int_obj_N;
                            r_int_obj.GetGeometry().ShapeFunctionsValues(int_obj_N, local_coords);
                            for (unsigned int i_node = 0; i_node < r_int_obj.GetGeometry().PointsNumber(); ++i_node) {
                                i_edge_val += r_int_obj.GetGeometry()[i_node].FastGetSolutionStepValue(rVariable) * int_obj_N[i_node];
                            }
                        }
					}

					// Check if the edge is intersected
					if (n_int_obj != 0) {
						// Update the element intersected edges counter
						n_int_edges++;
						// Add the average edge value (there might exist cases in where
						// more than one geometry intersects the edge of interest).
						it_elem->GetValue(rEmbeddedVariable) += i_edge_val / n_int_obj;
					}
				}

				// Average between all the intersected edges
				if (n_int_edges != 0) {
					it_elem->GetValue(rEmbeddedVariable) /= n_int_edges;
				}
			}
		}
	};

    ///@}

}; // Class CalculateDiscontinuousEdgeDistanceToSkinProcess

///@}

///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    CalculateDiscontinuousEdgeDistanceToSkinProcess<>& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const CalculateDiscontinuousEdgeDistanceToSkinProcess<>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CALCULATE_DISCONTINUOUS_EDGE_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED  defined
