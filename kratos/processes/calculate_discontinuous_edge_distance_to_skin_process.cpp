//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand, Ruben Zorrilla, Franziska Wahl
//

// System includes

// External includes

// Project includes
#include "geometries/plane_3d.h"
#include "processes/calculate_discontinuous_edge_distance_to_skin_process.h"
#include "utilities/geometry_utilities.h"
#include "utilities/intersection_utilities.h"
#include "utilities/plane_approximation_utility.h"

namespace Kratos
{
	template<std::size_t TDim>
	CalculateDiscontinuousEdgeDistanceToSkinProcess<TDim>::CalculateDiscontinuousEdgeDistanceToSkinProcess(ModelPart& rVolumePart, ModelPart& rSkinPart)
		: CalculateDiscontinuousDistanceToSkinProcess<TDim>(rVolumePart, rSkinPart), mrSkinPart(rSkinPart), mrVolumePart(rVolumePart)
	{}

	template<std::size_t TDim>
	CalculateDiscontinuousEdgeDistanceToSkinProcess<TDim>::~CalculateDiscontinuousEdgeDistanceToSkinProcess()
	{}

	template<std::size_t TDim>
	void CalculateDiscontinuousEdgeDistanceToSkinProcess<TDim>::Initialize()
	{
		// Initialize the intersected objects process
		this->mFindIntersectedObjectsProcess.Initialize();

		// Initialize the elemental distances to the domain characteristic length
		const double initial_distance = this->CalculateCharacteristicLength();
		constexpr std::size_t num_nodes = TDim + 1;
		array_1d<double,num_nodes> init_dist_vect;
		for (unsigned int i_node = 0; i_node < num_nodes; ++i_node) {
			init_dist_vect[i_node] = initial_distance;
		}
		constexpr std::size_t num_edges = (TDim -1) * 3;
		array_1d<double,num_edges> init_edge_dist_vect;
		for (unsigned int i_node = 0; i_node < num_edges; ++i_node) {
			init_edge_dist_vect[i_node] = -1;
		}

		// Also initialize the embedded velocity of the fluid element and the TO_SPLIT flag.
		#pragma omp parallel for
		for (int k = 0; k< static_cast<int> (mrVolumePart.NumberOfElements()); ++k) {
			ModelPart::ElementsContainerType::iterator itElement = mrVolumePart.ElementsBegin() + k;
			itElement->Set(TO_SPLIT, false);
			itElement->SetValue(EMBEDDED_VELOCITY, ZeroVector(3));
			itElement->SetValue(ELEMENTAL_DISTANCES,init_dist_vect);
			itElement->SetValue(ELEMENTAL_EDGE_DISTANCES,init_edge_dist_vect);
		}
	}

	template<std::size_t TDim>
	void CalculateDiscontinuousEdgeDistanceToSkinProcess<TDim>::CalculateDistances(std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects)
	{
		const int num_elements = (this->mFindIntersectedObjectsProcess.GetModelPart1()).NumberOfElements();
		auto& r_elements = (this->mFindIntersectedObjectsProcess.GetModelPart1()).ElementsArray();

		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < num_elements; ++i) {
			this->CalculateEdgeDistances(*(r_elements[i]), rIntersectedObjects[i]);
			this->CalculateNodeDistances(*(r_elements[i]), rIntersectedObjects[i]);
		}
	}

	template<std::size_t TDim> // TODO remove method?
	bool CalculateDiscontinuousEdgeDistanceToSkinProcess<TDim>::CheckIfIncised(Element& rElement1)
	{
		// Get edge-based elemental distances of the element in order to calculate node-based ones
		const Vector r_edge_distances = rElement1.GetValue(ELEMENTAL_EDGE_DISTANCES);
		const uint8_t r_num_edges = rElement1.GetGeometry().EdgesNumber();

		// Check for edge intersections
		unsigned int num_cut_edges = 0;
		for (unsigned int i = 0; i < r_num_edges; i++) {
			if (r_edge_distances[i] >= 0){
				num_cut_edges++;
			}
		}
		return (num_cut_edges > 0);
	}

	template<std::size_t TDim> // TODO remove method?
	bool CalculateDiscontinuousEdgeDistanceToSkinProcess<TDim>::CheckIfIntersected(Element& rElement1)
	{
		// Get edge-based elemental distances of the element in order to calculate node-based ones
		const Vector r_edge_distances = rElement1.GetValue(ELEMENTAL_EDGE_DISTANCES);
		const uint8_t r_num_edges = rElement1.GetGeometry().EdgesNumber();

		// Check for edge intersections
		unsigned int num_cut_edges = 0;
		/*constexpr double epsilon = std::numeric_limits<double>::epsilon();
		unsigned int num_greater_epsilon = 0;*/
		for (unsigned int i = 0; i < r_num_edges; i++) {
			if (r_edge_distances[i] >= 0){
				num_cut_edges++;
				/*if (r_edge_distances[i] > epsilon){
					num_greater_epsilon++;
				}*/
			}
		}
		const bool is_intersection = (num_cut_edges < rElement1.GetGeometry().WorkingSpaceDimension()) ? false : true;
		//TODO: element can still not be completely intersected!?!
		return is_intersection /*&& num_greater_epsilon > 0*/;
	}

	/// Turn back information as a string.
	template<std::size_t TDim>
	std::string CalculateDiscontinuousEdgeDistanceToSkinProcess<TDim>::Info() const
	{
		return "CalculateDiscontinuousEdgeDistanceToSkinProcess";
	}

	template<std::size_t TDim>
	void CalculateDiscontinuousEdgeDistanceToSkinProcess<TDim>::CalculateEdgeDistances(
		Element& rElement1,
		PointerVector<GeometricalObject>& rIntersectedObjects)
	{
		if (rIntersectedObjects.empty()) {
			rElement1.Set(TO_SPLIT, false);
			return;
		}

		// This function assumes tetrahedra element and triangle intersected object as input at this moment
		constexpr int num_edges = (TDim - 1) * 3;
		constexpr double epsilon = std::numeric_limits<double>::epsilon();
		Vector& edge_distances = rElement1.GetValue(ELEMENTAL_EDGE_DISTANCES);
		//std::vector<double>&

		if(edge_distances.size() != num_edges){
			edge_distances.resize(num_edges, false);
		}

		// Compute the number of intersected edges
		std::vector<double> intersect_ratio_vector;
		const unsigned int num_cut_edges = ComputeEdgeIntersectionRatios(rElement1, rIntersectedObjects, intersect_ratio_vector);

		for (unsigned int i = 0; i < num_edges; i++) {
			edge_distances[i] = intersect_ratio_vector[i];
		}

		// Only complete intersection considered - 3 or more intersected edges for a tetrahedron, 2 or more for triangle
		const bool is_intersection = (num_cut_edges < rElement1.GetGeometry().WorkingSpaceDimension()) ? false : true;

		// TODO: TO_SPLIT only for completely intersected elements? epsilon??
		//       like this, we have almost the same definition of a split element as for CalculateDiscontinuousDistanceToSkinProcess??
		// Check if the element is split and set the TO_SPLIT flag accordingly
		unsigned int n_greater_epsilon = 0;
		for (unsigned int i = 0; i < num_edges; i++) {
			if (edge_distances[i] > epsilon){
				n_greater_epsilon++;
			}
		}
		rElement1.Set(TO_SPLIT, is_intersection && n_greater_epsilon > 0);
	}

	template<std::size_t TDim>
	unsigned int CalculateDiscontinuousEdgeDistanceToSkinProcess<TDim>::ComputeEdgeIntersectionRatios(
		Element& rElement1,
		const PointerVector<GeometricalObject>& rIntersectedObjects,
      	std::vector<double> &rIntersectionRatios)
	{
		auto &r_geometry = rElement1.GetGeometry();
		const auto r_edges_container = r_geometry.GenerateEdges();
		const std::size_t r_num_edges = r_geometry.EdgesNumber();

		// Initialize cut edges and points arrays
		unsigned int num_cut_edges = 0;
		rIntersectionRatios.clear();
		std::vector<unsigned int> cut_edges_vector = std::vector<unsigned int>(r_num_edges, 0);

		// Check wich edges are intersected
		for (std::size_t i_edge = 0; i_edge < r_num_edges; ++i_edge){
			array_1d<double,3> avg_pt = ZeroVector(3);
			std::vector<array_1d<double,3> > aux_pts;
			// Check against all candidates to count the number of current edge intersections
			for (const auto &r_int_obj : rIntersectedObjects){
				// Call the compute intersection method
				Point intersect_pt;
				const auto &r_int_obj_geom = r_int_obj.GetGeometry();
				const int intersect_id = ComputeEdgeIntersection(r_int_obj_geom, r_edges_container[i_edge][0], r_edges_container[i_edge][1], intersect_pt);

				// There is intersection
				if (intersect_id == 1){
					// Check if there is a close intersection (repeated intersection point)
					bool is_repeated = false;
					for (auto aux_pt : aux_pts){
						const double aux_dist = norm_2(intersect_pt - aux_pt);
						const double tol_edge = 1e-2*norm_2(r_edges_container[i_edge][0] - r_edges_container[i_edge][1]);
						if (aux_dist < tol_edge){
							is_repeated = true;
							break;
						}
					}

					// If the intersection pt. is not repeated, consider it
					if (!is_repeated){
						// Add the intersection pt. to the aux array pts.
						aux_pts.push_back(intersect_pt);
						// Increase the edge intersections counter
						cut_edges_vector[i_edge] += 1;
						// Save the intersection point for computing the average
						avg_pt += intersect_pt;
					}
				}
			}

			// No intersection if the edge is intersected a pair number of times
			// It is assumed that the skin enters and leaves the element
			// if (cut_edges_vector[i_edge] % 2 != 0){
			if (cut_edges_vector[i_edge] != 0){
				avg_pt /= cut_edges_vector[i_edge];
				// NEW calculate ratio
				double dist_avg_pt = norm_2(r_edges_container[i_edge][0] - avg_pt);
				double edge_length = norm_2(r_edges_container[i_edge][0] - r_edges_container[i_edge][1]);
				double avg_ratio =  dist_avg_pt / edge_length;
				rIntersectionRatios.push_back(avg_ratio);
				num_cut_edges++;
			} else {
				rIntersectionRatios.push_back(-1);
			}
		}
		return num_cut_edges;
	}

	template<std::size_t TDim>
	int CalculateDiscontinuousEdgeDistanceToSkinProcess<TDim>::ComputeEdgeIntersection(
		const Element::GeometryType& rIntObjGeometry,
		const Element::NodeType& rEdgePoint1,
		const Element::NodeType& rEdgePoint2,
		Point& rIntersectionPoint)
	{
		int intersection_flag = 0;
		const auto work_dim = rIntObjGeometry.WorkingSpaceDimension();
		if (work_dim == 2){
			intersection_flag = IntersectionUtilities::ComputeLineLineIntersection<Element::GeometryType>(
				rIntObjGeometry, rEdgePoint1.Coordinates(), rEdgePoint2.Coordinates(), rIntersectionPoint.Coordinates());
		} else if (work_dim == 3){
			intersection_flag = IntersectionUtilities::ComputeTriangleLineIntersection<Element::GeometryType>(
				rIntObjGeometry, rEdgePoint1.Coordinates(), rEdgePoint2.Coordinates(), rIntersectionPoint.Coordinates());
		} else {
			KRATOS_ERROR << "Working space dimension value equal to " << work_dim << ". Check your skin geometry implementation." << std::endl;
		}

		return intersection_flag;
	}

	template<std::size_t TDim>
	void CalculateDiscontinuousEdgeDistanceToSkinProcess<TDim>::CalculateNodeDistances(
		Element& rElement1,
		const PointerVector<GeometricalObject>& rIntersectedObjects)
	{
		if (rIntersectedObjects.empty()) {
			rElement1.Set(TO_SPLIT, false);
			return;
		}

		// Get edge-based elemental distances of the element in order to calculate node-based ones
		const Vector r_edge_distances = rElement1.GetValue(ELEMENTAL_EDGE_DISTANCES);
		//KRATOS_WATCH(r_edge_distances);
		constexpr int num_edges = (TDim - 1) * 3;

		// Check for valid (complete) intersections
		unsigned int num_cut_edges = 0;
		for (unsigned int i = 0; i < num_edges; i++) {
			if (r_edge_distances[i] >= 0){
				num_cut_edges++;
			}
		}
		const bool is_intersection = (num_cut_edges < rElement1.GetGeometry().WorkingSpaceDimension()) ? false : true;

		// Get reference to node-based elemental distances variable to change it after calculations
		// This function assumes tetrahedra element and triangle intersected object as input at this moment
		constexpr int num_nodes = TDim + 1;
		Vector& node_distances = rElement1.GetValue(ELEMENTAL_DISTANCES);

		if(node_distances.size() != num_nodes){
			node_distances.resize(num_nodes, false);
		}

		// Calculate node-based elemental distances for the element if there is a valid intersection
		if (is_intersection){
			// Calculate points from nodes of edges and length ratio of intersections
			std::vector<array_1d <double,3> > intsect_pts_vector;
			ComputeIntersectPtsFromRatios(rElement1, r_edge_distances, intsect_pts_vector);
			//KRATOS_WATCH(intsect_pts_vector);

			// If there are more than 3 intersected edges, compute the least squares plane approximation
			// by using the ComputePlaneApproximation utility. Otherwise, the distance is computed using
			// the plane defined by the 3 intersection points.
			auto &r_geometry = rElement1.GetGeometry();
			const bool do_plane_approx = (num_cut_edges == r_geometry.WorkingSpaceDimension()) ? false : true;

			if (do_plane_approx){
				// Call the plane optimization utility
				array_1d<double,3> base_pt, normal;
				ComputePlaneApproximation(rElement1, intsect_pts_vector, base_pt, normal);

				// Compute the distance to the approximation plane
				Plane3D approximation_plane(normal, Point{base_pt});
				for (int i = 0; i < num_nodes; i++) {
					node_distances[i] = approximation_plane.CalculateSignedDistance(r_geometry[i]);
				}
			} else {
				// Create a plane with the 3 intersection points (or 2 in 2D)
				Plane3D plane = this->SetIntersectionPlane(intsect_pts_vector);

				// Compute the distance to the intersection plane
				for (int i = 0; i < num_nodes; i++) {
					node_distances[i] = plane.CalculateSignedDistance(r_geometry[i]);
				}
			}

			// Correct the distance values orientation
			CorrectDistanceOrientation(r_geometry, rIntersectedObjects, node_distances);
			//KRATOS_WATCH(node_distances);
		}
	}

	template<std::size_t TDim>
	void CalculateDiscontinuousEdgeDistanceToSkinProcess<TDim>::ComputeIntersectPtsFromRatios(
		Element& rElement1,
		const Vector &rEdgeDistances,
      	std::vector<array_1d <double,3> > &rIntersectionPointsArray)
	{
		auto &r_geometry = rElement1.GetGeometry();
		const auto r_edges_container = r_geometry.GenerateEdges();
		const std::size_t r_num_edges = r_geometry.EdgesNumber();

		// Initialize point vector and points arrays
		array_1d<double,3> avg_pt = ZeroVector(3);
		rIntersectionPointsArray.clear();

		for (std::size_t i_edge = 0; i_edge < r_num_edges; ++i_edge){
			if (rEdgeDistances[i_edge] >= 0){
				avg_pt = r_edges_container[i_edge][0] + rEdgeDistances[i_edge] * (r_edges_container[i_edge][1] - r_edges_container[i_edge][0]);
				rIntersectionPointsArray.push_back(avg_pt);
			}
		}
	}

	template<std::size_t TDim>
	void CalculateDiscontinuousEdgeDistanceToSkinProcess<TDim>::ComputeIntersectionNormal(
		Element::GeometryType& rGeometry,
		const Vector& rElementalDistances,
		array_1d<double,3>& rNormal)
	{
		double volume;
		array_1d<double,TDim+1> N;
		BoundedMatrix<double,TDim+1,TDim> DN_DX;
		GeometryUtils::CalculateGeometryData(rGeometry, DN_DX, N, volume);

		rNormal = ZeroVector(3);
		for (std::size_t comp = 0; comp < TDim; ++comp){
			for (std::size_t i_node = 0; i_node < rGeometry.PointsNumber(); ++i_node){
				rNormal(comp) += DN_DX(i_node,comp)*rElementalDistances[i_node];
			}
		}
		rNormal /= norm_2(rNormal);
	}

	template<std::size_t TDim>
	void CalculateDiscontinuousEdgeDistanceToSkinProcess<TDim>::ComputePlaneApproximation(
		const Element& rElement1,
		const std::vector< array_1d<double,3> >& rPointsCoord,
		array_1d<double,3>& rPlaneBasePointCoords,
		array_1d<double,3>& rPlaneNormal)
	{
		const auto work_dim = rElement1.GetGeometry().WorkingSpaceDimension();
		if (work_dim == 2){
			PlaneApproximationUtility<2>::ComputePlaneApproximation(rPointsCoord, rPlaneBasePointCoords, rPlaneNormal);
		} else if (work_dim == 3){
			PlaneApproximationUtility<3>::ComputePlaneApproximation(rPointsCoord, rPlaneBasePointCoords, rPlaneNormal);
		} else {
			KRATOS_ERROR << "Working space dimension value equal to " << work_dim << ". Check your skin geometry implementation." << std::endl;
		}
	}

	template<std::size_t TDim>
	void CalculateDiscontinuousEdgeDistanceToSkinProcess<TDim>::CorrectDistanceOrientation(
		Element::GeometryType& rGeometry,
		const PointerVector<GeometricalObject>& rIntersectedObjects,
		Vector& rElementalDistances)
	{
		// Check the obtained intersection orientation (normal as distance gradient)
		array_1d<double,3> distance_normal;
		ComputeIntersectionNormal(rGeometry, rElementalDistances, distance_normal);

		// Vote the intersection orientation using the intersecting entities normals
		unsigned int n_pos = 0;
		unsigned int n_neg = 0;

		for (const auto &r_int_obj : rIntersectedObjects){
			const auto &r_int_obj_geom = r_int_obj.GetGeometry();

			array_1d<double, 3> r_int_obj_normal;
			ComputeIntersectionNormalFromGeometry(r_int_obj_geom, r_int_obj_normal);
			r_int_obj_normal /= norm_2(r_int_obj_normal);

			if (inner_prod(r_int_obj_normal, distance_normal) < 0.0){
				n_neg++;
			} else {
				n_pos++;
			}
		}

		// Negative votes win. Switch the distance values
		if (n_neg > n_pos){
			for (std::size_t i_node = 0; i_node < TDim + 1; ++i_node){
				rElementalDistances[i_node] *= -1.0;
			}
		}
	}

	template<>
	void inline CalculateDiscontinuousEdgeDistanceToSkinProcess<2>::ComputeIntersectionNormalFromGeometry(
		const Element::GeometryType &rGeometry,
		array_1d<double,3> &rIntObjNormal)
	{
		rIntObjNormal[0] = rGeometry[0].Y() - rGeometry[1].Y();
		rIntObjNormal[1] = rGeometry[1].X() - rGeometry[0].X();
		rIntObjNormal[2] = 0.0;
	}

	template<>
	void inline CalculateDiscontinuousEdgeDistanceToSkinProcess<3>::ComputeIntersectionNormalFromGeometry(
		const Element::GeometryType &rGeometry,
		array_1d<double,3> &rIntObjNormal)
	{
		MathUtils<double>::CrossProduct(rIntObjNormal, rGeometry[1]-rGeometry[0], rGeometry[2]-rGeometry[0]);
	}

	template class Kratos::CalculateDiscontinuousEdgeDistanceToSkinProcess<2>;
	template class Kratos::CalculateDiscontinuousEdgeDistanceToSkinProcess<3>;

}  // namespace Kratos.
