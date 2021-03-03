//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


// System includes


// External includes


// Project includes
#include "includes/define.h"

#include "mor_application.h"
#include "mor_application_variables.h"

#include "geometries/line_2d_2.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/point_3d.h"
#include "geometries/triangle_3d_3.h"


namespace Kratos {

  typedef Node<3> NodeType;

KratosMORApplication::KratosMORApplication():
    KratosApplication("MORApplication"),
      mAcousticElement2D2N(0, Element::GeometryType::Pointer(new Line2D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
      mAcousticElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mAcousticElement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mAcousticElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mAcousticElement3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mAcousticPMLElement2D2N(0, Element::GeometryType::Pointer(new Line2D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
      mAcousticPMLElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mAcousticPMLElement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mAcousticPMLElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mAcousticPMLElement3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mPorousElement2D2N(0, Element::GeometryType::Pointer(new Line2D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
      mPorousElement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mPorousElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mPorousElement3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),


      // conditions
      mAcousticLoadCondition2D1N(0, Condition::GeometryType::Pointer(new Point2D<NodeType >(Condition::GeometryType::PointsArrayType(1)))),
      mAcousticLoadCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mAcousticLoadCondition3D1N(0, Condition::GeometryType::Pointer(new Point3D<NodeType >(Condition::GeometryType::PointsArrayType(1)))),
      mAcousticLoadCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mAcousticLoadCondition3D4N(0, Condition::GeometryType::Pointer(new Quadrilateral3D4<NodeType >(Condition::GeometryType::PointsArrayType(4)))),
      mAcousticRobinCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mAcousticRobinCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mAcousticRobinCondition3D4N(0, Condition::GeometryType::Pointer(new Quadrilateral3D4<NodeType >(Condition::GeometryType::PointsArrayType(4)))),
      mAcousticStructureCouplingCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mAcousticStructureCouplingCondition3D4N(0, Condition::GeometryType::Pointer(new Quadrilateral3D4<NodeType >(Condition::GeometryType::PointsArrayType(4)))),
      mAcousticStructureCouplingCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mAcousticStructureMappingCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mAcousticStructureMappingCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mAcousticStructureMappingCondition3D4N(0, Condition::GeometryType::Pointer(new Quadrilateral3D4<NodeType >(Condition::GeometryType::PointsArrayType(4)))),
      mDisplacementOutputCondition3D1N(0, Condition::GeometryType::Pointer(new Point3D<NodeType >(Condition::GeometryType::PointsArrayType(1)))),
      mPressureOutputCondition3D1N(0, Condition::GeometryType::Pointer(new Point3D<NodeType >(Condition::GeometryType::PointsArrayType(1))))
    {}

void KratosMORApplication::Register()
{
     // calling base class register to register Kratos components
     KratosApplication::Register();
     KRATOS_INFO("") << "Initializing KratosMORApplication..." << std::endl;

  KRATOS_REGISTER_VARIABLE( FREQUENCY )
  KRATOS_REGISTER_VARIABLE( ADMITTANCE )
  KRATOS_REGISTER_VARIABLE( ACOUSTIC_LOAD )
  KRATOS_REGISTER_VARIABLE( BUILD_LEVEL )
  KRATOS_REGISTER_VARIABLE( USE_FREQUENCY_DEPENDENT_MATERIAL )
  KRATOS_REGISTER_VARIABLE( SCALAR_OUTPUT )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( COMPONENT_OUTPUT )

  // complex dof variables
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( REAL_DISPLACEMENT )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( IMAG_DISPLACEMENT )
  KRATOS_REGISTER_VARIABLE( REAL_PRESSURE )
  KRATOS_REGISTER_VARIABLE( IMAG_PRESSURE )

  // Biot theory variables
  KRATOS_REGISTER_VARIABLE( LAMBDA_SOLID )
  KRATOS_REGISTER_VARIABLE( MUE_SOLID )
  KRATOS_REGISTER_VARIABLE( DAMPING_SOLID )
  KRATOS_REGISTER_VARIABLE( DENSITY_SOLID )
  KRATOS_REGISTER_VARIABLE( DENSITY_FLUID )
  KRATOS_REGISTER_VARIABLE( VISCOSITY_FLUID )
  KRATOS_REGISTER_VARIABLE( STANDARD_PRESSURE_FLUID )
  KRATOS_REGISTER_VARIABLE( HEAT_CAPACITY_FLUID )
  KRATOS_REGISTER_VARIABLE( PRANDTL_NUMBER_FLUID )
  KRATOS_REGISTER_VARIABLE( THERMAL_LENGTH )
  KRATOS_REGISTER_VARIABLE( TORTUOSITY )
  KRATOS_REGISTER_VARIABLE( FLOW_RESISTIVITY )
  KRATOS_REGISTER_VARIABLE( VISCOUS_LENGTH )

  // PML variables
  KRATOS_REGISTER_VARIABLE( PML_PRESCRIBED_POTENTIAL )
  KRATOS_REGISTER_VARIABLE( PML_FACTOR )
  KRATOS_REGISTER_VARIABLE( PML_LOCAL_WIDTH )
  KRATOS_REGISTER_VARIABLE( PML_TUNING_FREQUENCY )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( PML_DIRECTION )

  // elements
  KRATOS_REGISTER_ELEMENT("AcousticElement2D2N", mAcousticElement2D2N)
  KRATOS_REGISTER_ELEMENT("AcousticElement2D3N", mAcousticElement2D3N)
  KRATOS_REGISTER_ELEMENT("AcousticElement2D4N", mAcousticElement2D4N)
  KRATOS_REGISTER_ELEMENT("AcousticElement3D4N", mAcousticElement3D4N)
  KRATOS_REGISTER_ELEMENT("AcousticElement3D8N", mAcousticElement3D8N)
  KRATOS_REGISTER_ELEMENT("AcousticPMLElement2D2N", mAcousticPMLElement2D2N)
  KRATOS_REGISTER_ELEMENT("AcousticPMLElement2D3N", mAcousticPMLElement2D3N)
  KRATOS_REGISTER_ELEMENT("AcousticPMLElement2D4N", mAcousticPMLElement2D4N)
  KRATOS_REGISTER_ELEMENT("AcousticPMLElement3D4N", mAcousticPMLElement3D4N)
  KRATOS_REGISTER_ELEMENT("AcousticPMLElement3D8N", mAcousticPMLElement3D8N)
  KRATOS_REGISTER_ELEMENT("PorousElement2D2N", mPorousElement2D2N)
  KRATOS_REGISTER_ELEMENT("PorousElement2D4N", mPorousElement2D4N)
  KRATOS_REGISTER_ELEMENT("PorousElement3D4N", mPorousElement3D4N)
  KRATOS_REGISTER_ELEMENT("PorousElement3D8N", mPorousElement3D8N)

  // conditions
  KRATOS_REGISTER_CONDITION("AcousticLoadCondition2D1N", mAcousticLoadCondition2D1N)
  KRATOS_REGISTER_CONDITION("AcousticLoadCondition2D2N", mAcousticLoadCondition2D2N)
  KRATOS_REGISTER_CONDITION("AcousticLoadCondition3D1N", mAcousticLoadCondition3D1N)
  KRATOS_REGISTER_CONDITION("AcousticLoadCondition3D3N", mAcousticLoadCondition3D3N)
  KRATOS_REGISTER_CONDITION("AcousticLoadCondition3D4N", mAcousticLoadCondition3D4N)
  KRATOS_REGISTER_CONDITION("AcousticRobinCondition2D2N", mAcousticRobinCondition2D2N)
  KRATOS_REGISTER_CONDITION("AcousticRobinCondition3D3N", mAcousticRobinCondition3D3N)
  KRATOS_REGISTER_CONDITION("AcousticRobinCondition3D4N", mAcousticRobinCondition3D4N)
  KRATOS_REGISTER_CONDITION("AcousticStructureCouplingCondition2D2N", mAcousticStructureCouplingCondition2D2N)
  KRATOS_REGISTER_CONDITION("AcousticStructureCouplingCondition3D4N", mAcousticStructureCouplingCondition3D4N)
  KRATOS_REGISTER_CONDITION("AcousticStructureCouplingCondition3D3N", mAcousticStructureCouplingCondition3D3N)
  KRATOS_REGISTER_CONDITION("AcousticStructureMappingCondition2D2N", mAcousticStructureMappingCondition2D2N)
  KRATOS_REGISTER_CONDITION("AcousticStructureMappingCondition3D3N", mAcousticStructureMappingCondition3D3N)
  KRATOS_REGISTER_CONDITION("AcousticStructureMappingCondition3D4N", mAcousticStructureMappingCondition3D4N)

  KRATOS_REGISTER_CONDITION("DisplacementOutputCondition3D1N", mDisplacementOutputCondition3D1N)
  KRATOS_REGISTER_CONDITION("PressureOutputCondition3D1N", mPressureOutputCondition3D1N)

}
}  // namespace Kratos.
