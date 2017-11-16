//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 0.0 $
//
//   Implementation of the Gauss-Seidel two step Updated Lagrangian Velocity-Pressure element
//     ( There is a ScalingConstant to multiply the mass balance equation for a number because i read it somewhere)
//

// System includes

// External includes
 
// Project includes
#include "custom_elements/two_step_updated_lagrangian_V_P_element.h"
#include "includes/cfd_variables.h"

namespace Kratos {

  /*
   * public TwoStepUpdatedLagrangianVPElement<TDim> functions
   */


  template< unsigned int TDim >
  TwoStepUpdatedLagrangianVPElement<TDim>::TwoStepUpdatedLagrangianVPElement(TwoStepUpdatedLagrangianVPElement  const& rOther)
  :Element(rOther)
  {
    KRATOS_TRY;
    KRATOS_CATCH("");
  }


  template< unsigned int TDim >
  Element::Pointer TwoStepUpdatedLagrangianVPElement<TDim>::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {
    KRATOS_TRY;

    TwoStepUpdatedLagrangianVPElement NewElement(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );
    return Element::Pointer( new TwoStepUpdatedLagrangianVPElement(NewElement) );
    
    KRATOS_CATCH("");

  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
								     VectorType& rRightHandSideVector,
								     ProcessInfo& rCurrentProcessInfo)
  { 
    KRATOS_TRY;

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
      {
      case 1:
	{
	  this->CalculateLocalMomentumEquations(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
	  break;
	}
      case 5:
      	{
      	  this->CalculateLocalContinuityEqForPressure(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
      	  break;
      	}

      default:
	{
	  KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for TWO_STEP_UPDATED_LAGRANGIAN_V_P_ELEMENT index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
	}
      }

    KRATOS_CATCH("");
  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::EquationIdVector(EquationIdVectorType& rResult,
								 ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY;

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
      {
      case 1:
	{
	  this->VelocityEquationIdVector(rResult,rCurrentProcessInfo);
	  break;
	}
      case 5:
	{
	  this->PressureEquationIdVector(rResult,rCurrentProcessInfo);
	  break;
	}
      case 6:
	{
	  this->VelocityEquationIdVector(rResult,rCurrentProcessInfo);
	  break;
	}
      default:
	{
	  KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
	}
      }

    KRATOS_CATCH("");
  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::GetDofList(DofsVectorType& rElementalDofList,
							   ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY;

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
      {
      case 1:
	{
	  this->GetVelocityDofList(rElementalDofList,rCurrentProcessInfo);
	  break;
	}
      case 5:
	{
	  this->GetPressureDofList(rElementalDofList,rCurrentProcessInfo);
	  break;
	}
      case 6:
	{
	  this->GetVelocityDofList(rElementalDofList,rCurrentProcessInfo);
	  break;
	}
      default:
	{
	  KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
	}
      }

    KRATOS_CATCH("");
  }

  // template< unsigned int TDim >
  // GeometryData::IntegrationMethod TwoStepUpdatedLagrangianVPElement<TDim>::GetIntegrationMethod() const
  // {
  //   return GeometryData::GI_GAUSS_4;
  // }
  template< unsigned int TDim >
  GeometryData::IntegrationMethod TwoStepUpdatedLagrangianVPElement<TDim>::GetIntegrationMethod() const
  {
    return GeometryData::GI_GAUSS_1;
  }



  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::CalculateLocalMomentumEquations(MatrixType& rLeftHandSideMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY; 

    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = TDim * NumNodes;

    MatrixType MassMatrix= ZeroMatrix(LocalSize,LocalSize);
    MatrixType StiffnessMatrix= ZeroMatrix(LocalSize,LocalSize);

    // Check sizes and initialize
    if( rLeftHandSideMatrix.size1() != LocalSize )
      rLeftHandSideMatrix.resize(LocalSize,LocalSize);

    rLeftHandSideMatrix = ZeroMatrix(LocalSize,LocalSize);

    if( rRightHandSideVector.size() != LocalSize )
      rRightHandSideVector.resize(LocalSize);

    rRightHandSideVector = ZeroVector(LocalSize);

    // Shape functions and integration points
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);
    const unsigned int NumGauss = GaussWeights.size();
    const double TimeStep=rCurrentProcessInfo[DELTA_TIME];

    double theta=this->GetThetaMomentum();

    ElementalVariables rElementalVariables;
    this->InitializeElementalVariables(rElementalVariables);

    double totalVolume=0;
    double MeanValueMass=0;
    double Density=0.0;
    double DeviatoricCoeff = 0;
    double VolumetricCoeff = 0;
    // this->ComputeMaterialParameters(Density,DeviatoricCoeff,VolumetricCoeff,TimeStep);
 
    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; g++)
      {
	const double GaussWeight = GaussWeights[g];
	totalVolume+=GaussWeight;
	const ShapeFunctionsType& N = row(NContainer,g);
	const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g]; 

	double Pressure=0;
	double OldPressure=0;

	this->EvaluateInPoint(Pressure,PRESSURE,N,0);

	this->EvaluateInPoint(OldPressure,PRESSURE,N,1);

	rElementalVariables.MeanPressure=OldPressure*(1-theta)+Pressure*theta;  

	bool computeElement=this->CalcMechanicsUpdated(rElementalVariables,rCurrentProcessInfo,rDN_DX,g);

	this->ComputeMaterialParameters(Density,DeviatoricCoeff,VolumetricCoeff,TimeStep,rElementalVariables);

	this->CalcElasticPlasticCauchySplitted(rElementalVariables,TimeStep,g);

	if(computeElement==true){
	  // Add integration point contribution to the local mass matrix
	  // double DynamicWeight=GaussWeight*Density;
	  // this->ComputeMassMatrix(MassMatrix,N,DynamicWeight,MeanValueMass);

	  this->AddExternalForces(rRightHandSideVector,Density,N,GaussWeight);

	  this->AddInternalForces(rRightHandSideVector,rDN_DX,rElementalVariables,GaussWeight);

	  // double lumpedDynamicWeight=GaussWeight*Density;
	  // this->ComputeLumpedMassMatrix(MassMatrix,lumpedDynamicWeight,MeanValueMass);   

	  // double MeanValueMaterial=0.0;
	  // this->ComputeMeanValueMaterialTangentMatrix(rElementalVariables,MeanValueMaterial,rDN_DX,DeviatoricCoeff,VolumetricCoeff,GaussWeight,MeanValueMass,TimeStep);    

	  // // Add viscous term
	  // this->ComputeCompleteTangentTerm(rElementalVariables,rLeftHandSideMatrix,rDN_DX,DeviatoricCoeff,VolumetricCoeff,theta,GaussWeight);
	  this->ComputeCompleteTangentTerm(rElementalVariables,StiffnessMatrix,rDN_DX,DeviatoricCoeff,VolumetricCoeff,theta,GaussWeight);
	}
      }

    double lumpedDynamicWeight=totalVolume*Density;
    this->ComputeLumpedMassMatrix(MassMatrix,lumpedDynamicWeight,MeanValueMass);    

    double BulkReductionCoefficient=1.0;
    double MeanValueStiffness=0.0;
    this->ComputeBulkReductionCoefficient(MassMatrix,StiffnessMatrix,MeanValueStiffness,BulkReductionCoefficient,TimeStep);
    if(BulkReductionCoefficient!=1.0){
      // VolumetricCoeff*=BulkReductionCoefficient;
      VolumetricCoeff*=MeanValueMass*2.0/(TimeStep*MeanValueStiffness);
      StiffnessMatrix= ZeroMatrix(LocalSize,LocalSize);

      for (unsigned int g = 0; g < NumGauss; g++)
    	{
    	  const double GaussWeight = GaussWeights[g];
    	  const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g]; 
    	  this->ComputeCompleteTangentTerm(rElementalVariables,StiffnessMatrix,rDN_DX,DeviatoricCoeff,VolumetricCoeff,theta,GaussWeight);
    	}
    }


    // Add residual of previous iteration to RHS
    VectorType VelocityValues = ZeroVector(LocalSize);
    VectorType AccelerationValues = ZeroVector(LocalSize);

    // //1st order 
    // this->GetVelocityValues(VelocityValues,0);
    // AccelerationValues = VelocityValues/TimeStep;
    // this->GetAccelerationValues(LastAccValues,0);
    // this->GetVelocityValues(VelocityValues,1);
    // AccelerationValues += -VelocityValues/TimeStep; 
    // noalias( rRightHandSideVector ) += -prod(MassMatrix,AccelerationValues);
    // noalias( rLeftHandSideMatrix ) +=  MassMatrix/TimeStep;

    //2nd order 
    this->GetAccelerationValues(AccelerationValues,0);
    this->GetVelocityValues(VelocityValues,0);
    noalias(AccelerationValues)+=-2.0*VelocityValues/TimeStep;
    this->GetVelocityValues(VelocityValues,1);
    noalias(AccelerationValues)+=2.0*VelocityValues/TimeStep;//these are negative accelerations
    noalias( rRightHandSideVector )+= prod(MassMatrix,AccelerationValues);
    noalias( rLeftHandSideMatrix ) +=  StiffnessMatrix + MassMatrix*2/TimeStep;


    // // Add residual of previous iteration to RHS
    // VectorType VelocityValues = ZeroVector(LocalSize);
    // VectorType UpdatedAccelerations = ZeroVector(LocalSize);
    // VectorType LastAccValues = ZeroVector(LocalSize);

    // // //1st order 
    // // this->GetVelocityValues(VelocityValues,0);
    // // UpdatedAccelerations = VelocityValues/TimeStep;
    // // this->GetAccelerationValues(LastAccValues,0);
    // // this->GetVelocityValues(VelocityValues,1);
    // // UpdatedAccelerations += -VelocityValues/TimeStep; 
    // // // UpdatedAccelerations =LastAccValues;
    // // noalias( rRightHandSideVector ) += -prod(MassMatrix,UpdatedAccelerations);
    // // noalias( rLeftHandSideMatrix ) +=  MassMatrix/TimeStep;

    // //2nd order 
    // this->GetVelocityValues(VelocityValues,0);
    // UpdatedAccelerations = 2.0*VelocityValues/TimeStep;
    // this->GetAccelerationValues(LastAccValues,0);
    // this->GetVelocityValues(VelocityValues,1);
    // UpdatedAccelerations += -2.0*VelocityValues/TimeStep - LastAccValues; 
    // noalias( rRightHandSideVector ) += -prod(MassMatrix,UpdatedAccelerations);
    // noalias( rLeftHandSideMatrix ) +=  StiffnessMatrix;
    // noalias( rLeftHandSideMatrix ) +=  MassMatrix*2/TimeStep;

    KRATOS_CATCH( "" );
 
  }


  template<>
  void TwoStepUpdatedLagrangianVPElement<2>::ComputeCompleteTangentTerm(ElementalVariables & rElementalVariables,
									MatrixType& rDampingMatrix,
									const ShapeFunctionDerivativesType& rDN_DX,
									const double secondLame,
									const double bulkModulus,
									const double theta,
									const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    SizeType FirstRow=0;
    SizeType FirstCol=0;

  for (SizeType j = 0; j < NumNodes; ++j)
      {
        for (SizeType i = 0; i < NumNodes; ++i)
	  {
	    double lagDNXi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,0);
	    double lagDNYi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,1);
	    double lagDNXj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,0);
	    double lagDNYj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,1);
	    // lagDNXi=rDN_DX(i,0);
	    // lagDNYi=rDN_DX(i,1);
	    // lagDNXj=rDN_DX(j,0);
	    // lagDNYj=rDN_DX(j,1);


            // First Row
            rDampingMatrix(FirstRow,FirstCol) += Weight * ( (FourThirds * secondLame + bulkModulus)*  lagDNXi * lagDNXj + lagDNYi * lagDNYj * secondLame ) *theta;
            rDampingMatrix(FirstRow,FirstCol+1) += Weight * ( (nTwoThirds* secondLame + bulkModulus) *  lagDNXi * lagDNYj + lagDNYi * lagDNXj * secondLame )*theta;

            // Second Row
            rDampingMatrix(FirstRow+1,FirstCol) += Weight * ( (nTwoThirds * secondLame + bulkModulus) * lagDNYi * lagDNXj + lagDNXi * lagDNYj *  secondLame )*theta;
            rDampingMatrix(FirstRow+1,FirstCol+1) += Weight * ( (FourThirds * secondLame + bulkModulus) * lagDNYi * lagDNYj + lagDNXi * lagDNXj * secondLame )*theta;

            // Update Counter
            FirstRow += 2;
	  }
        FirstRow = 0;
        FirstCol += 2;
      }
  }
 

  template<>
  void TwoStepUpdatedLagrangianVPElement<3>::ComputeCompleteTangentTerm(ElementalVariables & rElementalVariables,
									MatrixType& rDampingMatrix,
									const ShapeFunctionDerivativesType& rDN_DX,
									const double secondLame,
									const double bulkModulus,
									const double theta,
									const double Weight){
    
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    SizeType FirstRow=0;
    SizeType FirstCol=0;

    for (SizeType j = 0; j < NumNodes; ++j)
      {
        for (SizeType i = 0; i < NumNodes; ++i)
	  {
	    double lagDNXi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,0)+rDN_DX(i,2)*rElementalVariables.InvFgrad(2,0);
	    double lagDNYi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,1)+rDN_DX(i,2)*rElementalVariables.InvFgrad(2,1);
	    double lagDNZi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,2)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,2)+rDN_DX(i,2)*rElementalVariables.InvFgrad(2,2);
	    double lagDNXj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,0)+rDN_DX(j,2)*rElementalVariables.InvFgrad(2,0);
	    double lagDNYj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,1)+rDN_DX(j,2)*rElementalVariables.InvFgrad(2,1);
	    double lagDNZj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,2)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,2)+rDN_DX(j,2)*rElementalVariables.InvFgrad(2,2);	 
	    // lagDNXi=rDN_DX(i,0);
	    // lagDNYi=rDN_DX(i,1);
	    // lagDNZi=rDN_DX(i,2);
	    // lagDNXj=rDN_DX(j,0);
	    // lagDNYj=rDN_DX(j,1);
	    // lagDNZj=rDN_DX(j,2);

            // First Row
            rDampingMatrix(FirstRow,FirstCol)     += Weight * ( (FourThirds * secondLame + bulkModulus)*  lagDNXi * lagDNXj + (lagDNYi * lagDNYj +lagDNZi * lagDNZj) * secondLame ) *theta;
            rDampingMatrix(FirstRow,FirstCol+1)   += Weight * ( (nTwoThirds* secondLame + bulkModulus) *  lagDNXi * lagDNYj + lagDNYi * lagDNXj * secondLame )*theta;
            rDampingMatrix(FirstRow,FirstCol+2)   += Weight * ( (nTwoThirds* secondLame + bulkModulus) *  lagDNXi * lagDNZj + lagDNZi * lagDNXj * secondLame )*theta;

            // Second Row
            rDampingMatrix(FirstRow+1,FirstCol)   += Weight * ( (nTwoThirds * secondLame + bulkModulus) * lagDNYi * lagDNXj + lagDNXi * lagDNYj *  secondLame )*theta;
            rDampingMatrix(FirstRow+1,FirstCol+1) += Weight * ( (FourThirds * secondLame + bulkModulus) * lagDNYi * lagDNYj + (lagDNXi * lagDNXj + lagDNZi * lagDNZj) * secondLame )*theta;
            rDampingMatrix(FirstRow+1,FirstCol+2) += Weight * ( (nTwoThirds * secondLame + bulkModulus) * lagDNYi * lagDNZj + lagDNZi * lagDNYj *  secondLame )*theta;

            // Third Row
            rDampingMatrix(FirstRow+2,FirstCol)   += Weight * ( (nTwoThirds * secondLame + bulkModulus) * lagDNZi * lagDNXj + lagDNXi * lagDNZj *  secondLame )*theta;
            rDampingMatrix(FirstRow+2,FirstCol+1) += Weight * ( (nTwoThirds* secondLame + bulkModulus)  * lagDNZi * lagDNYj + lagDNYi * lagDNZj *  secondLame )*theta;
            rDampingMatrix(FirstRow+2,FirstCol+2) += Weight * ( (FourThirds * secondLame + bulkModulus) * lagDNZi * lagDNZj + (lagDNXi * lagDNXj + lagDNYi * lagDNYj) * secondLame )*theta;
	   
	    // Update Counter
            FirstRow += 3;
	  }
        FirstRow = 0;
        FirstCol += 3;
      }
  }




  template< unsigned int TDim >
  int TwoStepUpdatedLagrangianVPElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo)
  {
    KRATOS_TRY;

    // Base class checks for positive Jacobian and Id > 0
    int ierr = Element::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Check that all required variables have been registered
    if(VELOCITY.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check that the application was correctly registered.","");
    if(ACCELERATION.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"ACCELERATION Key is 0. Check that the application was correctly registered.","");
    if(PRESSURE.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"PRESSURE Key is 0. Check that the application was correctly registered.","");
    if(BODY_FORCE.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"BODY_FORCE Key is 0. Check that the application was correctly registered.","");
    if(DENSITY.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY Key is 0. Check that the application was correctly registered.","");
    if(VISCOSITY.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"VISCOSITY Key is 0. Check that the application was correctly registered.","");
    if(DELTA_TIME.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"DELTA_TIME Key is 0. Check that the application was correctly registered.","");

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
      {
        if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(BODY_FORCE) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing BODY_FORCE variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(DENSITY) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing DENSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(VISCOSITY) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing VISCOSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
           this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
           this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE component degree of freedom on node ",this->GetGeometry()[i].Id());
      }
    
    // If this is a 2D problem, check that nodes are in XY plane
    if (this->GetGeometry().WorkingSpaceDimension() == 2)
      {
        for (unsigned int i=0; i<this->GetGeometry().size(); ++i)
	  {
            if (this->GetGeometry()[i].Z() != 0.0)
	      KRATOS_THROW_ERROR(std::invalid_argument,"Node with non-zero Z coordinate found. Id: ",this->GetGeometry()[i].Id());
	  }
      }

    return ierr;

    KRATOS_CATCH("");
  }



  template<>
  void TwoStepUpdatedLagrangianVPElement<2>::VelocityEquationIdVector(EquationIdVectorType& rResult,
								      ProcessInfo& rCurrentProcessInfo)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = NumNodes*2;

    SizeType LocalIndex = 0;

    if (rResult.size() != LocalSize)
      rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
      }

  }

  template<>
  void TwoStepUpdatedLagrangianVPElement<3>::VelocityEquationIdVector(EquationIdVectorType& rResult,
								      ProcessInfo& rCurrentProcessInfo)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3*NumNodes;

    SizeType LocalIndex = 0;

    if (rResult.size() != LocalSize)
      rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Z,xpos+2).EquationId();
      }
  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::PressureEquationIdVector(EquationIdVectorType& rResult,
									 ProcessInfo& rCurrentProcessInfo)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rResult.size() != NumNodes)
      rResult.resize(NumNodes);

    const unsigned int pos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);

    for (SizeType i = 0; i < NumNodes; ++i)
      rResult[i] = rGeom[i].GetDof(PRESSURE,pos).EquationId();
  }

  template<>
  void TwoStepUpdatedLagrangianVPElement<2>::GetVelocityDofList(DofsVectorType& rElementalDofList,
								ProcessInfo& rCurrentProcessInfo)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2*NumNodes;

    if (rElementalDofList.size() != LocalSize)
      rElementalDofList.resize(LocalSize);

    SizeType LocalIndex = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
      }
  }

  template<>
  void TwoStepUpdatedLagrangianVPElement<3>::GetVelocityDofList(DofsVectorType& rElementalDofList,
								ProcessInfo& rCurrentProcessInfo)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3*NumNodes;

    if (rElementalDofList.size() != LocalSize)
      rElementalDofList.resize(LocalSize);

    SizeType LocalIndex = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Z);
      }
  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::GetPressureDofList(DofsVectorType& rElementalDofList,
								   ProcessInfo& rCurrentProcessInfo)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rElementalDofList.size() != NumNodes)
      rElementalDofList.resize(NumNodes);

    SizeType LocalIndex = 0;
    for (SizeType i = 0; i < NumNodes; ++i)
      {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
      }
	
  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::GetPressureValues(Vector& rValues,
								  const int Step)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rValues.size() != NumNodes) rValues.resize(NumNodes);

    for (SizeType i = 0; i < NumNodes; ++i){
      rValues[i] = rGeom[i].FastGetSolutionStepValue(PRESSURE,Step);

      if(rGeom[i].Is(FREE_SURFACE)){
	rGeom[i].FastGetSolutionStepValue(FREESURFACE) = 1;

      }else{
      	rGeom[i].FastGetSolutionStepValue(FREESURFACE) = 0;
      }

    }
  }



  template< unsigned int TDim>
  bool TwoStepUpdatedLagrangianVPElement<TDim>::CalcMechanicsUpdated(ElementalVariables & rElementalVariables,
								     const ProcessInfo& rCurrentProcessInfo,
								     const ShapeFunctionDerivativesType& rDN_DX,
								     unsigned int g)
  {
  
    double theta=this->GetThetaMomentum();
    // bool computeElement=this->CalcStrainRate(rElementalVariables,rCurrentProcessInfo,rDN_DX,theta);
    bool computeElement=this->CalcCompleteStrainRate(rElementalVariables,rCurrentProcessInfo,rDN_DX,theta);
    return computeElement;
  
  }


  template<  unsigned int TDim>
  void TwoStepUpdatedLagrangianVPElement<TDim>::CalcMeanPressure(double &meanPressure,
								 const int Step)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    meanPressure=0;
    double coeff=0.0;
    for (SizeType i = 0; i < NumNodes; ++i)
      {
	meanPressure+=rGeom[i].GetSolutionStepValue(PRESSURE,Step);
	coeff+=1.0;
      }
    meanPressure*=1.0/coeff;
  }



  template< >
  void TwoStepUpdatedLagrangianVPElement<2>::CalcMeanVelocity(double &meanVelocity,
							      const int Step)
  {

    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    //SizeType Index = 0;
   
    double velX=0;
    double velY=0;
    double coeff=0.0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        velX += rGeom[i].FastGetSolutionStepValue(VELOCITY_X,Step);
        velY += rGeom[i].FastGetSolutionStepValue(VELOCITY_Y,Step);
	coeff+=1.0;
      }

    meanVelocity=velX*velX+velY*velY;
    meanVelocity=sqrt(meanVelocity)/coeff;

  }



  template< >
  void TwoStepUpdatedLagrangianVPElement<3>::CalcMeanVelocity(double &meanVelocity,
							      const int Step)
  {

    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    //SizeType Index = 0;
   
    double velX=0;
    double velY=0;
    double velZ=0;
    double coeff=0.0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        velX += rGeom[i].FastGetSolutionStepValue(VELOCITY_X,Step);
        velY += rGeom[i].FastGetSolutionStepValue(VELOCITY_Y,Step);
        velZ += rGeom[i].FastGetSolutionStepValue(VELOCITY_Z,Step);
	coeff+=1.0;
      }

    meanVelocity=velX*velX+velY*velY+velZ*velZ;
    meanVelocity=sqrt(meanVelocity)/coeff;

  }





template< unsigned int TDim>
void TwoStepUpdatedLagrangianVPElement<TDim>::CalculateDeltaPosition(Matrix & rDeltaPosition)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = this->GetGeometry().PointsNumber();
    unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

    rDeltaPosition = zero_matrix<double>( number_of_nodes , dimension);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        array_1d<double, 3 > & CurrentDisplacement  = this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        array_1d<double, 3 > & PreviousDisplacement = this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rDeltaPosition(i,j) = CurrentDisplacement[j]-PreviousDisplacement[j];
        }
    }

    KRATOS_CATCH( "" )
}

  template<>
  void TwoStepUpdatedLagrangianVPElement<2>::GetDisplacementValues(Vector& rValues,
								   const int Step)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2*NumNodes;

    if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_X,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_Y,Step);
      }
  }


  template<>
  void TwoStepUpdatedLagrangianVPElement<3>::GetDisplacementValues(Vector& rValues,
								   const int Step)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3*NumNodes;

    if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_X,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_Y,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_Z,Step);
      }
  }



  template<>
  void TwoStepUpdatedLagrangianVPElement<2>::GetVelocityValues(Vector& rValues,
							       const int Step)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2*NumNodes;

    if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_X,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y,Step);
      }
  }


  template<>
  void TwoStepUpdatedLagrangianVPElement<3>::GetVelocityValues(Vector& rValues,
							       const int Step)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3*NumNodes;

    if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_X,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Z,Step);
      }
  }



  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::GetElementalAcceleration(Vector& meanAcceleration,
									 const int Step,
									 const double TimeStep)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    double count=0;
    for (SizeType i = 0; i < NumNodes; ++i)
      {
	meanAcceleration+= 0.5/TimeStep*(rGeom[i].FastGetSolutionStepValue(VELOCITY,0)-rGeom[i].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[i].FastGetSolutionStepValue(ACCELERATION,1); 
	// meanAcceleration+=rGeom[i].FastGetSolutionStepValue(ACCELERATION,step);
	count+=1.0;
      }
    meanAcceleration*=1.0/count;
  }

  template< >
  void TwoStepUpdatedLagrangianVPElement<2>::GetAccelerationValues(Vector& rValues,
								   const int Step)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2*NumNodes;

    if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(ACCELERATION_X,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(ACCELERATION_Y,Step);
      }
  }

  template< >
  void TwoStepUpdatedLagrangianVPElement<3>::GetAccelerationValues(Vector& rValues,
								   const int Step)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3*NumNodes;

    if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(ACCELERATION_X,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(ACCELERATION_Y,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(ACCELERATION_Z,Step);
      }
  }



  // template< unsigned int TDim >
  // void TwoStepUpdatedLagrangianVPElement<TDim>::CalculateGeometryData(ShapeFunctionDerivativesArrayType &rDN_DX,
  // 								      Matrix &NContainer,
  // 								      Vector &rGaussWeights)
  // {
  //   const GeometryType& rGeom = this->GetGeometry();
  //   Vector DetJ;
  //   rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX,DetJ,GeometryData::GI_GAUSS_4);
  //   NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_4);
  //   const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_4);

  //   rGaussWeights.resize(rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_4),false);

  //   for (unsigned int g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_4); g++){
  //     rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
  //     if(rGaussWeights[g]<0)
  // 	std::cout<<"NEGATIVE GAUSS WEIGHT "<<rGaussWeights[g]<<std::endl;
  
  //   }

  // }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::CalculateGeometryData(ShapeFunctionDerivativesArrayType &rDN_DX,
  								      Matrix &NContainer,
  								      Vector &rGaussWeights)
  {
    const GeometryType& rGeom = this->GetGeometry();
    Vector DetJ;
    rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX,DetJ,GeometryData::GI_GAUSS_1);
    NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_1);

    rGaussWeights.resize(rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_1),false);

    for (unsigned int g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_1); g++){
      // rGaussWeights[g] = fabs(DetJ[g] * IntegrationPoints[g].Weight());
      rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
      if(rGaussWeights[g]<0)
    	std::cout<<"NEGATIVE GAUSS WEIGHT "<<rGaussWeights[g]<<std::endl;
    }
  }

  template< unsigned int TDim >
  double TwoStepUpdatedLagrangianVPElement<TDim>::ElementSize()
  {
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();


    double ElemSize =0;
    array_1d<double,3> Edge(3,0.0);
    noalias(Edge) = rGeom[1].Coordinates() - rGeom[0].Coordinates();
    double Length = Edge[0]*Edge[0];
    for (SizeType d = 1; d < TDim; d++)
      Length += Edge[d]*Edge[d];
    ElemSize+=sqrt(Length);

    double count=1.0;
    for (SizeType i = 2; i < NumNodes; i++){
      for(SizeType j = 0; j < i; j++)
        {
    	  noalias(Edge) = rGeom[i].Coordinates() - rGeom[j].Coordinates();
    	  Length = Edge[0]*Edge[0];
    	  for (SizeType d = 1; d < TDim; d++){
    	    Length += Edge[d]*Edge[d];
    	  }
    	  ElemSize+=sqrt(Length);
    	  count+=1.0;
        }
    }
    ElemSize*=1.0/count; 
    return ElemSize;

    // // calculate minimum element length (used in stabilization Tau)
    // array_1d<double,3> Edge(3,0.0);
    // Edge = rGeom[1].Coordinates() - rGeom[0].Coordinates();
    // double ElemSize = Edge[0]*Edge[0];
    // for (SizeType d = 1; d < TDim; d++)
    //   ElemSize += Edge[d]*Edge[d];

    // for (SizeType i = 2; i < NumNodes; i++)
    //   for(SizeType j = 0; j < i; j++)
    //     {
    // 	  Edge = rGeom[i].Coordinates() - rGeom[j].Coordinates();
    // 	  //for computing minimum distance
    // 	  double Length = Edge[0]*Edge[0];
    // 	  for (SizeType d = 1; d < TDim; d++){
    // 	    Length += Edge[d]*Edge[d];
    // 	  }
    // 	  if (Length < ElemSize) ElemSize = Length;
    //     }

    // return sqrt(ElemSize);

  }


template < > 
bool TwoStepUpdatedLagrangianVPElement<2>::CalcCompleteStrainRate(ElementalVariables & rElementalVariables,
								  const ProcessInfo &rCurrentProcessInfo,
								  const ShapeFunctionDerivativesType& rDN_DX,
								  const double theta)
{
  bool computeElement=true;
  unsigned int dimension=this->GetGeometry().WorkingSpaceDimension();
  GeometryType& rGeom = this->GetGeometry();
  const SizeType NumNodes = rGeom.PointsNumber();
  const SizeType LocalSize = dimension*NumNodes;
  VectorType  NodePosition= ZeroVector(LocalSize);
  VectorType VelocityValues = ZeroVector(LocalSize);
  VectorType RHSVelocities = ZeroVector(LocalSize);
  this->GetPositions(NodePosition,rCurrentProcessInfo,theta);
  this->GetVelocityValues(RHSVelocities,0); 
  RHSVelocities*=theta;
  this->GetVelocityValues(VelocityValues,1);
  RHSVelocities+=VelocityValues*(1.0-theta);

  rElementalVariables.Fgrad=ZeroMatrix(dimension,dimension);
  rElementalVariables.FgradVel=ZeroMatrix(dimension,dimension);
  for (SizeType i = 0; i < dimension; i++)
    {
      for (SizeType j = 0; j < dimension; j++)
	{
	  for (SizeType k = 0; k < NumNodes; k++)
	    {
	      rElementalVariables.Fgrad(i,j)+= NodePosition[dimension*k+i]*rDN_DX(k,j);
	      rElementalVariables.FgradVel(i,j)+= RHSVelocities[dimension*k+i]*rDN_DX(k,j);
	    }
	}
    }

  //Inverse
  rElementalVariables.InvFgrad=ZeroMatrix(dimension,dimension);
  rElementalVariables.DetFgrad=1;
  MathUtils<double>::InvertMatrix(rElementalVariables.Fgrad, 
				  rElementalVariables.InvFgrad, 
				  rElementalVariables.DetFgrad);

  rElementalVariables.InvFgradVel=ZeroMatrix(dimension,dimension);
  rElementalVariables.DetFgradVel=1;
  MathUtils<double>::InvertMatrix(rElementalVariables.FgradVel, 
				  rElementalVariables.InvFgradVel, 
				  rElementalVariables.DetFgradVel);


  //it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
  rElementalVariables.SpatialVelocityGrad.resize(dimension,dimension);
  rElementalVariables.SpatialVelocityGrad=prod(rElementalVariables.FgradVel,rElementalVariables.InvFgrad);
  
  rElementalVariables.VolumetricDefRate=0;
  for (SizeType i = 0; i < dimension; i++)
    {
      rElementalVariables.VolumetricDefRate+=rElementalVariables.SpatialVelocityGrad(i,i);
    }
  
  // //it checks whether tr(l) == div(v)
  // CheckStrain1(rElementalVariables.VolumetricDefRate,rElementalVariables.SpatialVelocityGrad);
  // CheckStrain2(rElementalVariables.SpatialVelocityGrad,rElementalVariables.Fgrad,ElementalVariables.FgradVel);
  // //it computes Material time Derivative of Green Lagrange strain tensor in MATERIAL configuration --> [D(E)/Dt]
  // // x-component
  // rElementalVariables.MDGreenLagrangeMaterial[0]=rElementalVariables.FgradVel(0,0)*rElementalVariables.Fgrad(0,0) + 
  //   rElementalVariables.FgradVel(1,0)*rElementalVariables.Fgrad(1,0);
  // // y-component
  // rElementalVariables.MDGreenLagrangeMaterial[1]=rElementalVariables.FgradVel(1,1)*rElementalVariables.Fgrad(1,1) + 
  //   rElementalVariables.FgradVel(0,1)*rElementalVariables.Fgrad(0,1);
  // // xy-component
  // rElementalVariables.MDGreenLagrangeMaterial[2]=(rElementalVariables.FgradVel(0,0)*rElementalVariables.Fgrad(0,1) + 
  // 						  rElementalVariables.FgradVel(1,0)*rElementalVariables.Fgrad(1,1) +
  // 						  rElementalVariables.FgradVel(0,1)*rElementalVariables.Fgrad(0,0) + 
  // 						  rElementalVariables.FgradVel(1,1)*rElementalVariables.Fgrad(1,0))*0.5;
  // //it computes Material time Derivative of Green Lagrange strain tensor in SPATIAL configuration  --> [d]
  // // x-component
  // rElementalVariables.SpatialDefRate[0]= rElementalVariables.InvFgrad(0,0)*rElementalVariables.MDGreenLagrangeMaterial[0]*rElementalVariables.InvFgrad(0,0) + 
  //   rElementalVariables.InvFgrad(1,0)*rElementalVariables.MDGreenLagrangeMaterial[2]*rElementalVariables.InvFgrad(0,0)*2 +
  //   rElementalVariables.InvFgrad(1,0)*rElementalVariables.MDGreenLagrangeMaterial[1]*rElementalVariables.InvFgrad(1,0);
  // // y-component
  // rElementalVariables.SpatialDefRate[1]= rElementalVariables.InvFgrad(0,1)*rElementalVariables.MDGreenLagrangeMaterial[0]*rElementalVariables.InvFgrad(0,1) + 
  //   rElementalVariables.InvFgrad(0,1)*rElementalVariables.MDGreenLagrangeMaterial[2]*rElementalVariables.InvFgrad(1,1)*2 +
  //   rElementalVariables.InvFgrad(1,1)*rElementalVariables.MDGreenLagrangeMaterial[1]*rElementalVariables.InvFgrad(1,1);
  // // xy-component
  // rElementalVariables.SpatialDefRate[2]=rElementalVariables.InvFgrad(0,0)*rElementalVariables.MDGreenLagrangeMaterial[0]*rElementalVariables.InvFgrad(0,1) + 
  //   rElementalVariables.InvFgrad(0,0)*rElementalVariables.MDGreenLagrangeMaterial[2]*rElementalVariables.InvFgrad(1,1) +
  //   rElementalVariables.InvFgrad(1,0)*rElementalVariables.MDGreenLagrangeMaterial[2]*rElementalVariables.InvFgrad(0,1) +
  //   rElementalVariables.InvFgrad(1,0)*rElementalVariables.MDGreenLagrangeMaterial[1]*rElementalVariables.InvFgrad(1,1);
 
  // computeElement=CheckStrain3(rElementalVariables.SpatialDefRate,rElementalVariables.SpatialVelocityGrad);


  rElementalVariables.SpatialDefRate[0]=rElementalVariables.SpatialVelocityGrad(0,0);
  rElementalVariables.SpatialDefRate[1]=rElementalVariables.SpatialVelocityGrad(1,1);
  rElementalVariables.SpatialDefRate[2]=0.5*(rElementalVariables.SpatialVelocityGrad(1,0)+rElementalVariables.SpatialVelocityGrad(0,1));

  double aThird=1.0/3.0;
  double dev_X=rElementalVariables.SpatialDefRate[0]-
    (rElementalVariables.SpatialDefRate[0]+rElementalVariables.SpatialDefRate[1])*aThird;
  double dev_Y=rElementalVariables.SpatialDefRate[1]-
    (rElementalVariables.SpatialDefRate[0]+rElementalVariables.SpatialDefRate[1])*aThird;
  rElementalVariables.DeviatoricInvariant=sqrt(2*(dev_X*dev_X + dev_Y*dev_Y +
						  rElementalVariables.SpatialDefRate[2]*rElementalVariables.SpatialDefRate[2]));

  rElementalVariables.EquivalentStrainRate=sqrt((2.0*rElementalVariables.SpatialDefRate[0]*rElementalVariables.SpatialDefRate[0] +
  			 2.0*rElementalVariables.SpatialDefRate[1]*rElementalVariables.SpatialDefRate[1] +
  			 4.0*rElementalVariables.SpatialDefRate[2]*rElementalVariables.SpatialDefRate[2]));

  return computeElement;

}  

template < > 
bool TwoStepUpdatedLagrangianVPElement<3>::CalcCompleteStrainRate(ElementalVariables & rElementalVariables,
								  const ProcessInfo &rCurrentProcessInfo,
								  const ShapeFunctionDerivativesType& rDN_DX,
								  const double theta)
{

  bool computeElement=true;
  unsigned int dimension=this->GetGeometry().WorkingSpaceDimension();
  GeometryType& rGeom = this->GetGeometry();
  const SizeType NumNodes = rGeom.PointsNumber();
  const SizeType LocalSize = dimension*NumNodes;
  VectorType  NodePosition= ZeroVector(LocalSize);
  VectorType VelocityValues = ZeroVector(LocalSize);
  VectorType RHSVelocities = ZeroVector(LocalSize);
  this->GetPositions(NodePosition,rCurrentProcessInfo,theta);
  this->GetVelocityValues(RHSVelocities,0); 
  RHSVelocities*=theta;
  this->GetVelocityValues(VelocityValues,1);
  RHSVelocities+=VelocityValues*(1.0-theta);

  rElementalVariables.Fgrad=ZeroMatrix(dimension,dimension);
  rElementalVariables.FgradVel=ZeroMatrix(dimension,dimension);
  for (SizeType i = 0; i < dimension; i++)
    {
      for (SizeType j = 0; j < dimension; j++)
	{
	  for (SizeType k = 0; k < NumNodes; k++)
	    {
	      rElementalVariables.Fgrad(i,j)+= NodePosition[dimension*k+i]*rDN_DX(k,j);
	      rElementalVariables.FgradVel(i,j)+= RHSVelocities[dimension*k+i]*rDN_DX(k,j);
	    }
	}
    }

  //Inverse
  rElementalVariables.InvFgrad=ZeroMatrix(dimension,dimension);
  rElementalVariables.DetFgrad=1;
  MathUtils<double>::InvertMatrix(rElementalVariables.Fgrad, 
				  rElementalVariables.InvFgrad, 
				  rElementalVariables.DetFgrad);

  rElementalVariables.InvFgradVel=ZeroMatrix(dimension,dimension);
  rElementalVariables.DetFgradVel=1;
  MathUtils<double>::InvertMatrix(rElementalVariables.FgradVel, 
				  rElementalVariables.InvFgradVel, 
				  rElementalVariables.DetFgradVel);


  //it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
  rElementalVariables.SpatialVelocityGrad.resize(dimension,dimension);
  rElementalVariables.SpatialVelocityGrad=prod(rElementalVariables.FgradVel,rElementalVariables.InvFgrad);
  
  rElementalVariables.VolumetricDefRate=0;
  for (SizeType i = 0; i < dimension; i++)
    {
      rElementalVariables.VolumetricDefRate+=rElementalVariables.SpatialVelocityGrad(i,i);
    }
  
  // //it checks whether tr(l) == div(v)
  // CheckStrain1(rElementalVariables.VolumetricDefRate,rElementalVariables.SpatialVelocityGrad);
  // CheckStrain2(rElementalVariables.SpatialVelocityGrad,rElementalVariables.Fgrad,ElementalVariables.FgradVel);
 
  // //it computes Material time Derivative of Green Lagrange strain tensor in MATERIAL configuration --> [D(E)/Dt]
  // MatrixType MatrixA= ZeroMatrix(3,3);
  // MatrixType MatrixB= ZeroMatrix(3,3);
  // MatrixType Matrix1= ZeroMatrix(3,3);
  // MatrixType Matrix2= ZeroMatrix(3,3);

  // MatrixA=rElementalVariables.Fgrad;
  // MatrixA(0,1)=rElementalVariables.Fgrad(1,0);
  // MatrixA(0,2)=rElementalVariables.Fgrad(2,0);
  // MatrixA(1,0)=rElementalVariables.Fgrad(0,1);
  // MatrixA(1,2)=rElementalVariables.Fgrad(2,1);
  // MatrixA(2,0)=rElementalVariables.Fgrad(0,2);
  // MatrixA(2,1)=rElementalVariables.Fgrad(1,2);

  // MatrixB=rElementalVariables.FgradVel;
  // MatrixB(0,1)=rElementalVariables.FgradVel(1,0);
  // MatrixB(0,2)=rElementalVariables.FgradVel(2,0);
  // MatrixB(1,0)=rElementalVariables.FgradVel(0,1);
  // MatrixB(1,2)=rElementalVariables.FgradVel(2,1);
  // MatrixB(2,0)=rElementalVariables.FgradVel(0,2);
  // MatrixB(2,1)=rElementalVariables.FgradVel(1,2);

  // noalias(Matrix1)=prod(MatrixB,rElementalVariables.Fgrad);
  // noalias(Matrix2)=prod(MatrixA,rElementalVariables.FgradVel);

  // rElementalVariables.MDGreenLagrangeMaterial[0]= ( Matrix1(0,0) + Matrix2(0,0) ) * 0.5;  //xx-component
  // rElementalVariables.MDGreenLagrangeMaterial[1]= ( Matrix1(1,1) + Matrix2(1,1) ) * 0.5;  //yy-component
  // rElementalVariables.MDGreenLagrangeMaterial[2]= ( Matrix1(2,2) + Matrix2(2,2) ) * 0.5;  //zz-component
  // rElementalVariables.MDGreenLagrangeMaterial[3]= ( Matrix1(0,1) + Matrix2(0,1) ) * 0.5;  //xy-component
  // rElementalVariables.MDGreenLagrangeMaterial[4]= ( Matrix1(0,2) + Matrix2(0,2) ) * 0.5;  //xz-component
  // rElementalVariables.MDGreenLagrangeMaterial[5]= ( Matrix1(1,2) + Matrix2(1,2) ) * 0.5;  //yz-component


  // //it computes Material time Derivative of Green Lagrange strain tensor in SPATIAL configuration  --> [d]
  // MatrixA=rElementalVariables.InvFgrad;
  // MatrixA(0,1)=rElementalVariables.InvFgrad(1,0);
  // MatrixA(0,2)=rElementalVariables.InvFgrad(2,0);
  // MatrixA(1,0)=rElementalVariables.InvFgrad(0,1);
  // MatrixA(1,2)=rElementalVariables.InvFgrad(2,1);
  // MatrixA(2,0)=rElementalVariables.InvFgrad(0,2);
  // MatrixA(2,1)=rElementalVariables.InvFgrad(1,2);

  // MatrixB(0,0)=rElementalVariables.MDGreenLagrangeMaterial[0];  //XX-component;
  // MatrixB(1,1)=rElementalVariables.MDGreenLagrangeMaterial[1];  //YY-component;
  // MatrixB(2,2)=rElementalVariables.MDGreenLagrangeMaterial[2];  //ZZ-component;
  // MatrixB(0,1)=rElementalVariables.MDGreenLagrangeMaterial[3];  //XY-component;
  // MatrixB(1,0)=rElementalVariables.MDGreenLagrangeMaterial[3];  //XY-component;
  // MatrixB(0,2)=rElementalVariables.MDGreenLagrangeMaterial[4];  //ZX-component;
  // MatrixB(2,0)=rElementalVariables.MDGreenLagrangeMaterial[4];  //ZX-component;
  // MatrixB(1,2)=rElementalVariables.MDGreenLagrangeMaterial[5];  //YZ-component;
  // MatrixB(2,1)=rElementalVariables.MDGreenLagrangeMaterial[5];  //YZ-component;

  // noalias(Matrix1)=prod(MatrixB,rElementalVariables.InvFgrad);
  // noalias(Matrix2)=prod(MatrixA,Matrix1);
 
  // rElementalVariables.SpatialDefRate[0]=Matrix2(0,0);
  // rElementalVariables.SpatialDefRate[1]=Matrix2(1,1);
  // rElementalVariables.SpatialDefRate[2]=Matrix2(2,2);
  // rElementalVariables.SpatialDefRate[3]=Matrix2(0,1);
  // rElementalVariables.SpatialDefRate[4]=Matrix2(0,2);
  // rElementalVariables.SpatialDefRate[5]=Matrix2(1,2);


  rElementalVariables.SpatialDefRate[0]=rElementalVariables.SpatialVelocityGrad(0,0);
  rElementalVariables.SpatialDefRate[1]=rElementalVariables.SpatialVelocityGrad(1,1);
  rElementalVariables.SpatialDefRate[2]=rElementalVariables.SpatialVelocityGrad(2,2);
  rElementalVariables.SpatialDefRate[3]=0.5*(rElementalVariables.SpatialVelocityGrad(1,0)+rElementalVariables.SpatialVelocityGrad(0,1));
  rElementalVariables.SpatialDefRate[4]=0.5*(rElementalVariables.SpatialVelocityGrad(2,0)+rElementalVariables.SpatialVelocityGrad(0,2));
  rElementalVariables.SpatialDefRate[5]=0.5*(rElementalVariables.SpatialVelocityGrad(2,1)+rElementalVariables.SpatialVelocityGrad(1,2));
  // computeElement=CheckStrain3(rElementalVariables.SpatialDefRate,rElementalVariables.SpatialVelocityGrad);

  double aThird=1.0/3.0;
  double dev_X=rElementalVariables.SpatialDefRate[0]-
    (rElementalVariables.SpatialDefRate[0]+rElementalVariables.SpatialDefRate[1]+rElementalVariables.SpatialDefRate[2])*aThird;
  double dev_Y=rElementalVariables.SpatialDefRate[1]-
    (rElementalVariables.SpatialDefRate[0]+rElementalVariables.SpatialDefRate[1]+rElementalVariables.SpatialDefRate[2])*aThird;
  double dev_Z=rElementalVariables.SpatialDefRate[2]-
    (rElementalVariables.SpatialDefRate[0]+rElementalVariables.SpatialDefRate[1]+rElementalVariables.SpatialDefRate[2])*aThird;
  rElementalVariables.DeviatoricInvariant=sqrt(2*(dev_X*dev_X + dev_Y*dev_Y + dev_Z*dev_Z +
						  rElementalVariables.SpatialDefRate[3]*rElementalVariables.SpatialDefRate[3] +
						  rElementalVariables.SpatialDefRate[4]*rElementalVariables.SpatialDefRate[4] +
						  rElementalVariables.SpatialDefRate[5]*rElementalVariables.SpatialDefRate[5]));

  rElementalVariables.EquivalentStrainRate=sqrt(2.0*(rElementalVariables.SpatialDefRate[0]*rElementalVariables.SpatialDefRate[0] +
  						     rElementalVariables.SpatialDefRate[1]*rElementalVariables.SpatialDefRate[1] +
  						     rElementalVariables.SpatialDefRate[2]*rElementalVariables.SpatialDefRate[2] +
  						     2.0*rElementalVariables.SpatialDefRate[3]*rElementalVariables.SpatialDefRate[3] +
  						     2.0*rElementalVariables.SpatialDefRate[4]*rElementalVariables.SpatialDefRate[4] +
  						     2.0*rElementalVariables.SpatialDefRate[5]*rElementalVariables.SpatialDefRate[5]));


  return computeElement;

}  




template< unsigned int TDim>
bool TwoStepUpdatedLagrangianVPElement<TDim>::CalcStrainRate(ElementalVariables & rElementalVariables,
							     const ProcessInfo &rCurrentProcessInfo,
							     const ShapeFunctionDerivativesType& rDN_DX,
							     const double theta)
{


  bool computeElement=true;

  this->CalcFGrad(rDN_DX,
		  rElementalVariables.Fgrad,
		  rElementalVariables.InvFgrad,
		  rElementalVariables.DetFgrad,
		  rCurrentProcessInfo,
		  theta);

  //it computes the material time derivative of the deformation gradient and its jacobian and inverse
  this->CalcVelDefGrad(rDN_DX,
		       rElementalVariables.FgradVel,
		       rElementalVariables.InvFgradVel,
		       rElementalVariables.DetFgradVel,
		       theta);

  //it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
  this->CalcSpatialVelocityGrad(rElementalVariables.InvFgrad,
				rElementalVariables.FgradVel,
				rElementalVariables.SpatialVelocityGrad);
  
  this->CalcVolDefRateFromSpatialVelGrad(rElementalVariables.VolumetricDefRate,
  					 rElementalVariables.SpatialVelocityGrad);
  

  // this->CalcVolumetricDefRate(rDN_DX,
  // 			      rElementalVariables.VolumetricDefRate,
  // 			      rElementalVariables.InvFgrad,
  //                          theta );

  // //it checks whether tr(l) == div(v)
  // CheckStrain1(rElementalVariables.VolumetricDefRate,
  // 	       rElementalVariables.SpatialVelocityGrad);
      
  // CheckStrain2(rElementalVariables.SpatialVelocityGrad,
  // 	       rElementalVariables.Fgrad,
  // 	       rElementalVariables.FgradVel);
 
  //it computes Material time Derivative of Green Lagrange strain tensor in MATERIAL configuration --> [D(E)/Dt]
  this->CalcMDGreenLagrangeMaterial(rElementalVariables.Fgrad,
				    rElementalVariables.FgradVel,
				    rElementalVariables.MDGreenLagrangeMaterial);
      
  //it computes Material time Derivative of Green Lagrange strain tensor in SPATIAL configuration  --> [d]
  this->CalcSpatialDefRate(rElementalVariables.MDGreenLagrangeMaterial,
			   rElementalVariables.InvFgrad,
			   rElementalVariables.SpatialDefRate);



  // computeElement=CheckStrain3(rElementalVariables.SpatialDefRate,
  // 			      rElementalVariables.SpatialVelocityGrad);

  this->CalcDeviatoricInvariant(rElementalVariables.SpatialDefRate,
				rElementalVariables.DeviatoricInvariant);

  this->CalcEquivalentStrainRate(rElementalVariables.SpatialDefRate,
				 rElementalVariables.EquivalentStrainRate);
  return computeElement;

}  


  

template <unsigned int TDim > 
void TwoStepUpdatedLagrangianVPElement<TDim>::CalcFGrad(const ShapeFunctionDerivativesType& rDN_DX,
							MatrixType &Fgrad,
							MatrixType &invFgrad,
							double &FJacobian,
							const ProcessInfo& rCurrentProcessInfo,
							const double theta)
{
  GeometryType& rGeom = this->GetGeometry();
  const SizeType NumNodes = rGeom.PointsNumber();
  const SizeType LocalSize = TDim*NumNodes;
  VectorType  NodePosition= ZeroVector(LocalSize);
  this->GetPositions(NodePosition,rCurrentProcessInfo,theta);

  Fgrad.resize(TDim,TDim);

  Fgrad=ZeroMatrix(TDim,TDim);
  for (SizeType i = 0; i < TDim; i++)
    {
      for (SizeType j = 0; j < TDim; j++)
	{
	  for (SizeType k = 0; k < NumNodes; k++)
	    {
	      Fgrad(i,j)+= NodePosition[TDim*k+i]*rDN_DX(k,j);
	    }
	}
    }


 //Inverse
  invFgrad.resize(TDim,TDim);
  invFgrad=ZeroMatrix(TDim,TDim);
  FJacobian=1;


  MathUtils<double>::InvertMatrix( Fgrad, invFgrad, FJacobian );

  // Fgrad.resize(2,2);

  // Fgrad(0,0)= NodePosition[0]*rDN_DX(0,0)+NodePosition[2]*rDN_DX(1,0)+NodePosition[4]*rDN_DX(2,0);
  // Fgrad(0,1)= NodePosition[0]*rDN_DX(0,1)+NodePosition[2]*rDN_DX(1,1)+NodePosition[4]*rDN_DX(2,1);
  // Fgrad(1,0)= NodePosition[1]*rDN_DX(0,0)+NodePosition[3]*rDN_DX(1,0)+NodePosition[5]*rDN_DX(2,0);
  // Fgrad(1,1)= NodePosition[1]*rDN_DX(0,1)+NodePosition[3]*rDN_DX(1,1)+NodePosition[5]*rDN_DX(2,1);

  // //Determinant of the material time derivative of the deformation gradient tensor
  // FJacobian= Fgrad(0,0)*Fgrad(1,1)-Fgrad(0,1)*Fgrad(1,0);
   
  // //Inverse
  // invFgrad.resize(2,2);
  // if(FJacobian==0){
  //   FJacobian=1;
  // }else{
  //   invFgrad(0,0)=  (Fgrad(1,1)/FJacobian);
  //   invFgrad(0,1)= -(Fgrad(0,1)/FJacobian);
  //   invFgrad(1,0)= -(Fgrad(1,0)/FJacobian);
  //   invFgrad(1,1)=  (Fgrad(0,0)/FJacobian); 
  // }

    // std::cout<<"             "<< FJacobian;
    // std::cout<<"__ "<< invFgrad(0,0);
    // std::cout<<"__ "<< invFgrad(1,1);
    // std::cout<<"__ "<< invFgrad(0,1);
    // std::cout<<"__ "<< invFgrad(1,0)<<std::endl;
    // std::cout<<"   :: "<< Fgrad(0,0);
    // std::cout<<" "<< Fgrad(0,1);
    // std::cout<<" "<< Fgrad(1,0);
    // std::cout<<" "<< Fgrad(1,1);


}

  

template < unsigned int TDim> 
void TwoStepUpdatedLagrangianVPElement<TDim>::CalcVelDefGrad(const ShapeFunctionDerivativesType& rDN_DX,
							     MatrixType &FgradVel,
							     MatrixType &invFgradVel,
							     double &FVelJacobian,
							     const double theta)
{
  GeometryType& rGeom = this->GetGeometry();
  const SizeType NumNodes = rGeom.PointsNumber();
  const SizeType LocalSize = TDim*NumNodes;
  VectorType VelocityValues = ZeroVector(LocalSize);
  VectorType RHSVelocities = ZeroVector(LocalSize);
  this->GetVelocityValues(RHSVelocities,0); 
  RHSVelocities*=theta;
  this->GetVelocityValues(VelocityValues,1);
  RHSVelocities+=VelocityValues*(1.0-theta);

  FgradVel.resize(TDim,TDim);

  FgradVel=ZeroMatrix(TDim,TDim);
  for (SizeType i = 0; i < TDim; i++)
    {
      for (SizeType j = 0; j < TDim; j++)
	{
	  for (SizeType k = 0; k < NumNodes; k++)
	    {
	      FgradVel(i,j)+= RHSVelocities[TDim*k+i]*rDN_DX(k,j);
	    }
	}
    }


 //Inverse
  invFgradVel.resize(TDim,TDim);
  invFgradVel=ZeroMatrix(TDim,TDim);
  FVelJacobian=1;

  MathUtils<double>::InvertMatrix( FgradVel, invFgradVel, FVelJacobian );


}


template < > 
void TwoStepUpdatedLagrangianVPElement<2>::CalcVolumetricDefRate(const ShapeFunctionDerivativesType& rDN_DX,
								 double &volumetricDefRate,
								 MatrixType &invGradDef,
								 const double theta)
{


  GeometryType& rGeom = this->GetGeometry();
  const SizeType NumNodes = rGeom.PointsNumber();
  const SizeType LocalSize = 2*NumNodes;
  VectorType VelocityValues = ZeroVector(LocalSize);
  VectorType RHSVelocities = ZeroVector(LocalSize);
  this->GetVelocityValues(RHSVelocities,0); 
  RHSVelocities*=theta;
  this->GetVelocityValues(VelocityValues,1);
  RHSVelocities+=VelocityValues*(1.0-theta);
  // RHSVelocities=UpdatedVelocities;
  // RHSVelocities=CurrentVelocities;

  // double lagDNX0=rDN_DX(0,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(0,1)*rElementalVariables.InvFgrad(1,0);
  // double lagDNX1=rDN_DX(1,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(1,1)*rElementalVariables.InvFgrad(1,0);
  // double lagDNX2=rDN_DX(2,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(2,1)*rElementalVariables.InvFgrad(1,0);
  // double lagDNY0=rDN_DX(0,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(0,1)*rElementalVariables.InvFgrad(1,1);
  // double lagDNY1=rDN_DX(1,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(1,1)*rElementalVariables.InvFgrad(1,1);
  // double lagDNY2=rDN_DX(2,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(2,1)*rElementalVariables.InvFgrad(1,1);
  double lagDNX0=rDN_DX(0,0)*invGradDef(0,0)+rDN_DX(0,1)*invGradDef(1,0);
  double lagDNX1=rDN_DX(1,0)*invGradDef(0,0)+rDN_DX(1,1)*invGradDef(1,0);
  double lagDNX2=rDN_DX(2,0)*invGradDef(0,0)+rDN_DX(2,1)*invGradDef(1,0);
  double lagDNY0=rDN_DX(0,0)*invGradDef(0,1)+rDN_DX(0,1)*invGradDef(1,1);
  double lagDNY1=rDN_DX(1,0)*invGradDef(0,1)+rDN_DX(1,1)*invGradDef(1,1);
  double lagDNY2=rDN_DX(2,0)*invGradDef(0,1)+rDN_DX(2,1)*invGradDef(1,1);


  volumetricDefRate= lagDNX0*RHSVelocities[0] + lagDNX1*RHSVelocities[2] + lagDNX2*RHSVelocities[4];
  volumetricDefRate+=lagDNY0*RHSVelocities[1] + lagDNY1*RHSVelocities[3] + lagDNY2*RHSVelocities[5];

}

template < unsigned int TDim > 
void TwoStepUpdatedLagrangianVPElement<TDim>::CalcSpatialVelocityGrad(MatrixType &invFgrad,
								      MatrixType &VelDefgrad,
								      MatrixType &SpatialVelocityGrad)
{
  SpatialVelocityGrad.resize(TDim,TDim);
  
  SpatialVelocityGrad=prod(VelDefgrad,invFgrad);

  // SpatialVelocityGrad(0,0)=VelDefgrad(0,0)*invFgrad(0,0) + VelDefgrad(0,1)*invFgrad(1,0);
  // SpatialVelocityGrad(0,1)=VelDefgrad(0,0)*invFgrad(0,1) + VelDefgrad(0,1)*invFgrad(1,1);
  // SpatialVelocityGrad(1,0)=VelDefgrad(1,0)*invFgrad(0,0) + VelDefgrad(1,1)*invFgrad(1,0);
  // SpatialVelocityGrad(1,1)=VelDefgrad(1,0)*invFgrad(0,1) + VelDefgrad(1,1)*invFgrad(1,1);

}


template < unsigned int TDim > 
void TwoStepUpdatedLagrangianVPElement<TDim>::CalcVolDefRateFromSpatialVelGrad(double &volumetricDefRate,
									       MatrixType &SpatialVelocityGrad)
{
  volumetricDefRate=0;
  for (SizeType i = 0; i < TDim; i++)
    {
      volumetricDefRate+=SpatialVelocityGrad(i,i);
    }
}



template < > 
void TwoStepUpdatedLagrangianVPElement<2>::CheckStrain1(double &VolumetricDefRate,
							MatrixType &SpatialVelocityGrad)
{
  double trace_l=SpatialVelocityGrad(0,0)+SpatialVelocityGrad(1,1);
  if(fabs(trace_l-VolumetricDefRate)<0.0000001){
  }else{
    std::cout<<" ERROR IN CHECKSTRAIN(1) -> ";
    std::cout<<"trace_l= "<<trace_l<<" VolDefRate"<<VolumetricDefRate<<std::endl;
  }
}

template < > 
void TwoStepUpdatedLagrangianVPElement<2>::CalcMDGreenLagrangeMaterial(MatrixType &Fgrad,
								       MatrixType &VelDefgrad, 
								       VectorType &MDGreenLagrangeMaterial)
{

  // x-component
  MDGreenLagrangeMaterial[0]=VelDefgrad(0,0)*Fgrad(0,0) + VelDefgrad(1,0)*Fgrad(1,0);
  // y-component
  MDGreenLagrangeMaterial[1]=VelDefgrad(1,1)*Fgrad(1,1) + VelDefgrad(0,1)*Fgrad(0,1);
  // xy-component
  MDGreenLagrangeMaterial[2]=(VelDefgrad(0,0)*Fgrad(0,1) + VelDefgrad(1,0)*Fgrad(1,1) +
  			      VelDefgrad(0,1)*Fgrad(0,0) + VelDefgrad(1,1)*Fgrad(1,0))*0.5;
}



template < > 
void TwoStepUpdatedLagrangianVPElement<3>::CalcMDGreenLagrangeMaterial(MatrixType &Fgrad,
								       MatrixType &VelDefgrad, 
								       VectorType &MDGreenLagrangeMaterial)
{
  MatrixType FgradTransp= ZeroMatrix(3,3);
  MatrixType VelDefgradTransp= ZeroMatrix(3,3);
  MatrixType part1= ZeroMatrix(3,3);
  MatrixType part2= ZeroMatrix(3,3);

  FgradTransp=Fgrad;
  FgradTransp(0,1)=Fgrad(1,0);
  FgradTransp(0,2)=Fgrad(2,0);
  FgradTransp(1,0)=Fgrad(0,1);
  FgradTransp(1,2)=Fgrad(2,1);
  FgradTransp(2,0)=Fgrad(0,2);
  FgradTransp(2,1)=Fgrad(1,2);

  VelDefgradTransp=VelDefgrad;
  VelDefgradTransp(0,1)=VelDefgrad(1,0);
  VelDefgradTransp(0,2)=VelDefgrad(2,0);
  VelDefgradTransp(1,0)=VelDefgrad(0,1);
  VelDefgradTransp(1,2)=VelDefgrad(2,1);
  VelDefgradTransp(2,0)=VelDefgrad(0,2);
  VelDefgradTransp(2,1)=VelDefgrad(1,2);

  part1=prod(VelDefgradTransp,Fgrad);
  part2=prod(FgradTransp,VelDefgrad);

  MDGreenLagrangeMaterial[0]= ( part1(0,0) + part2(0,0) ) * 0.5;  //xx-component
  MDGreenLagrangeMaterial[1]= ( part1(1,1) + part2(1,1) ) * 0.5;  //yy-component
  MDGreenLagrangeMaterial[2]= ( part1(2,2) + part2(2,2) ) * 0.5;  //zz-component
  MDGreenLagrangeMaterial[3]= ( part1(0,1) + part2(0,1) ) * 0.5;  //xy-component
  MDGreenLagrangeMaterial[4]= ( part1(0,2) + part2(0,2) ) * 0.5;  //xz-component
  MDGreenLagrangeMaterial[5]= ( part1(1,2) + part2(1,2) ) * 0.5;  //yz-component

}




template < > 
void TwoStepUpdatedLagrangianVPElement<2>::CalcSpatialDefRate(VectorType &MDGreenLagrangeMaterial,
							      MatrixType &invFgrad,
							      VectorType &SpatialDefRate)
{
  // x-component
  SpatialDefRate[0]= invFgrad(0,0)*MDGreenLagrangeMaterial[0]*invFgrad(0,0) + 
    invFgrad(1,0)*MDGreenLagrangeMaterial[2]*invFgrad(0,0)*2 +
    invFgrad(1,0)*MDGreenLagrangeMaterial[1]*invFgrad(1,0);
  // y-component
  SpatialDefRate[1]= invFgrad(0,1)*MDGreenLagrangeMaterial[0]*invFgrad(0,1) + 
    invFgrad(0,1)*MDGreenLagrangeMaterial[2]*invFgrad(1,1)*2 +
    invFgrad(1,1)*MDGreenLagrangeMaterial[1]*invFgrad(1,1);
  // xy-component
  SpatialDefRate[2]=invFgrad(0,0)*MDGreenLagrangeMaterial[0]*invFgrad(0,1) + 
    invFgrad(0,0)*MDGreenLagrangeMaterial[2]*invFgrad(1,1) +
    invFgrad(1,0)*MDGreenLagrangeMaterial[2]*invFgrad(0,1) +
    invFgrad(1,0)*MDGreenLagrangeMaterial[1]*invFgrad(1,1);
}


template < > 
void TwoStepUpdatedLagrangianVPElement<3>::CalcSpatialDefRate(VectorType &MDGreenLagrangeMaterial,
							      MatrixType &invFgrad,
							      VectorType &SpatialDefRate)
{
  MatrixType MDGLM= ZeroMatrix(3,3);
  MatrixType invFgradTransp= ZeroMatrix(3,3);
  MatrixType part1= ZeroMatrix(3,3);
  MatrixType totalMatrix= ZeroMatrix(3,3);

  invFgradTransp=invFgrad;
  invFgradTransp(0,1)=invFgrad(1,0);
  invFgradTransp(0,2)=invFgrad(2,0);
  invFgradTransp(1,0)=invFgrad(0,1);
  invFgradTransp(1,2)=invFgrad(2,1);
  invFgradTransp(2,0)=invFgrad(0,2);
  invFgradTransp(2,1)=invFgrad(1,2);

  MDGLM(0,0)=MDGreenLagrangeMaterial[0];  //XX-component;
  MDGLM(1,1)=MDGreenLagrangeMaterial[1];  //YY-component;
  MDGLM(2,2)=MDGreenLagrangeMaterial[2];  //ZZ-component;
  MDGLM(0,1)=MDGreenLagrangeMaterial[3];  //XY-component;
  MDGLM(1,0)=MDGreenLagrangeMaterial[3];  //XY-component;
  MDGLM(0,2)=MDGreenLagrangeMaterial[4];  //ZX-component;
  MDGLM(2,0)=MDGreenLagrangeMaterial[4];  //ZX-component;
  MDGLM(1,2)=MDGreenLagrangeMaterial[5];  //YZ-component;
  MDGLM(2,1)=MDGreenLagrangeMaterial[5];  //YZ-component;

  part1=prod(MDGLM,invFgrad);

  totalMatrix=prod(invFgradTransp,part1);
 
  SpatialDefRate[0]=totalMatrix(0,0);
  SpatialDefRate[1]=totalMatrix(1,1);
  SpatialDefRate[2]=totalMatrix(2,2);
  SpatialDefRate[3]=totalMatrix(0,1);
  SpatialDefRate[4]=totalMatrix(0,2);
  SpatialDefRate[5]=totalMatrix(1,2);
}



template < > 
void TwoStepUpdatedLagrangianVPElement<2>::CalcDeviatoricInvariant(VectorType &SpatialDefRate,
									double &DeviatoricInvariant)
{
  double trace_d=SpatialDefRate[0]+SpatialDefRate[1];
  double dev_X=SpatialDefRate[0]-trace_d/3.0;
  double dev_Y=SpatialDefRate[1]-trace_d/3.0;
  DeviatoricInvariant=sqrt(2*(dev_X*dev_X+SpatialDefRate[2]*SpatialDefRate[2]+ dev_Y*dev_Y));

}


template < > 
void TwoStepUpdatedLagrangianVPElement<2>::CalcEquivalentStrainRate(VectorType &SpatialDefRate,
								    double &EquivalentStrainRate)
{
  EquivalentStrainRate=sqrt(2.0*(SpatialDefRate[0]*SpatialDefRate[0] +
				 SpatialDefRate[1]*SpatialDefRate[1] +
				 2.0*SpatialDefRate[2]*SpatialDefRate[2]));
}


template < > 
void TwoStepUpdatedLagrangianVPElement<3>::CalcDeviatoricInvariant(VectorType &SpatialDefRate,
								   double &DeviatoricInvariant)
{
  double trace_d=SpatialDefRate[0]+SpatialDefRate[1]+SpatialDefRate[2];
  double dev_X=SpatialDefRate[0]-trace_d/3.0;
  double dev_Y=SpatialDefRate[1]-trace_d/3.0;
  double dev_Z=SpatialDefRate[2]-trace_d/3.0;
  DeviatoricInvariant=sqrt(2*(dev_X*dev_X+dev_Y*dev_Y+dev_Z*dev_Z+
			      SpatialDefRate[3]*SpatialDefRate[3]+
			      SpatialDefRate[4]*SpatialDefRate[4]+
			      SpatialDefRate[5]*SpatialDefRate[5]));
}


template < > 
void TwoStepUpdatedLagrangianVPElement<3>::CalcEquivalentStrainRate(VectorType &SpatialDefRate,
								    double &EquivalentStrainRate)
{
  EquivalentStrainRate=sqrt(2.0*(SpatialDefRate[0]*SpatialDefRate[0] +
				 SpatialDefRate[1]*SpatialDefRate[1] +
				 SpatialDefRate[2]*SpatialDefRate[2] +
				 2.0*SpatialDefRate[3]*SpatialDefRate[3] +
				 2.0*SpatialDefRate[4]*SpatialDefRate[4] +
				 2.0*SpatialDefRate[5]*SpatialDefRate[5]));
}

template < > 
double TwoStepUpdatedLagrangianVPElement<2>::CalcNormalProjectionDefRate(VectorType &SpatialDefRate)
{
  
  double NormalProjSpatialDefRate=0;
  GeometryType& rGeom = this->GetGeometry();
  const SizeType NumNodes = rGeom.PointsNumber();
  array_1d<double, 3>  NormalMean(3,0.0);

  for (SizeType i = 0; i < (NumNodes-1); i++)
      {
        for (SizeType j = (i+1); j < NumNodes; j++)
	  {
	    if(rGeom[i].Is(FREE_SURFACE) && rGeom[j].Is(FREE_SURFACE)){
	      noalias(NormalMean) += (rGeom[i].FastGetSolutionStepValue(NORMAL) +
				      rGeom[j].FastGetSolutionStepValue(NORMAL))*0.5;
	    }
	  }

      }

  // if(rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE)){
  //   NormalMean = (rGeom[0].FastGetSolutionStepValue(NORMAL) + rGeom[1].FastGetSolutionStepValue(NORMAL))*0.5;
  // }else if(rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE)){
  //   NormalMean = (rGeom[0].FastGetSolutionStepValue(NORMAL) + rGeom[2].FastGetSolutionStepValue(NORMAL))*0.5;
  // }else if(rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE)){
  //   NormalMean = (rGeom[1].FastGetSolutionStepValue(NORMAL) + rGeom[2].FastGetSolutionStepValue(NORMAL))*0.5;
  // }

 NormalProjSpatialDefRate=NormalMean[0]*SpatialDefRate[0]*NormalMean[0]+
    NormalMean[1]*SpatialDefRate[1]*NormalMean[1]+
    2*NormalMean[0]*SpatialDefRate[2]*NormalMean[1];

 return NormalProjSpatialDefRate;
}


template < > 
double TwoStepUpdatedLagrangianVPElement<3>::CalcNormalProjectionDefRate(VectorType &SpatialDefRate)
{
  double NormalProjSpatialDefRate=0;
  GeometryType& rGeom = this->GetGeometry();
  const SizeType NumNodes = rGeom.PointsNumber();
  array_1d<double, 3>  NormalMean(3,0.0);

    for (SizeType i = 0; i < (NumNodes-2); i++)
      {
        for (SizeType j = (i+1); j < (NumNodes-1); j++)
	  {
	    for (SizeType k = (j+1); k < NumNodes; k++)
	      {
		if(rGeom[i].Is(FREE_SURFACE) && rGeom[j].Is(FREE_SURFACE) && rGeom[k].Is(FREE_SURFACE)){
		  noalias(NormalMean) += (rGeom[i].FastGetSolutionStepValue(NORMAL) +
					  rGeom[j].FastGetSolutionStepValue(NORMAL) +
					  rGeom[k].FastGetSolutionStepValue(NORMAL)) / 3.0;
		}
	      }
	  }

      }


  // if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)){
  //   NormalMean = (rGeom[0].FastGetSolutionStepValue(NORMAL) + 
  // 		  rGeom[1].FastGetSolutionStepValue(NORMAL) + 
  // 		  rGeom[2].FastGetSolutionStepValue(NORMAL))/3.0;
  // }else if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
  //   NormalMean = (rGeom[0].FastGetSolutionStepValue(NORMAL) +
  // 		  rGeom[1].FastGetSolutionStepValue(NORMAL) +
  // 		  rGeom[3].FastGetSolutionStepValue(NORMAL))/3.0;
  // }else if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
  //   NormalMean = (rGeom[0].FastGetSolutionStepValue(NORMAL) +
  // 		  rGeom[2].FastGetSolutionStepValue(NORMAL) +
  // 		  rGeom[3].FastGetSolutionStepValue(NORMAL))/3.0;
  // }else if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){      
  //   NormalMean = (rGeom[1].FastGetSolutionStepValue(NORMAL) +
  // 		  rGeom[2].FastGetSolutionStepValue(NORMAL) +
  // 		  rGeom[3].FastGetSolutionStepValue(NORMAL))/3.0;
  // }

    NormalProjSpatialDefRate=NormalMean[0]*SpatialDefRate[0]*NormalMean[0]+
      NormalMean[1]*SpatialDefRate[1]*NormalMean[1]+
      NormalMean[2]*SpatialDefRate[2]*NormalMean[2]+
      2*NormalMean[0]*SpatialDefRate[3]*NormalMean[1]+
      2*NormalMean[0]*SpatialDefRate[4]*NormalMean[2]+
      2*NormalMean[1]*SpatialDefRate[5]*NormalMean[2];

    return NormalProjSpatialDefRate;

}


template < > 
void TwoStepUpdatedLagrangianVPElement<2>::CheckStrain2(MatrixType &SpatialVelocityGrad,
							MatrixType &Fgrad,
							MatrixType &VelDefgrad)
{
  if(fabs(VelDefgrad(0,0)-SpatialVelocityGrad(0,0)*Fgrad(0,0)-SpatialVelocityGrad(0,1)*Fgrad(1,0))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(2a) VDg(0,0)="<<VelDefgrad(0,0)<<" SVG: "<<std::endl;
    std::cout<<SpatialVelocityGrad(0,0)<<" "<<SpatialVelocityGrad(0,1)<<" __ ";
    std::cout<<SpatialVelocityGrad(1,0)<<" "<<SpatialVelocityGrad(1,1)<<" __ "<<std::endl;
  }
  if(fabs(VelDefgrad(0,1)-SpatialVelocityGrad(0,0)*Fgrad(0,1)-SpatialVelocityGrad(0,1)*Fgrad(1,1))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(2b) VDg(0,0)="<<VelDefgrad(0,0)<<" SVG: "<<std::endl;
    std::cout<<SpatialVelocityGrad(0,0)<<" "<<SpatialVelocityGrad(0,1)<<" __ ";
    std::cout<<SpatialVelocityGrad(1,0)<<" "<<SpatialVelocityGrad(1,1)<<" __ "<<std::endl;
  }
  if(fabs(VelDefgrad(1,0)-SpatialVelocityGrad(1,0)*Fgrad(0,0)-SpatialVelocityGrad(1,1)*Fgrad(1,0))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(2c) VDg(0,0)="<<VelDefgrad(0,0)<<" SVG: "<<std::endl;
    std::cout<<SpatialVelocityGrad(0,0)<<" "<<SpatialVelocityGrad(0,1)<<" __ ";
    std::cout<<SpatialVelocityGrad(1,0)<<" "<<SpatialVelocityGrad(1,1)<<" __ "<<std::endl;
  }
  if(fabs(VelDefgrad(1,1)-SpatialVelocityGrad(1,0)*Fgrad(0,1)-SpatialVelocityGrad(1,1)*Fgrad(1,1))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(2d) VDg(0,0)="<<VelDefgrad(0,0)<<" SVG: "<<std::endl;
    std::cout<<SpatialVelocityGrad(0,0)<<" "<<SpatialVelocityGrad(0,1)<<" __ ";
    std::cout<<SpatialVelocityGrad(1,0)<<" "<<SpatialVelocityGrad(1,1)<<" __ "<<std::endl;
  }
}


template < > 
bool TwoStepUpdatedLagrangianVPElement<2>::CheckStrain3(VectorType &SpatialDefRate,
							MatrixType &SpatialVelocityGrad)
{
  bool computeElement=true;
  if(fabs(SpatialDefRate[0]-SpatialVelocityGrad(0,0))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(3a) Sdf[0]="<<SpatialDefRate[0]<<" SVG: "<<SpatialVelocityGrad(0,0)<<std::endl;
    computeElement=false;
  }
  if(fabs(SpatialDefRate[1]-SpatialVelocityGrad(1,1))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(3b) Sdf[0]="<<SpatialDefRate[1]<<" SVG: "<<SpatialVelocityGrad(1,1)<<std::endl;
    computeElement=false;
  }
  if(fabs(SpatialDefRate[2]-0.5*(SpatialVelocityGrad(1,0)+SpatialVelocityGrad(0,1)))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(3c) Sdf[0]="<<SpatialDefRate[2]<<" SVG10: "<<SpatialVelocityGrad(1,0)<<" SVG: "<<SpatialVelocityGrad(0,1)<<std::endl;
    computeElement=false;
  }
  return computeElement;
}

template < > 
bool TwoStepUpdatedLagrangianVPElement<3>::CheckStrain3(VectorType &SpatialDefRate,
							MatrixType &SpatialVelocityGrad)
{
  bool computeElement=true;
  if(fabs(SpatialDefRate[0]-SpatialVelocityGrad(0,0))<0.0000001){
   }else{
    std::cout<<"ERROR IN CHECKSTRAIN(3a)";
    computeElement=false;
  }
  if(fabs(SpatialDefRate[1]-SpatialVelocityGrad(1,1))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(3b)";
    computeElement=false;
  }
  if(fabs(SpatialDefRate[2]-SpatialVelocityGrad(2,2))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(3c)";
    computeElement=false;
  }
  if(fabs(SpatialDefRate[3]-0.5*(SpatialVelocityGrad(1,0)+SpatialVelocityGrad(0,1)))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(3d)";
    computeElement=false;
  }
 if(fabs(SpatialDefRate[4]-0.5*(SpatialVelocityGrad(2,0)+SpatialVelocityGrad(0,2)))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(3e)";
    computeElement=false;
  }
 if(fabs(SpatialDefRate[5]-0.5*(SpatialVelocityGrad(2,1)+SpatialVelocityGrad(1,2)))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(3f)";
    computeElement=false;
  }
  return computeElement;
}


template < > 
void TwoStepUpdatedLagrangianVPElement<3>::CalcVolumetricDefRate(const ShapeFunctionDerivativesType& rDN_DX,
								 double &volumetricDefRate,
								 MatrixType &invGradDef,
								 const double theta)
{
  std::cout<<"TO BE IMPLEMENTED ------- CalcVolumetricDefRate -------"<<std::endl;
  //you can compute the volumetric deformation rate using CalcVolDefRateFromSpatialVelGrad
}


template < > 
void TwoStepUpdatedLagrangianVPElement<3>::CheckStrain1(double &VolumetricDefRate,
							     MatrixType &SpatialVelocityGrad)
{
  std::cout<<"TO BE IMPLEMENTED ------- CheckStrain1 -------"<<std::endl;
}


template < > 
void TwoStepUpdatedLagrangianVPElement<3>::CheckStrain2(MatrixType &SpatialVelocityGrad,
							MatrixType &Fgrad,
							MatrixType &VelDefgrad)
{
  std::cout<<"TO BE IMPLEMENTED ------- CheckStrain2 -------"<<std::endl;
}


  template< unsigned int TDim >
  double TwoStepUpdatedLagrangianVPElement<TDim>::EquivalentStrainRate(const ShapeFunctionDerivativesType &rDN_DX) const
  {
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    // Calculate Symetric gradient
    MatrixType S = ZeroMatrix(TDim,TDim);
    for (unsigned int n = 0; n < NumNodes; ++n)
      {
        const array_1d<double,3>& rVel = rGeom[n].FastGetSolutionStepValue(VELOCITY,1); // OLD VELOCITY (which is incompressible, unlike the fractional step one)
        for (unsigned int i = 0; i < TDim; ++i)
	  for (unsigned int j = 0; j < TDim; ++j)
	    S(i,j) += 0.5 * ( rDN_DX(n,j) * rVel[i] + rDN_DX(n,i) * rVel[j] );
      }

    // Norm of symetric gradient
    double NormS = 0.0;
    for (unsigned int i = 0; i < TDim; ++i)
      for (unsigned int j = 0; j < TDim; ++j)
	NormS += S(i,j) * S(i,j);

    return std::sqrt(2.0*NormS);
  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::ComputeLumpedMassMatrix(Matrix& rMassMatrix,
									const double Weight,
									double & MeanValue)
  {

    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    double Count=0;

    if((NumNodes==3 && TDim==2) || (NumNodes==4 && TDim==3)){
      double Coeff=1.0+TDim;
      for (SizeType i = 0; i < NumNodes; ++i)
	{
	  double Mij = Weight/Coeff;

	  for ( unsigned int j = 0; j <  TDim; j++ )
	    {
	      unsigned int index = i * TDim + j;
	      rMassMatrix( index, index ) += Mij;
	      Count+=1.0;
	      MeanValue+=Mij;
	    }

	}
    }
    else if(NumNodes==6 && TDim==2){
      double Mij = Weight/57.0;
      double consistent=1.0;
      for (SizeType i = 0; i < NumNodes; ++i)
	{
	  if(i<3){
	    consistent=3.0;
	  }else{
	    consistent=16.0;
	  }
	  for ( unsigned int j = 0; j <  TDim; j++ )
	    {
	      unsigned int index = i * TDim + j;
	      rMassMatrix( index, index ) += Mij*consistent;
	      Count+=1.0;
	      MeanValue+=Mij;
	    }

	}

    }else{
      std::cout<<"ComputeLumpedMassMatrix 3D quadratic not yet implemented!"<<std::endl;
    }
      MeanValue*=1.0/Count;

  }




  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::ComputeMassMatrix(Matrix& rMassMatrix,
								  const ShapeFunctionsType& rN,
								  const double Weight,
								  double & MeanValue)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    IndexType FirstRow = 0;
    IndexType FirstCol = 0;
    double Count=0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        for (SizeType j = 0; j < NumNodes; ++j)
	  {
            const double Mij = Weight * rN[i] * rN[j];
		
            for (SizeType d =  0; d < TDim; ++d){
	      rMassMatrix(FirstRow+d,FirstCol+d) += Mij;
	      MeanValue+=fabs(Mij);
	      Count+=1.0;
	    }
            FirstCol += TDim;
	  }
        FirstRow += TDim;
        FirstCol = 0;
      }
    MeanValue*=1.0/Count;

  }

 
  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::AddExternalForces(Vector& rRHSVector,
								       const double Density,
								       const ShapeFunctionsType& rN,
								       const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    SizeType FirstRow = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {

	if( this->GetGeometry()[i].SolutionStepsDataHas(VOLUME_ACCELERATION) ){ // it must be checked once at the begining only
	  array_1d<double, 3 >& VolumeAcceleration = this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);
	  // Build RHS
	  for (SizeType d = 0; d < TDim; ++d)
	    {
	      // Body force
	      rRHSVector[FirstRow+d] += Weight * Density * rN[i] * VolumeAcceleration[d];
	    }

	}

        FirstRow += TDim;

      }
  }



  template< >
  void TwoStepUpdatedLagrangianVPElement<2>::AddInternalForces(Vector& rRHSVector,
								    const ShapeFunctionDerivativesType& rDN_DX,
								    ElementalVariables& rElementalVariables,
								    const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    SizeType FirstRow = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
	double lagDNXi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,0);
	double lagDNYi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,1);
	// lagDNXi=rDN_DX(i,0);
	// lagDNYi=rDN_DX(i,1);

	rRHSVector[FirstRow]   += -Weight*(lagDNXi*rElementalVariables.UpdatedTotalCauchyStress[0]+
					   lagDNYi*rElementalVariables.UpdatedTotalCauchyStress[2]);
	
	rRHSVector[FirstRow+1] += -Weight*(lagDNYi*rElementalVariables.UpdatedTotalCauchyStress[1]+
					   lagDNXi*rElementalVariables.UpdatedTotalCauchyStress[2]);


	FirstRow += 2;
      }
  }

  template< >
  void TwoStepUpdatedLagrangianVPElement<3>::AddInternalForces(Vector& rRHSVector,
							       const ShapeFunctionDerivativesType& rDN_DX,
							       ElementalVariables& rElementalVariables,
							       const double Weight)
  {

    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    SizeType FirstRow = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
	double lagDNXi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,0)+rDN_DX(i,2)*rElementalVariables.InvFgrad(2,0);
	double lagDNYi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,1)+rDN_DX(i,2)*rElementalVariables.InvFgrad(2,1);
	double lagDNZi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,2)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,2)+rDN_DX(i,2)*rElementalVariables.InvFgrad(2,2);
	// lagDNXi=rDN_DX(i,0);
	// lagDNYi=rDN_DX(i,1);
	// lagDNZi=rDN_DX(i,2);

	rRHSVector[FirstRow]   += -Weight*(lagDNXi*rElementalVariables.UpdatedTotalCauchyStress[0]+
					   lagDNYi*rElementalVariables.UpdatedTotalCauchyStress[3]+
					   lagDNZi*rElementalVariables.UpdatedTotalCauchyStress[4]);

	rRHSVector[FirstRow+1] += -Weight*(lagDNYi*rElementalVariables.UpdatedTotalCauchyStress[1]+
					   lagDNXi*rElementalVariables.UpdatedTotalCauchyStress[3]+
					   lagDNZi*rElementalVariables.UpdatedTotalCauchyStress[5]);

	rRHSVector[FirstRow+2] += -Weight*(lagDNZi*rElementalVariables.UpdatedTotalCauchyStress[2]+
					   lagDNXi*rElementalVariables.UpdatedTotalCauchyStress[4]+
					   lagDNYi*rElementalVariables.UpdatedTotalCauchyStress[5]);

	FirstRow += 3;
      }


  }


  /*
   * Template class definition (this should allow us to compile the desired template instantiations)
   */

  template class TwoStepUpdatedLagrangianVPElement<2>;
  template class TwoStepUpdatedLagrangianVPElement<3>;

}
