//
// Author: Joaquín Irazábal jirazabal@cimne.upc.edu
//

#if !defined(KRATOS_BEAM_PARTICLE_H_INCLUDED)
#define KRATOS_BEAM_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include "spheric_continuum_particle.h"


namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) BeamParticle : public SphericContinuumParticle {

        public:

        /// Pointer definition of BeamParticle
        KRATOS_CLASS_POINTER_DEFINITION(BeamParticle);

        typedef SphericContinuumParticle BaseType;
        typedef BaseType::ParticleDataBuffer BaseBufferType;
        typedef std::unique_ptr<BaseType::ParticleDataBuffer> BaseBufferPointerType;

        /// Default constructor.
        BeamParticle();
        BeamParticle( IndexType NewId, GeometryType::Pointer pGeometry);
        BeamParticle( IndexType NewId, NodesArrayType const& ThisNodes);
        BeamParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);
        BeamParticle(Element::Pointer p_continuum_spheric_particle);

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

        /// Destructor
        virtual ~BeamParticle(){}

        BeamParticle& operator=(const BeamParticle& rOther);

        /// Turn back information as a string
        std::string Info() const override
        {
            std::stringstream buffer;
            buffer << "BeamParticle" ;
            return buffer.str();
        }

        /// Print information about this object
        void PrintInfo(std::ostream& rOStream) const override {rOStream << "BeamParticle";}

        /// Print object's data
        void PrintData(std::ostream& rOStream) const override {}

        virtual void Initialize(const ProcessInfo& r_process_info) override;

        virtual void InitializeSolutionStep(const ProcessInfo& r_process_info) override;

        virtual void ComputeBallToBallContactForce(SphericParticle::ParticleDataBuffer &,
                                                   const ProcessInfo& r_process_info,
                                                   array_1d<double, 3>& rElasticForce,
                                                   array_1d<double, 3>& rContactForce,
                                                   double& RollingResistance) override;

        void Move(const double delta_t, const bool rotation_option, const double force_reduction_factor, const int StepFlag) override;

        virtual void FinalizeSolutionStep(ProcessInfo& r_process_info) override;

        virtual double GetParticleInitialCohesion();
        void   SetParticleInitialCohesionFromProperties(double* particle_initial_cohesion);
        virtual double GetAmountOfCohesionFromStress();
        void   SetAmountOfCohesionFromStressFromProperties(double* amount_of_cohesion_from_stress);
        virtual double GetParticleConicalDamageContactRadius();
        void   SetParticleConicalDamageContactRadiusFromProperties(double* particle_contact_radius);
        virtual double GetParticleConicalDamageMaxStress();
        void   SetParticleConicalDamageMaxStressFromProperties(double* particle_max_stress);
        virtual double GetParticleConicalDamageGamma();
        void   SetParticleConicalDamageGammaFromProperties(double* particle_gamma);
        virtual double GetLevelOfFouling();
        void   SetLevelOfFoulingFromProperties(double* level_of_fouling);

        double SlowGetParticleInitialCohesion();
        double SlowGetAmountOfCohesionFromStress();
        double SlowGetParticleConicalDamageContactRadius();
        double SlowGetParticleConicalDamageMaxStress();
        double SlowGetParticleConicalDamageGamma();
        double SlowGetLevelOfFouling();

        std::vector<double> mNeighbourContactRadius;
        std::vector<double> mNeighbourRigidContactRadius;
        std::vector<double> mNeighbourIndentation;
        std::vector<double> mNeighbourRigidIndentation;
        std::vector<double> mNeighbourTgOfStatFriAng;
        std::vector<double> mNeighbourTgOfDynFriAng;
        std::vector<double> mNeighbourRigidTgOfStatFriAng;
        std::vector<double> mNeighbourRigidTgOfDynFriAng;
        std::vector<double> mNeighbourContactStress;
        std::vector<double> mNeighbourRigidContactStress;
        std::vector<double> mNeighbourCohesion;
        std::vector<double> mNeighbourRigidCohesion;

        protected:

        class ParticleDataBuffer: public SphericParticle::ParticleDataBuffer
        {
            public:

            ParticleDataBuffer(SphericParticle* p_this_particle): SphericParticle::ParticleDataBuffer(p_this_particle){}

            virtual ~ParticleDataBuffer(){}

            bool SetNextNeighbourOrExit(int& i) override
            {
                while (i < int(mpThisParticle->mNeighbourElements.size()) && (mpThisParticle->mNeighbourElements[i]==NULL)){
                    i++;
                }

                if (i < int(mpThisParticle->mNeighbourElements.size())) {
                    SetCurrentNeighbour(mpThisParticle->mNeighbourElements[i]);
                    mpOtherParticleNode = &(mpOtherParticle->GetGeometry()[0]);
                    return true;
                }

                else { // other_neighbour is nullified upon exiting loop
                    mpOtherParticle = NULL;
                    mpOtherParticleNode = NULL;
                    return false;
                }
            }
        };

        std::unique_ptr<SphericParticle::ParticleDataBuffer> CreateParticleDataBuffer(SphericParticle* p_this_particle) override
        {
            return std::unique_ptr<SphericParticle::ParticleDataBuffer>(new ParticleDataBuffer(p_this_particle));
        }

        void ComputeNewNeighboursHistoricalData(DenseVector<int>& temp_neighbours_ids,
                                                std::vector<array_1d<double, 3> >& temp_neighbour_elastic_contact_forces) override;

        void ComputeNewRigidFaceNeighboursHistoricalData() override;

        private:

        void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericContinuumParticle);
        }

        void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericContinuumParticle);
        }

}; // Class BeamParticle

    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream, BeamParticle& rThis) {return rIStream;}

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, const BeamParticle& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

}  // namespace Kratos

#endif // KRATOS_BEAM_PARTICLE_H_INCLUDED defined
