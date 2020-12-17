from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.MORApplication as MOR
import KratosMultiphysics.KratosUnittest as KratosUnittest
import math

from KratosMultiphysics.MORApplication.frequency_dependent_material_process import FrequencyDependentMaterialProcess

import KratosMultiphysics.kratos_utilities as kratos_utils
linear_solvers_application_available = kratos_utils.CheckIfApplicationsAvailable("LinearSolversApplication")

class TestFrequencyDependentMaterialProcess(KratosUnittest.TestCase):

    @KratosUnittest.skipUnless(linear_solvers_application_available,"Missing required application: LinearSolversApplication")
    def test_frequency_dependent_biot_material(self):
        import KratosMultiphysics.LinearSolversApplication as LinearSolversApplication
        # create model
        model = KratosMultiphysics.Model()
        mp = model.CreateModelPart("domain")
        mp.GetProperties()[0].SetValue(MOR.DENSITY_FLUID,1.21)
        mp.GetProperties()[0].SetValue(MOR.VISCOSITY_FLUID,1.84e-5)
        mp.GetProperties()[0].SetValue(MOR.STANDARD_PRESSURE_FLUID,101e3)
        mp.GetProperties()[0].SetValue(MOR.HEAT_CAPACITY_FLUID,1.4)
        mp.GetProperties()[0].SetValue(MOR.PRANDTL_NUMBER_FLUID,0.71)
        mp.GetProperties()[0].SetValue(MOR.DENSITY_SOLID,750.0)
        mp.GetProperties()[0].SetValue(MOR.LAMBDA_SOLID,487500.0)
        mp.GetProperties()[0].SetValue(MOR.MUE_SOLID,325000.0)
        mp.GetProperties()[0].SetValue(MOR.DAMPING_SOLID,0.0)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.POROSITY,0.96)
        mp.GetProperties()[0].SetValue(MOR.TORTUOSITY,1.7)
        mp.GetProperties()[0].SetValue(MOR.FLOW_RESISTIVITY,32e3)
        mp.GetProperties()[0].SetValue(MOR.VISCOUS_LENGTH,90e-6)
        mp.GetProperties()[0].SetValue(MOR.THERMAL_LENGTH,165e-6)

        # define solver
        complex_linear_solver = LinearSolversApplication.ComplexSparseLUSolver()
        scheme = MOR.MatrixBuilderScheme()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(KratosMultiphysics.LinearSolver())
        strategy = MOR.FrequencyResponseAnalysisStrategy(
            mp,
            scheme,
            builder_and_solver,
            complex_linear_solver,
            False,
            False)

        from KratosMultiphysics.MORApplication.frequency_dependent_material_process import FrequencyDependentMaterialProcess
        process_settings = KratosMultiphysics.Parameters(
                    """
                    {
                        "model_part_name"   : "domain",
                        "contribution_type" : "biot"
                    }
                    """
                )
        process = FrequencyDependentMaterialProcess(model, strategy, process_settings)
        process.ExecuteInitialize()
        freq = 12
        mp.ProcessInfo[MOR.FREQUENCY] = freq
        process.ExecuteInitializeSolutionStep()

        self.assertAlmostEqual(process.settings['k1'].factor, 0.959999543740704 - 0.000453749531958j)
        self.assertAlmostEqual(process.settings['k2'].factor, 3.770737987069413e-07 + 3.749996131886144e-04j)
        self.assertAlmostEqual(process.settings['m1'].factor, 50.181040501380200 - 0.079061318429012j)
        self.assertAlmostEqual(process.settings['m2'].factor, -1.851188307065158e-09 - 7.767041061366900e-07j)
        self.assertAlmostEqual(process.settings['m3'].factor, 1.382399342986613e+02 - 6.533993260198417e-02j)

if __name__ == '__main__':
    KratosUnittest.main()
