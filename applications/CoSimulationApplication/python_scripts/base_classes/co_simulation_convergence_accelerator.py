from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics as KM

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

class CoSimulationConvergenceAccelerator(object):
    def __init__(self, settings, solver_wrapper):
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(self._GetDefaultSettings())

        self.interface_data = solver_wrapper.GetInterfaceData(self.settings["data_name"].GetString())

        self.echo_level = self.settings["echo_level"].GetInt()

    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def InitializeCouplingIteration(self):
        # Saving the previous data for the computation of the residual
        # and the computation of the solution update
        self.input_data = self.interface_data.GetData()

    def FinalizeCouplingIteration(self):
        pass

    def ComputeAndApplyUpdate(self):
        current_data = self.interface_data.GetData()
        residual = current_data - self.input_data
        updated_data = self.input_data + self.ComputeUpdate(residual, self.input_data)
        self.interface_data.SetData(updated_data)

    def PrintInfo(self):
        '''Function to print Info abt the Object
        Can be overridden in derived classes to print more information
        '''
        cs_tools.cs_print_info("Convergence Accelerator", colors.bold(self._Name()))

    def Check(self):
        print("ConvAcc does not yet implement Check")

    def ComputeUpdate( self, residual, previous_data ):
        raise NotImplementedError('"ComputeUpdate" has to be implemented in the derived class!')

    def _Name(self):
        return self.__class__.__name__

    @classmethod
    def _GetDefaultSettings(cls):
        return KM.Parameters("""{
            "type"       : "UNSPECIFIED",
            "solver"     : "UNSPECIFIED",
            "data_name"  : "UNSPECIFIED",
            "echo_level" : 0
        }""")
