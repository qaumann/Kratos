# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics import kratos_utilities
default_data_comm = KM.DataCommunicator.GetDefault()
if default_data_comm.IsDistributed():
    from KratosMultiphysics.mpi import GatherModelPartUtility

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_io import CoSimulationIO

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication as KratosCoSim
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities

# Other imports
import os

def Create(settings, model, solver_name):
    return EmpireIO(settings, model, solver_name)

communication_folder = ".EmpireIO" # hardcoded in C++

class EmpireIO(CoSimulationIO):
    """IO for the legacy EMPIRE_API
    """
    def __init__(self, settings, model, solver_name):
        super().__init__(settings, model, solver_name)
        self.rank = default_data_comm.Rank()
        self.is_distributed = default_data_comm.Rank()
        if self.rank == 0:
            # Note: calling "EMPIRE_API_Connect" is NOT necessary, it is replaced by the next two lines
            KratosCoSim.EMPIRE_API.EMPIRE_API_SetEchoLevel(self.echo_level)
            KratosCoSim.EMPIRE_API.EMPIRE_API_PrintTiming(self.settings["api_print_timing"].GetBool())

            # delete and recreate communication folder to avoid leftover files
            kratos_utilities.DeleteDirectoryIfExisting(communication_folder)
            os.mkdir(communication_folder)

        if self.is_distributed:
            # creating auxiliar Model for gathering ModelParts on rank 0
            self.aux_model = KM.Model()
            self.gather_utilities = {}

    def Finalize(self):
        kratos_utilities.DeleteDirectoryIfExisting(communication_folder)

    def __del__(self):
        # make sure no communication files are left even if simulation is terminated prematurely
        if os.path.isdir(communication_folder):
            kratos_utilities.DeleteDirectoryIfExisting(communication_folder)
            if self.echo_level > 0:
                cs_tools.cs_print_info(self._ClassName(), "Deleting communication folder in destructor")

    def ImportCouplingInterface(self, interface_config):
        model_part_name = interface_config["model_part_name"]
        comm_name = interface_config["comm_name"]

        if not self.model.HasModelPart(model_part_name):
            main_model_part_name, *sub_model_part_names = model_part_name.split(".")
            model_part_utilities.RecursiveCreateModelParts(self.model[main_model_part_name], ".".join(sub_model_part_names))

        model_part = self.model[model_part_name]
        if self.rank == 0:
            KratosCoSim.EMPIRE_API.EMPIRE_API_recvMesh(model_part, comm_name)

    def ExportCouplingInterface(self, interface_config):
        model_part_name = interface_config["model_part_name"]
        comm_name = interface_config["comm_name"]
        model_part_to_export = self.model[model_part_name]

        if model_part_to_export.IsDistributed():
            gathered_model_part = self.aux_model.CreateModelPart(model_part_to_export.Name+"_gather_mesh_export")
            GatherModelPartUtility(0, model_part_to_export, gathered_model_part)
            model_part_for_api = gathered_model_part
        else:
            model_part_for_api = model_part_to_export

        if self.rank == 0:
            KratosCoSim.EMPIRE_API.EMPIRE_API_sendMesh(model_part_for_api, comm_name)

        # cleanup, the gathered ModelPart is no longer needed
        if model_part_to_export.IsDistributed():
            self.aux_model.DeleteModelPart(model_part_to_export.Name+"_gather_mesh_export")

    def ImportData(self, data_config):
        data_type = data_config["type"]

        if data_type == "coupling_interface_data":
            interface_data = data_config["interface_data"]
            model_part = self.__GetModelPartForAPICalls(interface_data)

            if self.rank == 0:
                KratosCoSim.EMPIRE_API.EMPIRE_API_recvDataField(model_part, self.solver_name+"_"+interface_data.name, interface_data.variable)

            self.__ScatterData(interface_data)
        else:
            raise NotImplementedError('Importing interface data of type "{}" is not implemented for this IO: "{}"'.format(data_type, self._ClassName()))

    def ExportData(self, data_config):
        data_type = data_config["type"]

        if data_type == "coupling_interface_data":
            interface_data = data_config["interface_data"]

            self.__GatherData(interface_data)
            model_part = self.__GetModelPartForAPICalls(interface_data)

            if self.rank == 0:
                KratosCoSim.EMPIRE_API.EMPIRE_API_sendDataField(model_part, self.solver_name+"_"+interface_data.name, interface_data.variable)

        elif data_type == "convergence_signal":
            if self.rank == 0:
                KratosCoSim.EMPIRE_API.EMPIRE_API_sendConvergenceSignal(data_config["is_converged"], self.solver_name)
        else:
            raise NotImplementedError('Exporting interface data of type "{}" is not implemented for this IO: "{}"'.format(data_type, self._ClassName()))

    def PrintInfo(self):
        print("This is the EMPIRE-IO")

    def Check(self):
        pass

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "api_print_timing" : false
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())

        return this_defaults

    def __GetModelPartForAPICalls(self, interface_data):
        model_part = interface_data.GetModelPart()
        if not model_part.IsDistributed():
            return model_part

        aux_mp_name = self.__GetNameOfAuxModelPartFromInterfaceData(interface_data)
        if self.aux_model.HasModelPart(aux_mp_name):
            return self.aux_model[aux_mp_name]
        else:
            cs_tools.cs_print_info(self._ClassName(), 'Creating gathered modelpart for "{}"'.format(model_part.FullName()))
            gather_mp = self.aux_model.CreateModelPart(aux_mp_name)
            if interface_data.location != "node_historical":
                raise Exception("Currently only nodal historical values are supported!")
            gather_mp.AddNodalSolutionStepVariable(interface_data.variable)
            self.gather_utilities[aux_mp_name] = GatherModelPartUtility(0, model_part, gather_mp)
            return gather_mp

    def __ScatterData(self, interface_data):
        model_part = interface_data.GetModelPart()
        if not model_part.IsDistributed():
            # do nothing if ModelPart is not distributed
            return

        aux_mp_name = self.__GetNameOfAuxModelPartFromInterfaceData(interface_data)
        self.gather_utilities[aux_mp_name].ScatterFromMaster(interface_data.variable)


    def __GatherData(self, interface_data):
        model_part = interface_data.GetModelPart()
        if not model_part.IsDistributed():
            # do nothing if ModelPart is not distributed
            return

        aux_mp_name = self.__GetNameOfAuxModelPartFromInterfaceData(interface_data)
        self.gather_utilities[aux_mp_name].GatherOnMaster(interface_data.variable)

    def __GetNameOfAuxModelPartFromInterfaceData(self, interface_data):
        return "aux_" + interface_data.model_part_name.replace(".", "-") + "_var_" + interface_data.variable.Name()
