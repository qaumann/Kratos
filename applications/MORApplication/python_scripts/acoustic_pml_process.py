# Importing the Kratos Library
import KratosMultiphysics
# Import MORApplication
import KratosMultiphysics.MORApplication as MOR
from KratosMultiphysics.MORApplication.frequency_dependent_material_process import FrequencyDependentMaterialProcess

import sys

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return AcousticPMLProcess(Model, settings["Parameters"])

class AcousticPMLProcess(KratosMultiphysics.Process):
    def __init__(self, Model, linear_solver, strategy, settings):

        default_settings = KratosMultiphysics.Parameters(
            """
            {
                "help"                           : "This process performs all preprocessing steps for an acoustic PML.",
                "pml_model_part_name"            : "PML",
                "interface_model_part_name"      : "interface",
                "boundary_model_part_name"       : "boundary",
                "stretching_function_type"       : "type",
                "stretching_function_parameters" : {
                    "a" : 0,
                    "m" : 0
                }
            }
            """
        )

        settings.RecursivelyValidateAndAssignDefaults(default_settings)

        KratosMultiphysics.Process.__init__(self)
        self.Model = Model
        self.linear_solver = linear_solver
        self.strategy = strategy
        self.pml_process = []
        self.pml_model_part = Model[settings["pml_model_part_name"].GetString()]
        self.interface_model_part = Model[settings["interface_model_part_name"].GetString()]
        self.boundary_model_part = Model[settings["boundary_model_part_name"].GetString()]
        stretching_function_type = settings["stretching_function_type"].GetString()

        if settings["stretching_function_type"].GetString() == "polynomial":
            #"polynomial" is of the form (a*r**m)/(m*d**(m-1)), r is the distance between point and interface, d the thickness of the PML
            a = settings["stretching_function_parameters"]["a"].GetDouble()
            if a <= 0.0:
                err_msg  = 'Please provide a > 0.0 .'
                raise Exception(err_msg)
            m = settings["stretching_function_parameters"]["m"].GetInt()
            if m < 1:
                err_msg  = 'Please provide an integer m > 0 .'
                raise Exception(err_msg)
            self.stretching_function = lambda r,d: (a*r**m)/(m*d**(m-1))
        else:
            err_msg  = 'Unknown stretching function type "{}". '.format(stretching_function_type)
            err_msg += 'Possible choices are "polynomial".'
            raise Exception(err_msg)

    def ExecuteInitialize(self):
        # check if DOMAIN_SIZE is set. Otherwise no gradient is computed
        if self.pml_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 0:
            err_msg  = 'DOMAIN_SIZE must be defined for the PML modelpart.'
            raise Exception(err_msg)

        # compute direction, distance, and width in the PML layer
        process = MOR.AcousticPMLDirectionProcess(
            self.pml_model_part, \
            self.interface_model_part, \
            self.boundary_model_part, \
            self.linear_solver)

        process.Execute()

        # apply the stretching function
        for n in self.pml_model_part.Nodes:
            d = n.GetValue(MOR.PML_LOCAL_WIDTH)
            r = n.GetValue(MOR.PML_FACTOR)
            val = self.stretching_function(r,d)
            n.SetValue(MOR.PML_FACTOR, val)

        # define the frequency dependent process
        frequency_dependent_process_settings = KratosMultiphysics.Parameters()
        frequency_dependent_process_settings.AddString("model_part_name", self._GetFullModelPartPath())
        frequency_dependent_process_settings.AddString("contribution_type", "pml")
        self.pml_process = FrequencyDependentMaterialProcess(self.Model, self.strategy, frequency_dependent_process_settings)
        self.pml_process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        self.pml_process.ExecuteInitializeSolutionStep()

    def _GetFullModelPartPath(self):
        this_mp = self.pml_model_part
        mp_path = this_mp.Name
        while not this_mp == self.pml_model_part.GetRootModelPart():
            this_mp = this_mp.GetParentModelPart()
            mp_path = this_mp.Name + "." + mp_path

        return mp_path
