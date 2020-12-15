# Importing the Kratos Library
import KratosMultiphysics
# Import MORApplication
import KratosMultiphysics.MORApplication as MOR

import sys
from cmath import sqrt

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return FrequencyDependentMaterialProcess(Model, settings["Parameters"])

class FrequencyDependentMaterialProcess(KratosMultiphysics.Process):
    def __init__(self, Model, strategy, settings):

        default_settings = KratosMultiphysics.Parameters(
            """
            {
                "help"              : "This process adds and updates frequency dependent material properties.",
                "model_part_name"   : "Structure",
                "contribution_type" : "biot"
            }
            """
        )

        settings.ValidateAndAssignDefaults(default_settings)

        KratosMultiphysics.Process.__init__(self)
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.contribution_type = settings["contribution_type"].GetString()
        self.strategy = strategy
        self.settings = {}
        self.functions = {}

    def ExecuteInitialize(self):
        if self.contribution_type == "biot":
            self._SetUpBiot()
            pass
        elif self.contribution_type == "ghm":
            pass
        else:
            err_msg  = 'Unknown contribution type "{}".'.format(contribution_type)
            err_msg += 'Possible choices are "biot" and "ghm".'
            raise Exception(err_msg)

    def ExecuteInitializeSolutionStep(self):
        frequency = self.model_part.ProcessInfo[MOR.FREQUENCY]
        for key in self.settings.keys():
            self.settings[key].factor = self.functions[key](frequency)
            print(self.functions[key](frequency))

    def _SetUpBiot(self):
        lambda_solid = self._RetrieveProperty(MOR.LAMBDA_SOLID)
        damping_solid = self._RetrieveProperty(MOR.DAMPING_SOLID)
        density_solid = self._RetrieveProperty(MOR.DENSITY_SOLID)
        density_fluid = self._RetrieveProperty(MOR.DENSITY_FLUID)
        viscosity_fluid = self._RetrieveProperty(MOR.VISCOSITY_FLUID)
        standard_pressure_fluid = self._RetrieveProperty(MOR.STANDARD_PRESSURE_FLUID)
        heat_capacity_fluid = self._RetrieveProperty(MOR.HEAT_CAPACITY_FLUID)
        prandtl_number_fluid = self._RetrieveProperty(MOR.PRANDTL_NUMBER_FLUID)
        porosity = self._RetrieveProperty(KratosMultiphysics.POROSITY)
        thermal_length = self._RetrieveProperty(MOR.THERMAL_LENGTH)
        tortuosity = self._RetrieveProperty(MOR.TORTUOSITY)
        flow_resistivity = self._RetrieveProperty(MOR.FLOW_RESISTIVITY)
        viscous_length = self._RetrieveProperty(MOR.VISCOUS_LENGTH)

        viscous_drag = lambda omega: flow_resistivity * porosity**2 * \
            sqrt(1 + 4j*omega*tortuosity**2*viscosity_fluid*density_fluid / \
            (flow_resistivity**2 * viscous_length**2 * porosity**2))
        apparent_mass_density = porosity * density_fluid * (tortuosity-1)

        self.settings['k1'] = MOR.FrequencyDependentMaterialSettings(self.model_part, 91)
        self.functions['k1'] = \
            lambda omega: -porosity * (-apparent_mass_density + 1j*viscous_drag(omega)/omega) / \
            (porosity*density_fluid + apparent_mass_density - 1j*viscous_drag(omega)/omega)
        self.strategy.SetFrequencyDependentMaterial(self.settings['k1'])


        #omega = 12
        #k1 = 0.959999543740704 - 0.000453749531958i
        #k2 = 3.770737987069413e-07 + 3.749996131886144e-04i
        #m1 = 0.348479447926251 - 0.000549036933535i
        #m2 = -1.285547435461915e-11 - 5.393778514838125e-09i

        # VISCOUS_DRAG = FLOW_RESISTIVITY * POROSITY^2 * (1 + ...
        #     4 * 1j * omega * TORTUOSITY^2 * VISCOSITY_FLUID * DENSITY_FLUID / (FLOW_RESISTIVITY^2 * VISCOUS_LENGHT^2 * POROSITY^2))^0.5;
        # APPARENT_MASS_DENSITY = POROSITY * DENSITY_FLUID * (TORTUOSITY -1)

        # alpha = (1 + 8 * VISCOSITY_FLUID/(1j * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID)* ...
        #     (1 + 1i * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID /(16 * VISCOSITY_FLUID))^0.5);

        # k1=-POROSITY*(-APPARENT_MASS_DENSITY+1i*VISCOUS_DRAG/omega)/(POROSITY*DENSITY_FLUID+APPARENT_MASS_DENSITY-1i*VISCOUS_DRAG/omega);
        # k2=(POROSITY^2)/(POROSITY*DENSITY_FLUID+APPARENT_MASS_DENSITY-1i*VISCOUS_DRAG/omega);

        # m1=-1i*VISCOUS_DRAG/omega-(-APPARENT_MASS_DENSITY+1i*VISCOUS_DRAG/omega)^2/(POROSITY*DENSITY_FLUID+APPARENT_MASS_DENSITY-1i*VISCOUS_DRAG/omega);
        # m2=-(HEAT_CAPACITY_FLUID-1)*alpha^(-1)/(HEAT_CAPACITY_FLUID*STANDARD_PRESSURE_FLUID);

    def _RetrieveProperty(self, property_type):
        for prop in self.model_part.GetProperties():
                if prop.Has(property_type):
                    if prop.GetValue(property_type) > sys.float_info.epsilon:
                        value = prop.GetValue(property_type)
                        return value
        return 0.0
