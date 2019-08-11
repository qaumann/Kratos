from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ReplacePropertiesProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ReplacePropertiesProcess(KratosMultiphysics.Process):
    """This process replaces the properties in a given instant

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings ):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "help"                        : "This process replaces the properties in a given instant",
            "materials_filename"          : "",
            "interval"                    : [0.0, 1e30]
        }
        """)

        # Detect "End" as a tag and replace it by a large number
        if settings.Has("interval"):
            if settings["interval"][1].IsString():
                if settings["interval"][1].GetString() == "End":
                    settings["interval"][1].SetDouble(1e30) # = default_settings["interval"][1]
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:"+settings["interval"].PrettyPrintJsonString())

        settings.ValidateAndAssignDefaults(default_settings)

        # The model
        self.model = Model

        # Assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        # Materials filename
        self.materials_filename = settings["materials_filename"].GetString()

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        # We activate/deactivate conditions dependeding of interval
        if self.interval.IsInInterval(current_time):
            if self.materials_filename != "":
                # Add constitutive laws and material properties from json file to model parts.
                material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
                material_settings["Parameters"]["materials_filename"].SetString(self.materials_filename)
                KratosMultiphysics.ReadMaterialsUtility(material_settings, self.model)
