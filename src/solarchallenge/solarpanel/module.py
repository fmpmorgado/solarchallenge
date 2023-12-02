import pandas as pd

#TODO
PATH_MODULE_DATA = "src/database/CEC_modules.csv"

class SolarModule:

    def __init__(self, parameters : dict):
        self.parameters = parameters

    def list_modules(path : str = PATH_MODULE_DATA) -> None:
        df = pd.read_csv(path)

    @classmethod
    def from_database(cls, name : str, path : str = PATH_MODULE_DATA):
        df = pd.read_csv(path)
        
        try:
            parameters = df[df['Name'] == name].iloc[0].to_dict()
        except IndexError:
            raise IndexError("The module name is not on the database.")

        return cls(parameters)
