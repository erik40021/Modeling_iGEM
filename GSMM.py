from abc import abstractclassmethod


class GSMM:
    '''
    Genome-scale metabolic model super class
    '''
    def __init__(self) -> None:
        pass

    @abstractclassmethod
    def build_model(self):
        pass

    @abstractclassmethod   
    def run_medium_analysis(self):
        pass