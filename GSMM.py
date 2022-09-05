from abc import abstractclassmethod


class GSMM:
    '''
    Genome-scale metabolic model abstract super class
    '''
    def __init__(self) -> None:
        pass

    @abstractclassmethod
    def build_model(self):
        pass

    @abstractclassmethod   
    def run_media_analysis(self):
        pass