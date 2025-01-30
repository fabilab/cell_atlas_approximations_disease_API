class DiseaseNotFoundError(KeyError):
    def __init__(self, msg, disease):
        self.disease = disease
        super().__init__(self, msg)

class CellTypeNotFoundError(KeyError):
    def __init__(self, msg, cell_type):
        self.cell_type = cell_type
        super().__init__(self, msg)
        
class TissueNotFoundError(KeyError):
    def __init__(self, msg, tissue):
        self.tissue = tissue
        super().__init__(self, msg)
        
class FeatureNotFoundError(KeyError):
    def __init__(self, msg, feature):
        self.feature = feature
        super().__init__(self, msg)
        
class SomeFeaturesNotFoundError(KeyError):
    def __init__(self, msg, features):
        self.features = features
        super().__init__(self, msg)

class DevelopmentStageNotFoundError(KeyError):
    def __init__(self, msg, development_stage_general):
        self.development_stage_general = development_stage_general
        super().__init__(self, msg)

class NoMatchingDatasetsError(ValueError):
    def __init__(self, msg, filters):
        self.filters = filters
        super().__init__(msg)

class NoContrastingConditionsInADatasetError(ValueError):
    def __init__(self, msg, filters):
        self.filters = filters
        super().__init__(msg)