"""Main module for API v1"""

from api.v1.endpoints import get_api_endpoint
from api.v1.objects import (
    Metadata,
    DifferentialCellTypeAbundance,
    DifferentialGeneExpression,
    HighestMeasurement,
    Average,
    FractionDetected,
    DotPlot,
    # CombinedDifferentialGeneExpression
)

__all__ = ("api_dict",)

# Connect routes to api objects
api_dict = {
    "endpoint_handler": get_api_endpoint,
    "objects": {
        "metadata": Metadata,
        "differential_gene_expression": DifferentialGeneExpression,
        "differential_cell_type_abundance": DifferentialCellTypeAbundance,
        "highest_measurement": HighestMeasurement,
        "average": Average,
        "fraction_detected": FractionDetected,
        "dotplot": DotPlot,
        # "combined_differential": CombinedDifferentialGeneExpression
    },
}

