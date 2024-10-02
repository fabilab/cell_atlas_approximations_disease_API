"""Baseline values for differential analysis."""

_differential_baselines = {
    "disease": "normal",
    "sex": "female",
    "age": "adult",
}


def get_differential_baseline(differential_axis):
    """Get the baseline value for a given differential axis."""
    return _differential_baselines[differential_axis]
