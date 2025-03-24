from flask import request
from flask_restful import abort

from models.exceptions import (
    DiseaseNotFoundError,
    CellTypeNotFoundError,
    TissueNotFoundError,
    FeatureNotFoundError,
    SomeFeaturesNotFoundError,
    DevelopmentStageNotFoundError,
    NoMatchingDatasetsError,
    NoContrastingConditionsInADatasetError,
    ParamsConflictError,
    UniqueIdNotFoundError,
    InvalidParameterError
)

class FeatureStringFormatError(Exception):
    def __init__(self, msg, features):
        self.features = features
        super().__init__(self, msg)

def required_parameters(*required_args):
    """Function that checks if mandatory parameters are present."""
    
    def inner(wrapped):
        def func(*args_inner, **kwargs_inner):
            
            if request.method == 'POST':
                params = request.args.to_dict()
                params.update(request.form.to_dict())
                params.update(request.get_json(silent=True) or {})
            else:
                params = request.args

            missing_params = [
                arg for arg in required_args 
                if not params.get(arg)
            ]
            
            if missing_params:
                # Use return with the abort call
                return {
                    "message": f"Missing required parameter(s): {', '.join(missing_params)}. Please include them in the request.",
                    "error": {
                        "type": "missing_parameter",
                        "missing_parameters": missing_params
                    }
                }, 400  # HTTP 400 Bad Request
                
            return wrapped(*args_inner, **kwargs_inner)
        return func
    return inner

def model_exceptions(func):
    """Closure that deals with model exceptions at the API level."""
    
    def inner(*args_inner, **kwargs_inner):
        try:
            return func(*args_inner, **kwargs_inner)
        
        except DiseaseNotFoundError as exc:
            abort(
                400,
                message=f"No disease found that matches '{exc.disease}'. Please check the spelling or try a different query.",
                error={
                    "type": "invalid_parameter",
                    "invalid_parameter": "disease",
                    "invalid_value": exc.disease,
                },
            )
        except CellTypeNotFoundError as exc:
            abort(
                400,
                message=f"No cell type found that matches '{exc.cell_type}'.",
                error={
                    "type": "invalid_parameter",
                    "invalid_parameter": "cell_type",
                    "invalid_value": exc.cell_type,
                },
            )
        except TissueNotFoundError as exc:
            abort(
                400,
                message=f"No tissue found that matches '{exc.tissue}'.",
                error={
                    "type": "invalid_parameter",
                    "invalid_parameter": "tissue",
                    "invalid_value": exc.tissue,
                },
            )

        except FeatureNotFoundError as exc:
            abort(
                400,
                message=f"No feature found that matches '{exc.feature}'. Please ensure the feature name is correct.",
                error={
                    "type": "invalid_parameter",
                    "invalid_parameter": "feature",
                    "invalid_value": exc.feature,
                },
            )

        except SomeFeaturesNotFoundError as exc:
            abort(
                400,
                message=f"Some features could not be found: {', '.join(exc.features)}. Please remove it from the query",
                error={
                    "type": "invalid_parameter",
                    "invalid_parameter": "features",
                    "invalid_value": exc.features,
                },
            )
        except DevelopmentStageNotFoundError as exc:
            abort(
                400,
                message=f"No development stage found that matches '{exc.development_stage_general}'.",
                error={
                    "type": "invalid_parameter",
                    "invalid_parameter": "development_stage_general",
                    "invalid_value": exc.development_stage_general,
                },
            )
        
        except NoMatchingDatasetsError as exc:
            abort(
                404,
                message=exc.args[0],
                error={
                    "type": "no_results",
                    "filters_used": exc.filters
                },
            )
        
        except NoContrastingConditionsInADatasetError as exc:
            abort(
                404,
                message=exc.args[0],
                error={
                    "type": "no_results",
                    "filters_used": exc.filters
                },
            )

        except ParamsConflictError:
            abort(
                400,
                message="You can specify either unique_ids or metadata filters, not both.",
                error={
                    "type": "invalid_parameter_combination",
                    "invalid_parameters": ["unique_ids", "metadata filters"],
                },
            )

        except UniqueIdNotFoundError as exc:
            abort(
                400,
                message=f"No unique id found that matches '{exc.unique_ids}'.",
                error={
                    "type": "invalid_parameter",
                    "invalid_parameter": "unique_ids",
                    "invalid_value": exc.unique_ids,
                },
            )
        except InvalidParameterError as exc:
            abort(
                400,
                message=f"Invalid parameter name(s): {', '.join(exc.invalid_param_names)}. "
                "Please check the spelling or refer to the documentation for valid parameter names.",
                error={
                    "type": "invalid_parameter_name",
                    "invalid_param_names": exc.invalid_param_names,
                },
            )
    return inner