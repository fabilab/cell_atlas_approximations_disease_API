from flask import request
from flask_restful import abort

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
        except Exception as exc:
            return {
                "message": str(exc),
                "error": {
                    "type": "internal_error",
                    "details": str(exc)
                }
            }, 500  # HTTP 500 Internal Server Error
    return inner