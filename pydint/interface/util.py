"""
reads info
"""

def combine_params(param1, param2, rule='default'):
    """ perform a combining rule for two parameters
    """

    if rule == 'default':
        combined_param = (param1 + param2) / 2.0
    else:
        raise NotImplementedError

    return combined_param
