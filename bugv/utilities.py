"""
version 0.2.1 by jiajinbu 2021.01.04
"""

def update_dict(origin_config, new_config):
    """
    Use new_config to update origin_config.
    You can use . in the key in new_config, 
    Such as: 
    origin_config = {
        "color" : "red",
        "line" : { "color" : "blue"
        }
    }
    new_config = {
        "color" : "white",
        "line.color" : "red",
        "line.label" : "red line"
    }
    update_dict(origin_config, new_config)
    print(origin_config)
    origin_config = {
        "color" : "white",
        "line" : { "color" : "red",
                   "label" : "red line"
        }
    }
    Note:
    If the value in new_config is list, then will perform shallow copy. If 
    the value in new_config is dict, then will perform shallow copy and iter update the value in dict.
    If the dict is empty, it will update the value to empty.
    """
    
    def get_dict_value_by_list(config, keys):
        v = config
        for k in keys:
            if k not in v:
                v[k] = {}
            v = v[k]
        return(v) 
                  
    for k, v in new_config.items():
        kd = k.split(".")
        if len(kd) > 1:
            pointer_config = get_dict_value_by_list(origin_config, kd[:-1])
            k = kd[-1]
        else:
            pointer_config = origin_config
            k = kd[0]
        if isinstance(v, dict):
            if k not in pointer_config:
                pointer_config[k] = {}
            if v:
                update_dict(pointer_config[k], v)
            else:
                pointer_config[k] = {}
        elif isinstance(v,list):
            pointer_config[k] = v.copy()
        else:
            pointer_config[k] = v
