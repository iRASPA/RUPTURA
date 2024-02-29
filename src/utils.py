isothermMeta = {
    'Langmuir': {
        'labels': ['q_sat', 'b']
    },
    'Anti-Langmuir': {
        'labels': ['a', 'b']
    },
    'BET': {
        'labels': ['q_sat', 'b', 'c']
    },
    'Henry': {
        'labels': ['a']
    },
    'Freundlich': {
        'labels': ['a', 'nu']
    },
    'Sips': {
        'labels': ['q_sat', 'b', 'nu']
    },
    'Langmuir-Freundlich': {
        'labels': ['q_sat', 'b', 'nu']
    },
    'Redlich-Peterson': {
        'labels': ['a', 'b', 'nu']
    },
    'Toth': {
        'labels': ['q_sat', 'b', 'nu']
    },
    'Unilan': {
        'labels': ['q_sat', 'b', 'eta']
    },
    'OBrien-Myers': {
        'labels': ['q_sat', 'b', 'sigma']
    },
    'Quadratic': {
        'labels': ['q_sat', 'b', 'c']
    },
    'Temkin': {
        'labels': ['q_sat', 'b', 'c']
    },
    'Bingel-Walton': {
        'labels': ['q_sat', 'a', 'c']
    }
}

pressureScales = {"log": 0, "linear": 1}
markers = ["o", "+", "^", "D", "x", "*", "p", "s", "v"]
getMarker = lambda i: markers[i % len(markers)]

#################################################################
###### Functions to read input files from simulation.input ######
#################################################################

def is_float(value):
    """Check if the given value can be converted to a float."""
    try:
        float(value)
        return True
    except ValueError:
        return False

def format_value(value):
    """Format the input value based on its type."""
    if isinstance(value, list):
        # Process lists recursively, unpacking single-item lists
        if len(value) == 1:
            return format_value(value[0])
        else:
            return [format_value(item) for item in value]
    elif is_float(value):
        v = float(value)
        return v if v%1.0!=0.0 else int(v)
    elif value.lower() == "yes":
        return True
    elif value.lower() == "no": 
        return False
    else:
        return value
    
def parse_lines(lines):
    """Parse the lines from the configuration file and format them accordingly."""
    formatted_lines = [line.split("//")[0].strip() for line in lines if line.strip() and not line.strip().startswith("//")]
    return [line for line in formatted_lines if line]

def input_to_config(file_path):
    """Process the configuration file and generate the output configuration."""
    output_config = {"components": []}
    comp = None
    simulation_type = None
    
    with open(file_path, 'r') as file:
        lines = file.readlines()

    lines = parse_lines(lines)

    for line in lines:
        sp = line.split()
        if sp[0] == "Component":
            comp = 0 if comp is None else comp + 1
            output_config["components"].append({sp[2]: sp[3], "isotherms": []})
        elif sp[0] == "SimulationType":
            simulation_type = sp[1]
            output_config[simulation_type] = {}
        elif sp[0] in isothermMeta.keys():
            output_config["components"][comp]["isotherms"].append(format_value(sp))
        elif sp[0] == "NumberOfIsothermSites":
            continue
        else:
            if comp is None:
                output_config[simulation_type][sp[0]] = format_value(sp[1:])
            else:
                output_config["components"][comp][sp[0]] = format_value(sp[1:])
    return output_config