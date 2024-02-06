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
