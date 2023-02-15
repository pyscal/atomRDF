import json
import yaml
from pyscal_rdf.encoder import NumpyArrayEncoder

def write_file(outfile, data):
    """
    Write a given dict as json file
    
    Parameters
    ----------
    outfile: string
        name of output file. `.json` will be added to the given file name
    
    data: dict
        input data dict
    
    Returns
    -------
    None
    """
    with open(".".join([outfile, "json"]), "w") as fout:
        json.dump(data, fout, cls=NumpyArrayEncoder)
    #with open(".".join([outfile, "yaml"]), "w") as fout:
    #    yaml.unsafe_dump(convert_to_dict(sys), fout)