#!/scistor/tc/dtt741/BandFrag/venv/bin/python3.11
import sys
import subprocess
import os
from pathlib import Path
import re
import math
import pandas as pd

# The fragmenter has 3 parts:
# 1. Create bash files
# We read the configuration in the BandFrag input file, read the molecules we should fragment from the .ams.amv mxyz movie file, and then we copy
# the ams run file from the BandFrag input and copy that for every structure in the movie file and replace placeholders with our molecule, taking into account the BandFrag config
# You can track this in the scripts/ folder that is generated when it runs
# 2. Run bash files
# All the bash files we created from the input, we now just run in order
# After every completed frame we move folders and results around to keep things organized (technically this is done in the bash script)
# 3. Scrape and print results
# Now that everything has been calculated and organized, just read it all out and put it in a .ipynb for graphs and a .csv for raw results

graph_template_name = "template.ipynb"
graph_template_name_destination = "graphs.ipynb"
raw_results_name_destination = "results.csv"

our_folder = Path(os.path.dirname(os.path.realpath(__file__)))

# XYZ object with lattices for periodic structures
class LatticeXyz(object):
    atoms = list()
    lattice_vectors = list()

    def __init__(self, atoms, lattice_vectors):
        self.atoms = atoms
        self.lattice_vectors = lattice_vectors

class Results(object):
    Pauli:float = 0
    dispersion:float = 0
    elstat:float = 0
    orbital:float = 0
    interaction:float = 0

    fragment_1_energy:float = 0
    fragment_2_energy:float = 0

    def __init__(self, Pauli, dispersion, elstat, orbital, interaction, fragment_1_energy, fragment_2_energy):
        self.Pauli = float(Pauli)
        self.dispersion = float(dispersion)
        self.elstat = float(elstat)
        self.orbital = float(orbital)
        self.interaction = float(interaction)

        self.fragment_1_energy = float(fragment_1_energy)
        self.fragment_2_energy = float(fragment_2_energy)

# Read the fragment mapping, which can be 1,2,3-10, 17-20, 37 and convert it to a list of integers
def read_fragment_mapping(fragment_mapping_string):
    fragment_mapping = list()

    for string in fragment_mapping_string.split(","):
        if "-" in string:
            start, end = string.split("-")
            start = int(start)
            end = int(end)
            for i in range(start, end + 1):
                fragment_mapping.append(i)
        else:
            number = int(string)
            fragment_mapping.append(number)
    return fragment_mapping

# Holds which atoms belong to which fragment
class FragmentMapping(object):
    fragment_mapping_1: list[int] = None
    fragment_mapping_2: list[int] = None

    def __init__(self, fragment_mapping_1, fragment_mapping_2):
        self.fragment_mapping_1 = read_fragment_mapping(fragment_mapping_1)
        self.fragment_mapping_2 = read_fragment_mapping(fragment_mapping_2)

# Holds the info we obtain from the 'BandFrag' header of the input file
class Configuration(object):
    # Absolute path that leads to the IRC file
    IRC_path:str = None
    fragment_mapping: FragmentMapping = None

    strain_fragment_1:float = 0
    strain_fragment_2:float = 0

    bond_length_atom_1:int = 0
    bond_length_atom_2:int = 0

    # This is really only for making a nice graph
    bond_length:float = 0

    def __init__(self, IRC_path, fragment_mapping, strain_fragment_1, strain_fragment_2, bond_length_atom_1, bond_length_atom_2, bond_length):
        self.IRC_path = IRC_path
        self.fragment_mapping = fragment_mapping
        self.strain_fragment_1 = float(strain_fragment_1)
        self.strain_fragment_2 = float(strain_fragment_2)
        # These indexes are 1 based, and python is 0 based, so just take care of it here
        self.bond_length_atom_1 = int(bond_length_atom_1) - 1
        self.bond_length_atom_2 = int(bond_length_atom_2) - 1
        self.bond_length = float(bond_length)

# Read the the configuration header from the input file and put it in an object for easy access
def read_configuration(data, calculation_folder):
    start = "BandFrag"
    end = "End"
    splice = (data.split(start))[1].split(end)[0]
    config_list = splice.split(" ")

    IRC_path = read_parameter(config_list, "IRC_path", 1)
    fragment_mapping_1 = read_parameter(config_list, "fragment_mapping_1", 1)
    fragment_mapping_2 = read_parameter(config_list, "fragment_mapping_2", 1)
    
    strain_fragment_1 = read_parameter(config_list, "strain_fragment_1", 1)
    strain_fragment_2 = read_parameter(config_list, "strain_fragment_2", 1)

    bond_length_atom_1, bond_length_atom_2, bond_length = read_parameter(config_list, "bond_length", 3)

    fragment_mapping = FragmentMapping(fragment_mapping_1, fragment_mapping_2)

    return Configuration(calculation_folder + "/" + IRC_path, fragment_mapping, strain_fragment_1, strain_fragment_2, bond_length_atom_1, bond_length_atom_2, bond_length)

# Read a parameter from the config
def read_parameter(config_list, parameter_name, parameter_count):
    if parameter_name in config_list:
        return_values = []
        for i in range(parameter_count):
            return_values.append(config_list[config_list.index(parameter_name) + 1 + i].rstrip())
        return return_values if len(return_values) > 1 else return_values[0]
    else:
        return None

# Read an .ams.amv file with multiple .xyz frames (or one, I cant stop you), which is typically obtained
# from an IRC 
def get_xyz_list(IRC_path):
    with open(IRC_path, "r") as file:
        data = file.read()
    
    xyz_list: list[LatticeXyz] = list()

    # Every frame starts with this, lets hope they dont change that ;_;
    # TODO: find some other way to do this, I hate this
    header = "Geometry"
    data_list = data.split("\n")
    frames = list()

    # Split it into frames, each frame starts with "Geometry"
    for data in data_list:
        if header in data:
            frames.append(data_list.index(data))

    # Now go through the frames we found and extract the lattice vectors and atom coordinates
    for frame_index in frames:
        lattice_vectors = {}
        atom_coordinates = []

        for line in data_list[frame_index:]:
            arguments = line.split(" ")
            # Empty line, we hit. Go to next frame, we must
            # TODO: so if there's no empty line at the end of the file we just ignore it ughh
            if len(arguments) < 2:
                xyz_list.append(LatticeXyz(atom_coordinates, lattice_vectors))
                break
            elif arguments[0] == "VEC1":
                # VEC1 2.0 1.0 0.0
                lattice_vectors["VEC1"] = [float(arguments[1]), float(arguments[2]), float(arguments[3])]
            elif arguments[0] == "VEC2":
                # VEC2 1.0 2.0 0.0
                lattice_vectors["VEC2"] = [float(arguments[1]), float(arguments[2]), float(arguments[3])]
            elif len(arguments) == 4:
                # C 1.0 2.0 1.0
                atom_coordinates.append([arguments[0], float(arguments[1]), float(arguments[2]), float(arguments[3])])

    return xyz_list

# Get the ams script part of the input file (which we get from the <<eor> ... eor in bash)
def scrape_template_script(data: str):
    # This is always the same and at the top I don't care to scrape it manually
    bin_bash = "#!/bin/bash\n"
    splitter = "End\n"
    splice = data.split(splitter, maxsplit = 1)[1]
    return bin_bash + splice

# Generate the fragments atoms from the mapping
def generate_atoms_from_mapping(fragment_mapping_1: list[int], xyz: LatticeXyz):
    mapped_atoms_string = ""
    # this makes ams happy 
    decimals = 8

    for num in fragment_mapping_1:
        num -= 1 # The mapping is 1 indexed, but our list is 0 indexed
        mapped_atoms_string += f"{xyz.atoms[num][0]} {round(xyz.atoms[num][1], decimals)} {round(xyz.atoms[num][2], decimals)} {round(xyz.atoms[num][3], decimals)}\n"
    
    return mapped_atoms_string.rstrip()

# Generate the atom mapping (naming is confusing, fuck)
# But this is the one where it just says '1 72' where 1 is the index in our fragment and 72 the index
# inside the full fragment
def generate_atom_mappings(fragment_mapping: FragmentMapping, xyz: LatticeXyz):

    atom_mapping_1, atom_mapping_2 = "",""

    for i in range(len(xyz.atoms) + 1):
        if i in fragment_mapping.fragment_mapping_1:
            atom_mapping_1 += f"{fragment_mapping.fragment_mapping_1.index(i) + 1} {i}\n"
        if i in fragment_mapping.fragment_mapping_2:
            atom_mapping_2 += f"{fragment_mapping.fragment_mapping_2.index(i) + 1} {i}\n"

    return atom_mapping_1.rstrip(), atom_mapping_2.rstrip()

# Get the data input, clean it into a nice template bash file, and generate a bash file for every .xyz coordinate
def build_bash_scripts(data, configuration: Configuration, xyz_list: list[LatticeXyz]):
    bash_scripts = list()

    template_script: str = scrape_template_script(data)

    for i in range(len(xyz_list)):
        bash_script = template_script

        # Lattices
        lattices = (
        f"{str(xyz_list[i].lattice_vectors['VEC1'][0])} {xyz_list[i].lattice_vectors['VEC1'][1]} {xyz_list[i].lattice_vectors['VEC1'][2]}\n"
        f"{str(xyz_list[i].lattice_vectors['VEC2'][0])} {xyz_list[i].lattice_vectors['VEC2'][1]} {xyz_list[i].lattice_vectors['VEC2'][2]}"
        )
        # I'm not sure why I decided to even make these different
        bash_script = bash_script.replace("FRAGMENT_1_LATTICES", lattices)
        bash_script = bash_script.replace("FRAGMENT_2_LATTICES", lattices)
        bash_script = bash_script.replace("FRAGMENT_FULL_LATTICE", lattices)

        # Fragment_1_xyz
        atom_mapping_fragment_1 = generate_atoms_from_mapping(configuration.fragment_mapping.fragment_mapping_1, xyz_list[i])
        bash_script = bash_script.replace("FRAGMENT_1_XYZ", atom_mapping_fragment_1)

        # Fragment_2_xyz
        atom_mapping_fragment_2 = generate_atoms_from_mapping(configuration.fragment_mapping.fragment_mapping_2, xyz_list[i])
        bash_script = bash_script.replace("FRAGMENT_2_XYZ", atom_mapping_fragment_2)

        # Fragment_full_xyz
        atom_mapping_full = generate_atoms_from_mapping(configuration.fragment_mapping.fragment_mapping_1 + configuration.fragment_mapping.fragment_mapping_2, xyz_list[i])
        bash_script = bash_script.replace("FRAGMENT_FULL_XYZ", atom_mapping_full)

        # Atoms mapping
        atom_mapping_1, atom_mapping_2 = generate_atom_mappings(configuration.fragment_mapping, xyz_list[i])
        bash_script = bash_script.replace("FRAGMENT_1_ATOM_MAPPING", atom_mapping_1)
        bash_script = bash_script.replace("FRAGMENT_2_ATOM_MAPPING", atom_mapping_2)

        # This is a bit sketchy, but 'eor' delimits what to send to .py, but is also
        # the easiest way to controle ams job input, so we use a fake eor when writing the bash
        # file (eoj), and replace that later when we actually want to run ams jobs
        # TODO: dont fucking do this
        bash_script = bash_script.replace("eoj", "eor")

        bash_scripts.append(bash_script)
    return bash_scripts

# Run all the bash scripts we made in .ams (and be so fucking careful you only use this when you're slurming, 
# or you absolutely murder a login/interactive node)
def run_bash_scripts(bash_scripts: list[str], calculation_folder):
    frame_folders = list()

    for i in range(len(bash_scripts)):
        script_folder = f"{calculation_folder}/scripts"
        frame_name = f"frame_{i}"
        script_path = f"{script_folder}/{frame_name}.sh"

        os.makedirs(script_folder, exist_ok = True)
        with open(script_path, "w") as file:
            file.write(bash_scripts[i].replace("FRAME_DIR", frame_name))
        
        frame_folder = f"{calculation_folder}/{frame_name}"
        os.makedirs(frame_folder, exist_ok = True)
        frame_folders.append(frame_folder)

        print(f"Running fragment_{i}.sh")
        
        # Make sure to comment out if you're testing on local
        # It's probably fine because the .ams installation is hard to access by accident without slurm
        # But it still makes me reeeeallly nervous
        subprocess.run(["bash", script_path])
    
    return frame_folders

def scrape_results(frame_folders: list[str]) -> list[Results]:
    results = list()
    for path in frame_folders:
        with open(f"{path}/fragment_full.out", 'r') as file_full:
            data_full = file_full.read()

            # The stuff we need looks something like this
            # Let's hope they never even slightly change how they format this...
            # term      au      eV     kcal/mol kj/mol
            # E_Pauli | 2.05330 55.873 1288.47 5391.0
            # E_disp | -0.00437 -0.119 -2.74 -11.5
            # E_elstat | -0.43276 -11.776 -271.56 -1136.2
            # E_orb | -2.36459 -64.344 -1483.80 -6208.2
            # E_int | -0.74841 -20.365 -469.64 -1965.0
            # yes of course I did this completely myself
            pattern = r'^\s*(E_(?:Pauli|disp|elstat|orb|int))\s*\|\s*(?:\S+\s+){2}(-?\d+(?:\.\d+)?)\s+\S+\s*$'

            energies = dict(re.findall(pattern, data_full, re.MULTILINE))
            fragment_energies = list()

            # Grab the fragment energies
            for log_file in ["fragment_1.log", "fragment_2.log"]:
                with open(f"{path}/{log_file}", 'r') as file_fragment:
                    data_fragment = file_fragment.read()

                    m = list(re.finditer(r"ENERGY OF FORMATION", data_fragment))[-1]
                    subtext = data_fragment[m.end():]
                    if len(subtext) > 0:
                        m2 = re.search(r"([-+]?[0-9]*\.?[0-9]+)\s*KCAL/MOL", subtext, flags=re.IGNORECASE)
                        fragment_energies.append(float(m2.group(1)))

            result:Results = Results(energies["E_Pauli"], 
                                    energies["E_disp"], 
                                    energies["E_elstat"], 
                                    energies["E_orb"], 
                                    energies["E_int"], 
                                    round(fragment_energies[0] - configuration.strain_fragment_1, 4), 
                                    round(fragment_energies[1] - configuration.strain_fragment_2, 4),
                                )
            results.append(result)

    return results

def print_results(results: list[Results], configuration:Configuration, xyz_list:list[LatticeXyz], destination:str):
    distances = list()
    for result in results:
        index = results.index(result)
        xyz_lattice = xyz_list[index]

        # these are lists of ['C', 0.11, 2.22, 1.11]
        atom_1 = xyz_lattice.atoms[configuration.bond_length_atom_1]
        atom_2 = xyz_lattice.atoms[configuration.bond_length_atom_2]
        # sqrt((x2 - x1)^2 + (y2 - y1)^2 ...) - TS bond length
        distances.append(round(math.sqrt((atom_2[1] - atom_1[1])**2 + (atom_2[2] - atom_1[2])**2 + (atom_2[3] - atom_1[3])**2) - configuration.bond_length, 4))


    write_to_notebook(results, configuration, xyz_list, destination, distances)
    write_to_output(results, configuration, xyz_list, destination, distances)

def write_to_notebook(results: list[Results], configuration:Configuration, xyz_list:list[LatticeXyz], destination:str, distances):
    template = ""
    for path in our_folder.glob(graph_template_name):
        with open(path, 'r') as template:
            template = template.read()
        break

    template = template.replace("X_AXIS", str(distances))
    template = template.replace("Y_AXIS_PAULI", str([result.Pauli for result in results]))
    template = template.replace("Y_AXIS_DISPERSION", str([result.dispersion for result in results]))
    template = template.replace("Y_AXIS_ELSTAT", str([result.elstat for result in results]))
    template = template.replace("Y_AXIS_ORBITAL", str([result.orbital for result in results]))
    template = template.replace("Y_AXIS_INTERACTION", str([result.interaction for result in results]))

    template = template.replace("Y_AXIS_FRAGMENT_1", str([result.fragment_1_energy for result in results]))
    template = template.replace("Y_AXIS_FRAGMENT_2", str([result.fragment_2_energy for result in results]))

    # So we can get a nice "X-Y" label on the graph
    template = template.replace("LABEL_BONDS", f"{xyz_list[0].atoms[configuration.bond_length_atom_1][0]}—{xyz_list[0].atoms[configuration.bond_length_atom_2][0]}")

    with open(f"{destination}/{graph_template_name_destination}", 'w') as jupyter_notebook:
        jupyter_notebook.write(template)
        print("Written results to", graph_template_name_destination)

def write_to_output(results: list[Results], configuration:Configuration, xyz_list:list[LatticeXyz], destination:str, distances):
    info = list()

    for result in results:
        info.append([
            results.index(result),
            distances[results.index(result)],
            result.fragment_1_energy,
            result.fragment_2_energy,
            result.Pauli,
            result.dispersion,
            result.elstat,
            result.orbital,
            result.interaction,
        ])


    table = pd.DataFrame(info, columns = [
        "FRAME", 
        "DISTANCE",
        "STRAIN_FRAGMENT_1", 
        "STRAIN_FRAGMENT_2", 
        "PAULI", 
        "DISP", 
        "ELSTAT", 
        "ORBITAL", 
        "INTERACTION",
        ])

    pd.set_option('display.max_colwidth', None)
    pd.set_option('display.max_columns', None)

    table.to_csv(f"{destination}/{raw_results_name_destination}")
    print("Written to", raw_results_name_destination)

    print(table)

### The main code execution part

# Grab all the data from the <<eor> ... eor in the input file
data = sys.stdin.read()

# Get the calculation/slurm submit folder (the one where the input file is in, which slurm gives to our input file)
calculation_folder = sys.argv[1]

# Read the BandFrag configuration header so we can define how to behave
configuration: Configuration = read_configuration(data, calculation_folder)

# Turn the .amv into a list of xyz_lattice objects
xyz_list = get_xyz_list(configuration.IRC_path)

# Build the bash scripts from the input file, config and the list of xyz_lattice frames we need to do
bash_scripts = build_bash_scripts(data, configuration, xyz_list)

# Run the scripts we constructed in ams
frame_folders = run_bash_scripts(bash_scripts, calculation_folder)

# frame_folders = [f"/scistor/tc/dtt741/BandFrags/chlorobenzene/Pd_Pyrr_4NDG/Backward/frame_{i}" for i in range(0, 40 +1)]

results = scrape_results(frame_folders)

# Print the results to something so a human can view it (jupy? print? excel?)
print_results(results, configuration, xyz_list, calculation_folder)