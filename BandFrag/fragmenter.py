#!/scistor/tc/dtt741/BandFrag/venv/bin/python3.11

import sys
import subprocess
import os

class LatticeXyz(object):
    atoms = list()
    lattice_vectors = list()

    def __init__(self, atoms, lattice_vectors):
        self.atoms = atoms
        self.lattice_vectors = lattice_vectors

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

class FragmentMapping(object):
    fragment_mapping_1 = None
    fragment_mapping_2 = None

    def __init__(self, fragment_mapping_1, fragment_mapping_2):
        self.fragment_mapping_1 = read_fragment_mapping(fragment_mapping_1)
        self.fragment_mapping_2 = read_fragment_mapping(fragment_mapping_2)

class Configuration(object):
    IRC_path = None
    fragment_mapping: FragmentMapping = None

    strain_fragment_1 = 0
    strain_fragment_2 = 0

    bond_length_atom_1 = 0
    bond_length_atom_2 = 0

    bond_length = 0

    def __init__(self, IRC_path, fragment_mapping, strain_fragment_1, strain_fragment_2, bond_length_atom_1, bond_length_atom_2, bond_length):
        self.IRC_path = IRC_path
        self.fragment_mapping = fragment_mapping
        self.strain_fragment_1 = strain_fragment_1
        self.strain_fragment_2 = strain_fragment_2
        self.bond_length_atom_1 = bond_length_atom_1
        self.bond_length_atom_2 = bond_length_atom_2
        self.bond_length = bond_length

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

def read_parameter(config_list, parameter_name, parameter_count):
    if parameter_name in config_list:
        return_values = []
        for i in range(parameter_count):
            return_values.append(config_list[config_list.index(parameter_name)+1+i].rstrip())
        return return_values if len(return_values) > 1 else return_values[0]
    else:
        return None
    
def get_xyz_list(IRC_path):
    with open(IRC_path, "r") as file:
        data = file.read()
    
    xyz_list: list[LatticeXyz] = list()

    # Every frame starts with this, lets hope they dont change that ;_;
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

def scrape_template_script(data: str):
    # This is always the same and at the top I don't care to scrape it manually
    bin_bash = "#!/bin/bash\n"
    splitter = "End\n"
    splice = data.split(splitter, maxsplit = 1)[1]
    return bin_bash + splice

# Generate the fragments atoms from the mapping
def generate_atoms_from_mapping(fragment_mapping_1: list[int], xyz: LatticeXyz):
    mapped_atoms_string = ""

    for num in fragment_mapping_1:
        num -= 1 # The mapping is 1 indexed, but our list is 0 indexed
        mapped_atoms_string += f"{xyz.atoms[num][0]} {round(xyz.atoms[num][1], 8)} {round(xyz.atoms[num][2], 8)} {round(xyz.atoms[num][3], 8)}\n"
    
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

def run_bash_scripts(bash_scripts: list[str], calculation_folder):
    for i in range(len(bash_scripts)):
        script_folder = f"{calculation_folder}/scripts"
        frame_name = f"frame_{i}"
        script_path = f"{script_folder}/{frame_name}.sh"

        os.makedirs(script_folder, exist_ok = True)
        with open(script_path, "w") as file:
            file.write(bash_scripts[i].replace("FRAME_DIR", frame_name))
        
        os.makedirs(f"{calculation_folder}/{frame_name}", exist_ok = True)

        print(f"Running fragment_{i}.sh")
        
        subprocess.run(["bash", script_path])

data = sys.stdin.read()

calculation_folder = sys.argv[1]

configuration: Configuration = read_configuration(data, calculation_folder)

xyz_list = get_xyz_list(configuration.IRC_path)

bash_scripts = build_bash_scripts(data, configuration, xyz_list)

run_bash_scripts(bash_scripts, calculation_folder)