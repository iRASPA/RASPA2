"""
WRAPSA2, the python wrapper for RASPA2.

This is intended 1. to allow users to automate + simplify their workflows and
2. to enable scaling of simulations to millions of structures.
"""
from ctypes import cdll, c_void_p, c_char_p, c_bool, cast
import os
from textwrap import dedent
from multiprocessing import Process, Pipe
import argparse


from .__init__ import __version__
from RASPA2.output_parser import parse

raspa_dir = os.path.dirname(os.path.realpath(__file__))
libraspa_dir = os.path.join(raspa_dir, "simulations/lib")
libraspa_file = next(f for f in os.listdir(libraspa_dir) if "libraspa" in f)

try:
    import pybel
    PYBEL_LOADED = True
except ImportError:
    PYBEL_LOADED = False


def run(structure, molecule_name, temperature=273.15, pressure=101325,
        helium_void_fraction=1.0, unit_cells=(1, 1, 1),
        framework_name="streamed", simulation_type="MonteCarlo", cycles=2000,
        init_cycles="auto", forcefield="CrystalGenerator",
        input_file_type="cif"):
    """Runs a simulation with the specified parameters.

    Args:
        structure: The structure to test for adsorption, as a string of type
            `input_file_type` (default is "cif").
        molecule_name: The molecule to test for adsorption. A file of the same
            name must exist in `$RASPA_DIR/share/raspa/molecules/TraPPE`.
        temperature: (Optional) The temperature of the simulation, in Kelvin.
        pressure: (Optional) The pressure of the simulation, in Pascals.
        helium_void_fraction: (Optional) The helium void fraction of the input
            structure. Required for excess adsorption back-calculation.
        unit_cells: (Optional) The number of unit cells to use, by dimension.
        framework_name: (Optional) If not streaming, this will load the
            structure at `$RASPA_DIR/share/raspa/structures`. Ignored if
            streaming.
        simulation_type: (Optional) The type of simulation to run, defaults
            to "MonteCarlo".
        cycles: (Optional) The number of simulation cycles to run.
        init_cycles: (Optional) The number of initialization cycles to run.
            Defaults to the minimum of cycles / 2 and 10,000.
        forcefield: (Optional) The forcefield to use. Name must match a folder
            in `$RASPA_DIR/share/raspa/forcefield`, which contains the properly
            named `.def` files.
        input_file_type: (Optional) The type of input structure. Assumes cif.
    Returns:
        A string representing the contents of a simulation input file.

    The goal of this function is to mask the complexity of RASPA by limiting
    parameters and assuming sensible defaults. This should streamline common-
    case usage, but it means that this function won't work in all use cases.
    In these cases, look into loading your own simulation input file and
    passing it to `RASPA.run_script`.
    """
    return parse(run_script(create_script(**locals()), structure))


def run_script(input_script, structure=None, stream=True):
    """Runs RASPA on the inputted structure, returning simulation results.

    Args:
        input_script: A RASPA simulation script, as a loaded string.
        structure: (Optional) Data encoding a CIF, MOL, or CSSR file. If not
            specified, the program will look for the file specified by
            "FrameworkName" in `$RASPA_DIR/share/raspa/structures`.
        stream: (Optional) Reads arguments from strings and writes simulation
            results to stdout. NOTE: Only functional on Linux systems.
    Returns:
        A string representing the output data. (TODO This should be parsed
        and returned as a mutable data structure)

    If a structure is not specified, then nothing will be streamed and this
    will output files and folders. Otherwise, this function will return the
    output of RASPA, as a string.
    """
    # This supports `pybel.Molecule` objects with charge data by converting
    # them into RASPA-formatted cif strings
    if PYBEL_LOADED and isinstance(structure, pybel.Molecule):
        structure = pybel_to_raspa_cif(structure)
    # This supports python objects grabbed from json files or databases
    elif isinstance(structure, dict):
        structure = pybel_to_raspa_cif(json_to_pybel(structure))

    # RASPA leaks a lot of memory, which is an issue in high-throughput
    # screening. Calling each simulation in a subprocess adds complexity and
    # time, but allows still reachable memory to be garbage collected by the
    # OS and also protects against segfaults destablizing an entire engine.
    parent_conn, child_conn = Pipe()
    p = Process(target=_script_subprocess, args=(input_script, structure or "",
                                                 raspa_dir, stream,
                                                 child_conn))
    p.start()
    if stream:
        output = parent_conn.recv()
    p.join(1)
    p.terminate()

    if stream:
        return output


def _script_subprocess(input_script, structure, raspa_dir, stream, conn):
    """Loads and runs libraspa2. Called through multiprocessing.Process.

    Intended as a workaround for RASPA's memory leaks. Once RASPA is running
    valgrind pure, remove this approach and replace with a direct call.
    """
    libraspa = cdll.LoadLibrary(os.path.join(libraspa_dir, libraspa_file))
    libraspa.run.argtypes = (c_char_p, c_char_p, c_char_p, c_bool)
    libraspa.run.restype = c_void_p

    ptr = libraspa.run(input_script.encode("ascii"),
                       structure.encode("ascii"),
                       raspa_dir.encode("ascii"), stream)
    if stream:
        conn.send(cast(ptr, c_char_p).value[:].decode("utf-8"))
    conn.close()


def create_script(molecule_name, temperature=273.15, pressure=101325,
                  helium_void_fraction=1.0, unit_cells=(1, 1, 1),
                  simulation_type="MonteCarlo", cycles=2000,
                  init_cycles="auto", forcefield="CrystalGenerator",
                  input_file_type="cif", **kwargs):
    """Creates a RASPA simulation input file from parameters.

    Args:
        molecule_name: The molecule to test for adsorption. A file of the same
            name must exist in `$RASPA_DIR/share/raspa/molecules/TraPPE`.
        temperature: (Optional) The temperature of the simulation, in Kelvin.
        pressure: (Optional) The pressure of the simulation, in Pascals.
        helium_void_fraction: (Optional) The helium void fraction of the input
            structure. Required for excess adsorption back-calculation.
        unit_cells: (Optional) The number of unit cells to use, by dimension.
        simulation_type: (Optional) The type of simulation to run, defaults
            to "MonteCarlo".
        cycles: (Optional) The number of simulation cycles to run.
        init_cycles: (Optional) The number of initialization cycles to run.
            Defaults to the minimum of cycles / 2 and 10,000.
        forcefield: (Optional) The forcefield to use. Name must match a folder
            in `$RASPA_DIR/share/raspa/forcefield`, which contains the properly
            named `.def` files.
        input_file_type: (Optional) The type of input structure. Assumes cif.
        charged: (Optional) A boolean indicating whether or not to use Ewald
            charge parameters.
    Returns:
        A string representing the contents of a simulation input file.

    The goal of this function is to mask the complexity of RASPA by limiting
    parameters and assuming sensible defaults. This should streamline common-
    case usage, but it means that this function won't work in all use cases.
    In these cases, look into loading your own simulation input file and
    passing it to `RASPA.run_script`.
    """
    is_mol = "yes" if input_file_type.lower() == "mol" else "no"
    print_every = cycles // 10
    a, b, c = unit_cells
    if init_cycles == "auto":
        init_cycles = min(cycles // 2, 10000)

    return dedent("""
                  SimulationType                {simulation_type}
                  NumberOfCycles                {cycles}
                  NumberOfInitializationCycles  {init_cycles}
                  PrintEvery                    {print_every}
                  RestartFile                   no

                  Forcefield                    {forcefield}
                  CutOff                        12.8
                  ChargeMethod                  Ewald
                  EwaldPrecision                1e-6
                  UseChargesFromMOLFile         {is_mol}

                  Framework                     0
                  FrameworkName                 streamed
                  InputFileType                 {input_file_type}
                  UnitCells                     {a} {b} {c}
                  HeliumVoidFraction            {helium_void_fraction}
                  ExternalTemperature           {temperature}
                  ExternalPressure              {pressure}

                  Movies                        no
                  WriteMoviesEvery              100

                  Component 0 MoleculeName             {molecule_name}
                              StartingBead             0
                              MoleculeDefinition       TraPPE
                              IdealGasRosenbluthWeight 1.0
                              TranslationProbability   1.0
                              RotationProbability      1.0
                              ReinsertionProbability   1.0
                              SwapProbability          1.0
                              CreateNumberOfMolecules  0
                  """.format(**locals())).strip()


def run_mixture(structure, molecules, mol_fractions, temperature=273.15,
                pressure=101325, helium_void_fraction=1.0,
                unit_cells=(1, 1, 1), simulation_type="MonteCarlo",
                cycles=2000, init_cycles="auto", forcefield="CrystalGenerator",
                input_file_type="cif"):
    """Runs a simulation with mixture of gases.

    Args:
        structure: The structure to test for adsorption, as a string of type
            `input_file_type` (default is "cif").
        molecules: The molecules to test for adsorption. Files of the same
            names must exist in `$RASPA_DIR/share/raspa/molecules/TraPPE`.
        mol_fractions: The mol fractions of each gas that you want to separate.
            Corresponds to the `molecules` list.
        temperature: (Optional) The temperature of the simulation, in Kelvin.
        pressure: (Optional) The pressure of the simulation, in Pascals.
        helium_void_fraction: (Optional) The helium void fraction of the input
            structure. Required for excess adsorption back-calculation.
        unit_cells: (Optional) The number of unit cells to use, by dimension.
        simulation_type: (Optional) The type of simulation to run, defaults
            to "MonteCarlo".
        cycles: (Optional) The number of simulation cycles to run.
        init_cycles: (Optional) The number of initialization cycles to run.
            Defaults to the minimum of cycles / 2 and 10,000.
        forcefield: (Optional) The forcefield to use. Name must match a folder
            in `$RASPA_DIR/share/raspa/forcefield`, which contains the properly
            named `.def` files.
        input_file_type: (Optional) The type of input structure. Assumes cif.
    Returns:
        A string representing the contents of a simulation input file.

    The goal of this function is to mask the complexity of RASPA by limiting
    parameters and assuming sensible defaults. This should streamline common-
    case usage, but it means that this function won't work in all use cases.
    In these cases, look into loading your own simulation input file and
    passing it to `RASPA.run_script`.
    """
    is_mol = "yes" if input_file_type.lower() == "mol" else "no"
    print_every = cycles // 10
    a, b, c = unit_cells
    molecule_list = " ".join(str(n) for n in range(len(molecules)))
    molecule_count = len(molecules)
    if init_cycles == "auto":
        init_cycles = min(cycles // 2, 10000)

    script = dedent("""
                    SimulationType                {simulation_type}
                    NumberOfCycles                {cycles}
                    NumberOfInitializationCycles  {init_cycles}
                    PrintEvery                    {print_every}
                    RestartFile                   no

                    Forcefield                    {forcefield}
                    CutOff                        12.8
                    ChargeMethod                  Ewald
                    EwaldPrecision                1e-6
                    UseChargesFromMOLFile         {is_mol}

                    Framework                     0
                    FrameworkName                 streamed
                    InputFileType                 {input_file_type}
                    UnitCells                     {a} {b} {c}
                    HeliumVoidFraction            {helium_void_fraction}
                    ExternalTemperature           {temperature}
                    ExternalPressure              {pressure}

                    Movies                        no
                    WriteMoviesEvery              100
                    """.format(**locals())).strip()

    for i, (molecule, fraction) in enumerate(zip(molecules, mol_fractions)):
        script += dedent("""
                      Component {i} MoleculeName                 {molecule}
                                  StartingBead                 0
                                  MoleculeDefinition           TraPPE
                                  MolFraction                  {fraction}
                                  IdentityChangeProbability    1.0
                                  NumberOfIdentityChanges      {molecule_count}
                                  IdentityChangesList          {molecule_list}
                                  IdealGasRosenbluthWeight     1.0
                                  TranslationProbability       1.0
                                  RotationProbability          1.0
                                  ReinsertionProbability       1.0
                                  SwapProbability              1.0
                                  CreateNumberOfMolecules      0
                         """.format(**locals()))
    return parse(run_script(script, structure))


def get_geometric_surface_area(structure, unit_cells=(1, 1, 1), cycles=500,
                               input_file_type="cif", units="m^2/g",
                               forcefield="CrystalGenerator"):
    """Calculates the geometric surface area of an inputted structure.

    Args:
        structure: The structure to use, as a string of type
            `input_file_type` (default is "cif").
        input_file_type: (Optional) The type of input structure. Assumes cif.
        unit_cells: (Optional) The number of unit cells to use, by dimension.
        cycles: (Optional) The number of simulation cycles to run.
        units: (Optional) The units in which to return the surface area. Can be
            "m^2/g", "A^2", or "m^2/cm^3".
        forcefield: (Optional) The forcefield to use. Name must match a folder
            in `$RASPA_DIR/share/raspa/forcefield`, which contains the properly
            named `.def` files.
    Returns:
        The geometric surface area, as a float.
    """
    print_every = cycles // 10
    a, b, c = unit_cells

    script = dedent("""
                    SimulationType                MonteCarlo
                    NumberOfCycles                {cycles}
                    PrintEvery                    {print_every}
                    PrintPropertiesEvery          {print_every}

                    Forcefield                    {forcefield}
                    CutOff                        12.8

                    Framework                     0
                    FrameworkName                 streamed
                    InputFileType                 {input_file_type}
                    UnitCells                     {a} {b} {c}
                    SurfaceAreaProbeDistance      Sigma

                    Component 0 MoleculeName             N2
                                StartingBead             0
                                MoleculeDefinition       TraPPE
                                SurfaceAreaProbability   1.0
                                CreateNumberOfMolecules  0
                    """.format(**locals())).strip()

    output = parse(run_script(script, structure))
    return output["Average Surface Area"]["[{}]".format(units)][0]


def get_helium_void_fraction(structure, unit_cells=(1, 1, 1), cycles=2000,
                             input_file_type="cif",
                             forcefield="CrystalGenerator"):
    """Calculates the helium void fraction of the inputted structure.

    Args:
        structure: The structure to test for helium void fraction,
            as a string of type 'input_file_type` (default is "cif").
        unit_cells: (Optional) The number of unit cells to use, by dimension.
        cycles: (Optional) The number of simulation cycles to run.
        input_file_type: (Optional) The type of input structure. Assumes cif.
        forcefield: (Optional) The forcefield to use. Name must match a folder
            in `$RASPA_DIR/share/raspa/forcefield`, which contains the properly
            named `.def` files.
    Returns:
        The helium void fraction of the structure, as a float.
    """
    print_every = cycles // 10
    a, b, c = unit_cells

    script = dedent("""
                    SimulationType                MonteCarlo
                    NumberOfCycles                {cycles}
                    PrintEvery                    {print_every}
                    PrintPropertiesEvery          {print_every}

                    Forcefield                    {forcefield}
                    CutOff                        12.8

                    Framework                     0
                    FrameworkName                 streamed
                    InputFileType                 {input_file_type}
                    UnitCells                     {a} {b} {c}
                    ExternalTemperature           298.0

                    Component 0 MoleculeName             helium
                                MoleculeDefinition       TraPPE
                                WidomProbability         1.0
                                CreateNumberOfMolecules  0
                    """.format(**locals())).strip()

    output = parse(run_script(script, structure))
    return output["Average Widom Rosenbluth factor"]["Widom"][0]


def get_pore_size_distribution(structure, unit_cells=(1, 1, 1), cycles=500,
                               input_file_type="cif",
                               forcefield="CrystalGenerator",
                               bins=50):
    """Calculates the pore size distribution of the inputted structure.

    Args:
        structure: The structure to test for helium void fraction,
            as a string of type 'input_file_type` (default is "cif").
        unit_cells: (Optional) The number of unit cells to use, by dimension.
        cycles: (Optional) The number of simulation cycles to run.
        input_file_type: (Optional) The type of input structure. Assumes cif.
        forcefield: (Optional) The forcefield to use. Name must match a folder
            in `$RASPA_DIR/share/raspa/forcefield`, which contains the properly
            named `.def` files.
        bins: (Optional) The number of bins to use in the output histogram.
    Returns:
        A list of the form [[x1, x2, ..., xn], [y1, y2, ..., yn]], where x is
        the binned pore size (in Angstroms) and y is the partial pore volume
        (in cm^3 / g).
    """
    print_every = cycles // 10
    a, b, c = unit_cells

    script = dedent("""
                    SimulationType                PSD
                    NumberOfCycles                {cycles}
                    NumberOfInitializationCycles  0
                    PrintEvery                    {print_every}

                    PSDProbeDistance              Sigma
                    WritePSDHistogramEvery        {print_every}
                    PSDHistogramSize              {bins}
                    PSDRange                      12.5

                    Forcefield                    {forcefield}
                    CutOff                        12.8

                    Framework                     0
                    FrameworkName                 streamed
                    InputFileType                 {input_file_type}
                    UnitCells                     {a} {b} {c}
                    ExternalTemperature           100.0
                    """.format(**locals())).strip()

    output = run_script(script, structure)

    # PSD returns a different data structure which must be parsed separately.
    # Luckily, it's an easy parse. Skip # commented lines, split by whitespace.
    info = [[float(x) for x in r.split()] for r in output.strip().splitlines()
            if r[0] != "#"]
    return [[i[0] for i in info], [i[2] for i in info]]


def get_density(molecule_name, temperature=273.15, pressure=101325,
                cycles=5000, init_cycles="auto",
                forcefield="CrystalGenerator"):
    """Calculates the density of a gas through an NPT ensemble.

    Args:
        molecule_name: The molecule to test for adsorption. A file of the same
            name must exist in `$RASPA_DIR/share/raspa/molecules/TraPPE`.
        temperature: (Optional) The temperature of the simulation, in Kelvin.
        pressure: (Optional) The pressure of the simulation, in Pascals.
        cycles: (Optional) The number of simulation cycles to run.
        init_cycles: (Optional) The number of initialization cycles to run.
            Defaults to the minimum of cycles / 2 and 10,000.
        forcefield: (Optional) The forcefield to use. Name must match a folder
            in `$RASPA_DIR/share/raspa/forcefield`, which contains the properly
            named `.def` files.
    Returns:
        The density, as a float, in kg/m^3.
    """
    print_every = cycles // 10
    if init_cycles == "auto":
        init_cycles = min(cycles // 2, 10000)

    script = dedent("""
                    SimulationType                {simulation_type}
                    NumberOfCycles                {cycles}
                    NumberOfInitializationCycles  {init_cycles}
                    PrintEvery                    {print_every}

                    Forcefield                    {forcefield}

                    Box                           0
                    BoxLengths                    30 30 30
                    ExternalTemperature           {temperature}
                    ExternalPressure              {pressure}

                    VolumeChangeProbability       0.25

                    Component 0 MoleculeName             {molecule_name}
                                MoleculeDefinition       TraPPE
                                TranslationProbability   0.5
                                ReinsertionProbability   0.5
                                CreateNumberOfMolecules  256
                  """.format(**locals())).strip()

    output = parse(run_script(script))
    return output["Average Density"]["[kg/m^3]"][0]


def json_to_pybel(data):
    """Converts a python data structure to pybel.Molecule.

    The data structure is a plain python object of the form:

    ```
    {
        "atoms": [{"location": [0, 0, 0], "element": "H", "label": "H1",
                   "charge": 0}, ...],
        "bonds": [{"source": 0, "target": 0, "order": 1}, ...],
        "unitcell": [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    }
    ```

    It is referred to as "json" because this data structure is intended to
    be read from and written to json files or databases.

    As RASPA makes no use of bond information, this field is ignored.

    The labels are stripped and replaced with "MOF_{element}", in accordance
    with the CrystalGenerator forcefield notation. Therefore, labels are also
    ignored.

    Args:
        data: the molecule, as a python object
    Returns:
        An instance of `pybel.Molecule`
    """
    if not PYBEL_LOADED:
        raise ImportError("Open Babel not installed.")

    table = pybel.ob.OBElementTable()

    if "building_blocks" in data:
        data["atoms"] = [a for bb in data["building_blocks"]
                         for a in bb["atoms"]]
    obmol = pybel.ob.OBMol()
    obmol.BeginModify()
    for atom in data["atoms"]:
        obatom = obmol.NewAtom()
        obatom.SetAtomicNum(table.GetAtomicNum(str(atom["element"])))
        obatom.SetVector(*atom["location"])

    uc = pybel.ob.OBUnitCell()
    uc.SetData(*(pybel.ob.vector3(*v) for v in data["unitcell"]))
    uc.SetSpaceGroup("P1")
    obmol.CloneData(uc)

    obmol.EndModify()

    mol = pybel.Molecule(obmol)

    # Add partial charges
    if "charge" in data["atoms"][0]:
        mol.OBMol.SetPartialChargesPerceived()
        for atom, pyatom in zip(data["atoms"], mol.atoms):
            pyatom.OBAtom.SetPartialCharge(atom["charge"])

    return mol


def pybel_to_raspa_cif(structure):
    """Converts instances of `pybel.Molecule` to a RASPA charged cif format."""
    if not PYBEL_LOADED:
        raise ImportError("Open Babel not installed.")

    table = pybel.ob.OBElementTable()
    uc = structure.unitcell
    a, b, c = uc.GetA(), uc.GetB(), uc.GetC()
    alpha, beta, gamma = uc.GetAlpha(), uc.GetBeta(), uc.GetGamma()

    cif = dedent("""
                 data_I
                 _chemical_name_common ''
                 _cell_length_a {a:.4f}
                 _cell_length_b {b:.4f}
                 _cell_length_c {c:.4f}
                 _cell_angle_alpha {alpha:.0f}
                 _cell_angle_beta {beta:.0f}
                 _cell_angle_gamma {gamma:.0f}
                 _space_group_name_H-M_alt 'P 1'
                 _space_group_name_Hall 'P 1'
                 loop_
                     _symmetry_equiv_pos_as_xyz
                     x,y,z
                 loop_
                     _atom_site_label
                     _atom_site_type_symbol
                     _atom_site_fract_x
                     _atom_site_fract_y
                     _atom_site_fract_z
                     _atom_site_charge
                 """.format(**locals())).strip()

    for atom in structure:
        element = table.GetSymbol(atom.atomicnum)
        label = "Mof_" + element
        charge = atom.partialcharge
        c = uc.WrapFractionalCoordinate(uc.CartesianToFractional(atom.vector))
        x, y, z = c.GetX(), c.GetY(), c.GetZ()

        cif += ("\n    {label:<7s} {element:<4s} {x:.5f} {y:9.5f} {z:9.5f} "
                "{charge:7.3f}").format(**locals())
    cif += "\n_end\n"

    return cif


def run_command_line():
    """Called by the `simulate` command, enables CLI interface."""
    parser = argparse.ArgumentParser(description="A general purpose classical "
                                     "simulation package.")
    parser.add_argument("-v", "--version", action="version",
                        version=__version__)
    parser.add_argument("input", type=str, help="a RASPA simulation script")
    parser.add_argument("-s", "--stream", action="store_true", help="Rather "
                        "than inputting and outputting files, this option "
                        "reads strings and outputs to stdout. Can be used to "
                        "simplify shell scripts. NOTE: Currently only "
                        "functional on linux.")
    parser.add_argument("-c", "--crystal", type=str, default="", help="If "
                        "streaming is enabled, this can be used to specify an "
                        "input structure without changing the RASPA script.")
    args = parser.parse_args()

    # Support both filepaths and data-streamed strings
    if args.stream:
        try:
            with open(args.input) as in_file:
                args.input = in_file.read()
        except IOError:
            pass
        try:
            with open(args.crystal) as in_file:
                args.crystal = in_file.read()
        except IOError:
            pass

    if args.stream:
        print(run_script(args.input, args.crystal, args.stream))
    else:
        run_script(args.input, args.crystal, args.stream)


def get_raspa_dir():
    """Called by the `raspa-dir` command, enables easy directory switching."""
    parser = argparse.ArgumentParser(description="Returns the data directory "
                                     "for the RASPA simulation program. Use "
                                     "'raspa-dir' to view the path, or 'cd "
                                     "`raspa-dir`' to change directories.")
    parser.add_argument("-v", "--version", action="version",
                        version=__version__)
    print(os.path.join(raspa_dir, "share/raspa"))


if __name__ == "__main__":
    run_command_line()
