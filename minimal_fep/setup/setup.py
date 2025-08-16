from typing import Tuple, Any, Dict, List
from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
from rdkit import Chem
import openfe
from openfe import SmallMoleculeComponent, SolventComponent, ProteinComponent, ChemicalSystem
from openff.units import unit
from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
from gufe import LigandNetwork
from gufe.mapping.ligandatommapping import LigandAtomMapping


def load_ligands(ligands_path: str) -> list[SmallMoleculeComponent]:
    """
    Load ligand molecules from an SDF file using RDKit.

    Parameters
    ----------
    ligands_path : str
        Path to the SDF file containing ligand structures.

    Returns
    -------
    list of SmallMoleculeComponent
        A list of SmallMoleculeComponent objects created from the molecules in the SDF file.

    Notes
    -----
    This function uses RDKit's SDMolSupplier to read molecules from the SDF file.
    Hydrogens are not removed during loading.
    """
    # Load ligands using RDKit
    ligands_sdf = Chem.SDMolSupplier(ligands_path, removeHs=False)

    # Now pass these to form a list of Molecules
    ligand_mols = [SmallMoleculeComponent(sdf) for sdf in ligands_sdf]

    return ligand_mols

def create_ligand_network(
    ligand_mols: List[SmallMoleculeComponent],
    atom_mapper: str,
    network_type: str
) -> LigandNetwork:
    """
    Create a ligand network for free energy calculations using specified atom mapping and network type.
    Parameters
    ----------
    ligand_mols : list
        A list of ligand molecule objects to be included in the network.
    atom_mapper : str
        The atom mapping algorithm to use. Supported values are:
        - 'lomap': Use LomapAtomMapper from openfe.
        - 'kartograf': Use KartografAtomMapper from kartograf.
    network_type : str
        The type of ligand network to generate. Supported values are:
        - 'mst': Minimal spanning tree network.
        - 'lomap': Lomap network.
        - 'radial': Radial network with the first ligand as the central node.
    Returns
    -------
    network : object
        The generated ligand network object suitable for free energy calculations.
    Raises
    ------
    ValueError
        If an unknown atom_mapper or network_type is provided.
    """

    from openfe.setup import LomapAtomMapper
    if atom_mapper == 'lomap':
        from openfe.setup import LomapAtomMapper
        mapper = LomapAtomMapper()
    elif atom_mapper == 'kartograf':
        from kartograf import KartografAtomMapper
        mapper = KartografAtomMapper()
    else:
        raise ValueError(f"Unknown atom mapper: {atom_mapper}")
    
    if network_type == 'mst':
        from openfe.setup.ligand_network_planning import generate_minimal_spanning_network
        network = generate_minimal_spanning_network(
            ligands=ligand_mols, 
            scorer=openfe.lomap_scorers.default_lomap_score, 
            mappers=[mapper,]
        )

    elif network_type == 'lomap':
        from openfe.setup.ligand_network_planning import generate_lomap_network
        network = generate_lomap_network(
            molecules=ligand_mols,
            scorer=openfe.lomap_scorers.default_lomap_score,
            mappers=[mapper,]
        )
    elif network_type == 'radial':
        from openfe.setup.ligand_network_planning import generate_radial_network 
        network = generate_radial_network(
            ligands=ligand_mols[1:],
            central_ligand=ligand_mols[0],
            mappers=[mapper,]
        )
    else:
        raise ValueError(f"Unknown network type: {network_type}")
    return network

def create_chemical_system(
    protein_path: str,
    network: LigandNetwork,
    ligand_name: str,
    solvent_params: Dict[str, Any]
) -> Tuple[ChemicalSystem, ChemicalSystem, LigandAtomMapping]:
    """
    Creates two ChemicalSystem objects representing the complex of a protein with two ligand states (A and B),
    along with a solvent environment, for use in free energy perturbation simulations.
    Args:
        protein_path (str): Path to the protein PDB file.
        network (object): A network object containing edges with ligand components.
        ligand_name (str): The name of the ligand to select from the network.
        solvent_params (dict): Dictionary containing solvent parameters:
            - 'positive_ion' (str): Name of the positive ion.
            - 'negative_ion' (str): Name of the negative ion.
            - 'neutralize' (bool): Whether to neutralize the system.
            - 'ion_concentration' (float): Ion concentration in molar units.
    Returns:
        tuple: A tuple containing:
            - stateA_complex (ChemicalSystem): Chemical system with ligand state A.
            - stateB_complex (ChemicalSystem): Chemical system with ligand state B.
            - chosen_edge (object): The network edge corresponding to the selected ligand.
    """
    
    # Create the chemical system for simulation
    protein = ProteinComponent.from_pdb_file(protein_path)

    solvent = SolventComponent(positive_ion=solvent_params['positive_ion'], 
                               negative_ion=solvent_params['negative_ion'],
                               neutralize=solvent_params['neutralize'], 
                               ion_concentration=solvent_params['ion_concentration']*unit.molar)

    chosen_edge = [edge for edge in network.edges if edge.componentB.name == ligand_name][0]

    stateA_complex = ChemicalSystem({'ligand': chosen_edge.componentA,
                                  'solvent': solvent,
                                  'protein': protein,},
                               name=chosen_edge.componentA.name)
    stateB_complex = ChemicalSystem({'ligand': chosen_edge.componentB,
                                 'solvent': solvent,
                                 'protein': protein,},
                               name=chosen_edge.componentB.name)

    return stateA_complex, stateB_complex, chosen_edge


def get_simulation_settings_and_protocol(
    simulation_settings: Dict[str, float | int | str]
) -> RelativeHybridTopologyProtocol:
    """
    Configures and returns simulation settings and protocol for relative binding free energy (RBFE) calculations.
    Parameters
    ----------
    simulation_settings : dict
        Dictionary containing simulation parameters. Expected keys include:
            - 'equilibration_length' (float): Equilibration time in picoseconds.
            - 'production_length' (float): Production time in picoseconds.
            - 'minimization_steps' (int): Number of minimization steps.
            - 'n_replicas' (int): Number of simulation replicas.
            - 'real_time_analysis_interval' (float): Interval for real-time analysis in picoseconds.
            - 'real_time_analysis_minimum_time' (float): Minimum time for real-time analysis in picoseconds.
            - 'sampler_method' (str): Sampling method to use.
            - 'sams_flatness_criteria' (float): Flatness criteria for SAMS.
            - 'sams_gamma0' (float): Initial gamma value for SAMS.
            - 'time_per_iteration' (float): Simulation time per iteration in picoseconds.
    Returns
    -------
    rbfe_settings : RelativeHybridTopologyProtocol.Settings
        Configured RBFE simulation settings object.
    rbfe_protocol : RelativeHybridTopologyProtocol
        RBFE protocol object initialized with the configured settings.
    """

    # Get the simulation settings and protocol
    # Create the default settings
    rbfe_settings = RelativeHybridTopologyProtocol.default_settings()

    rbfe_settings.simulation_settings.equilibration_length = simulation_settings['equilibration_length'] * unit.picosecond
    rbfe_settings.simulation_settings.production_length = simulation_settings['production_length'] * unit.picosecond
    rbfe_settings.simulation_settings.minimization_steps = simulation_settings['minimization_steps']
    rbfe_settings.simulation_settings.n_replicas = simulation_settings['n_replicas']
    rbfe_settings.simulation_settings.real_time_analysis_interval = simulation_settings['real_time_analysis_interval'] * unit.picosecond
    rbfe_settings.simulation_settings.real_time_analysis_minimum_time = simulation_settings['real_time_analysis_minimum_time'] * unit.picosecond
    rbfe_settings.simulation_settings.sampler_method = simulation_settings['sampler_method']
    rbfe_settings.simulation_settings.sams_flatness_criteria = simulation_settings['sams_flatness_criteria']
    rbfe_settings.simulation_settings.sams_gamma0 = simulation_settings['sams_gamma0']
    rbfe_settings.simulation_settings.time_per_iteration = simulation_settings['time_per_iteration'] * unit.picosecond

    # Create RBFE Protocol class
    rbfe_protocol = RelativeHybridTopologyProtocol(
        settings=rbfe_settings
    )
    return rbfe_protocol