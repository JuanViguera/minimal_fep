import argparse
from minimal_fep.setup.setup import load_ligands, create_ligand_network, create_chemical_system, get_simulation_settings_and_protocol
from minimal_fep.compute.compute import compute_rbfe

def main():
    parser = argparse.ArgumentParser(description="Script for computing relative binding free energies using openfe.")

    parser.add_argument('--ligands_path', type=str, required=True, help="Path to an sdf file containing ligands.")
    parser.add_argument('--protein_path', type=str, required=True, help="Path to the pdb file containing the protein structure.")

    parser.add_argument('--atom_mapper', type=str, default='lomap', choices=['lomap', 'kartograf'], help="Atom mapper to use for mapping ligands.")
    parser.add_argument('--network_type', type=str, default='mst', choices=['mst', 'lomap', 'radial'], help="Type of ligand network to create.")
    parser.add_argument('--ligand_name', type=str, default="lig_ejm_47", help="Name of the ligand to use for the simulation.")
    #solvent components
    parser.add_argument('--positive_ion', type=str, default='Na', help="Positive ion to use in the solvent component.")
    parser.add_argument('--negative_ion', type=str, default='Cl', help="Negative ion to use in the solvent component.")
    parser.add_argument('--no_neutralize', action='store_true', help="Whether to neutralize the system with ions.")
    parser.add_argument('--ion_concentration', type=float, default=0.15, help="Ion concentration for neutralization (mol/L).")
    #simulation settings
    parser.add_argument('--early_termination_target_error', type=float, default=0.0, help="Early termination target error in kcal/mol.")
    parser.add_argument('--equilibration_length', type=float, default=10.0, help="Equilibration length in picoseconds.")
    parser.add_argument('--minimization_steps', type=int, default=5000, help="Number of minimization steps.")
    parser.add_argument('--n_replicas', type=int, default=11, help="Number of replicas for the simulation.")
    parser.add_argument('--production_length', type=float, default=50.0, help="Production length in picoseconds.")
    parser.add_argument('--real_time_analysis_interval', type=float, default=250.0, help="Real-time analysis interval in picoseconds.")
    parser.add_argument('--real_time_analysis_minimum_time', type=float, default=500.0, help="Minimum time for real-time analysis in picoseconds.")
    parser.add_argument('--sampler_method', type=str, default='repex', choices=['repex'], help="Sampler method to use.")
    parser.add_argument('--sams_flatness_criteria', type=str, default='logZ-flatness', help="SAMS flatness criteria.")
    parser.add_argument('--sams_gamma0', type=float, default=1.0, help="Initial gamma0 value for SAMS.")
    parser.add_argument('--time_per_iteration', type=float, default=1.0, help="Time per iteration in picoseconds.")

    args = parser.parse_args()
    solvent_params = {
        'positive_ion': args.positive_ion,
        'negative_ion': args.negative_ion,
        'neutralize': not args.no_neutralize,
        'ion_concentration': args.ion_concentration
    }
    simulation_settings = {
        'early_termination_target_error': args.early_termination_target_error,
        'equilibration_length': args.equilibration_length,
        'minimization_steps': args.minimization_steps,
        'n_replicas': args.n_replicas,
        'production_length': args.production_length,
        'real_time_analysis_interval': args.real_time_analysis_interval,
        'real_time_analysis_minimum_time': args.real_time_analysis_minimum_time,
        'sampler_method': args.sampler_method,
        'sams_flatness_criteria': args.sams_flatness_criteria,
        'sams_gamma0': args.sams_gamma0,
        'time_per_iteration': args.time_per_iteration
    }

    ligand_mols = load_ligands(args.ligands_path)
    ligand_network = create_ligand_network(ligand_mols, args.atom_mapper, args.network_type)
    stateA_complex, stateB_complex, edge_mapping = create_chemical_system(
        args.protein_path, 
        ligand_network, 
        ligand_name=args.ligand_name,
        solvent_params=solvent_params)
    protocol = get_simulation_settings_and_protocol(simulation_settings)

    compute_rbfe(stateA_complex, stateB_complex, edge_mapping, protocol)

if __name__ == "__main__":
    main()