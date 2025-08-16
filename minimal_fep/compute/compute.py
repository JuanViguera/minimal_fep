import openfe
from gufe.protocols.protocoldag import execute_DAG
from gufe import AtomMapper, Protocol
import pathlib

def compute_rbfe(
    stateA_complex: openfe.ChemicalSystem,
    stateB_complex: openfe.ChemicalSystem,
    mapping: AtomMapper,
    rbfe_protocol: Protocol
) -> None:
    """
    Computes the relative binding free energy (RBFE) between two states in a complex system.
    This function sets up and executes an RBFE transformation using the provided states, atom mapping,
    and protocol. It performs a dry-run of the protocol unit, creates necessary result directories,
    executes the transformation DAG, and gathers the results to estimate the free energy difference.
    Args:
        stateA_complex: The initial state of the complex (e.g., ligand-protein system).
        stateB_complex: The final state of the complex after transformation.
        mapping: An atom mapping object describing correspondence between stateA and stateB.
        rbfe_protocol: The protocol object defining the RBFE calculation procedure.
    Returns:
        None. Prints the estimated free energy difference (dG) and its uncertainty to stdout.
    Raises:
        Any exceptions raised during transformation setup, execution, or result gathering.
    """

    transformation_complex = openfe.Transformation(
            stateA=stateA_complex,
            stateB=stateB_complex,
            mapping=mapping,
            protocol=rbfe_protocol, 
            name=f"{stateA_complex.name}_{stateB_complex.name}_complex"
        )
    complex_dag = transformation_complex.create()
    # complex dry-run
    complex_unit = list(complex_dag.protocol_units)[0]

    complex_unit.run(dry=True, verbose=True)    

    # Finally we can run the simulations
    complex_path = pathlib.Path('./results/complex')
    complex_path.mkdir(parents=True, exist_ok=True)

    # First the complex transformation
    complex_dag_results = execute_DAG(complex_dag, scratch_basedir=complex_path, shared_basedir=complex_path)

    # Get the complex and solvent results
    complex_results = rbfe_protocol.gather([complex_dag_results])

    print(f"Complex dG: {complex_results.get_estimate()}, err {complex_results.get_uncertainty()}")