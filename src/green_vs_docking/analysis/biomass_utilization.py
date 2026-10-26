"""Biomass utilisation efficiency (BUE) scoring with NameRxn atom mapping."""

from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from .utils import canonicalize_smiles


@dataclass(frozen=True)
class BiomassComputation:
    """Final BUE value plus numerator and denominator masses for debugging."""

    utilization: float
    biomass_mass: float
    block_input_mass: float


@dataclass(frozen=True)
class MoleculeNode:
    """Minimal molecule-node representation extracted from a route tracker."""

    node_id: str
    smiles: str
    canonical: str | None
    depth: int


@dataclass(frozen=True)
class ReactionNode:
    """Minimal reaction-node representation extracted from a route tracker."""

    node_id: str
    smiles: str
    depth: int


def build_block_mass_lookup(block_smiles: set[str]) -> Dict[str, float]:
    """Return exact molar masses keyed by canonical block SMILES."""
    masses: Dict[str, float] = {}
    for smi in block_smiles:
        canonical = canonicalize_smiles(smi)
        if canonical is None or canonical in masses:
            continue
        mol = Chem.MolFromSmiles(canonical)
        if mol is None:
            continue
        masses[canonical] = float(rdMolDescriptors.CalcExactMolWt(mol))
    return masses


def _collect_nodes(route_entry: dict) -> tuple[Dict[str, MoleculeNode], Dict[str, ReactionNode]]:
    """Extract molecule and reaction nodes from one Syntheseus route entry."""
    mol_nodes: Dict[str, MoleculeNode] = {}
    rxn_nodes: Dict[str, ReactionNode] = {}
    for key, payload in route_entry.items():
        if not key.startswith('node_') or not isinstance(payload, dict):
            continue
        depth = int(payload.get('depth', 0))
        if payload.get('is_mol'):
            smiles = payload.get('mol_smiles') or ''
            mol_nodes[key] = MoleculeNode(
                node_id=key,
                smiles=smiles,
                canonical=canonicalize_smiles(smiles),
                depth=depth,
            )
        elif payload.get('is_rxn'):
            rxn_nodes[key] = ReactionNode(node_id=key, smiles=payload.get('rxn_smiles') or '', depth=depth)
    return mol_nodes, rxn_nodes


def _canonical_atom_order(mol: Chem.Mol) -> List[int]:
    """Return atom indices in RDKit canonical rank order."""
    ranks = list(Chem.CanonicalRankAtoms(mol))
    return [idx for idx, rank in sorted(enumerate(ranks), key=lambda t: (t[1], t[0]))]


def _canonicalize_without_maps(smiles: str) -> str | None:
    """Canonicalize mapped NameRxn SMILES after removing atom-map numbers."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol)


def _flags_for_mol(mol: Chem.Mol, canonical_flags: Sequence[bool]) -> List[bool]:
    """Convert canonical-order atom flags back to the atom order of `mol`."""
    mol_copy = Chem.Mol(mol)
    for atom in mol_copy.GetAtoms():
        atom.SetAtomMapNum(0)
    order = _canonical_atom_order(mol_copy)
    if len(order) != len(canonical_flags):
        raise ValueError('Canonical flag length does not match atom count.')
    idx_to_flag = {order[i]: canonical_flags[i] for i in range(len(order))}
    return [idx_to_flag.get(idx, False) for idx in range(mol.GetNumAtoms())]


def _collapse_flags(mol: Chem.Mol, atom_flags: Sequence[bool]) -> tuple[bool, ...]:
    """Store atom flags in canonical atom order for stable multistep reuse."""
    if len(atom_flags) != mol.GetNumAtoms():
        raise ValueError('Atom flag length mismatch.')
    mol_copy = Chem.Mol(mol)
    for atom in mol_copy.GetAtoms():
        atom.SetAtomMapNum(0)
    order = _canonical_atom_order(mol_copy)
    return tuple(atom_flags[idx] for idx in order)


def _parse_reaction_smiles(rxn_smiles: str) -> tuple[List[str], List[str]]:
    """Split a reaction SMILES into reactant and product fragment lists."""
    if '>>' not in rxn_smiles:
        return [], []
    left, right = rxn_smiles.split('>>', 1)
    reactants = [frag.strip() for frag in left.split('.') if frag.strip()]
    products = [frag.strip() for frag in right.split('.') if frag.strip()]
    return reactants, products


def _map_products_to_reactions(
    mol_nodes: Dict[str, MoleculeNode],
    rxn_nodes: Dict[str, ReactionNode],
) -> tuple[Dict[str, str], Dict[str, str]]:
    """Map each reaction node to the molecule node it produces."""
    product_of_reaction: Dict[str, str] = {}
    produced_by: Dict[str, str] = {}
    for rxn in rxn_nodes.values():
        reactants, products = _parse_reaction_smiles(rxn.smiles)
        if not products:
            continue
        target_depth = rxn.depth - 1
        matched = None
        for prod in products:
            canonical = canonicalize_smiles(prod)
            if canonical is None:
                continue
            for node in mol_nodes.values():
                if node.depth == target_depth and node.canonical == canonical:
                    matched = node.node_id
                    break
            if matched:
                break
        if matched is None:
            continue
        product_of_reaction[rxn.node_id] = matched
        produced_by[matched] = rxn.node_id
    return product_of_reaction, produced_by


def _map_reactants_to_reactions(
    mol_nodes: Dict[str, MoleculeNode],
    rxn_nodes: Dict[str, ReactionNode],
) -> Dict[str, List[str]]:
    """Map each reaction node to molecule nodes used as its reactants."""
    buckets: Dict[tuple[int, str], List[str]] = {}
    for node in mol_nodes.values():
        if node.canonical is None:
            continue
        key = (node.depth, node.canonical)
        buckets.setdefault(key, []).append(node.node_id)
    for value in buckets.values():
        value.sort()

    reactants_map: Dict[str, List[str]] = {}
    for rxn in rxn_nodes.values():
        reactants, _ = _parse_reaction_smiles(rxn.smiles)
        assigned: List[str] = []
        depth_key = rxn.depth + 1
        for reactant in reactants:
            canonical = canonicalize_smiles(reactant)
            if canonical is None:
                continue
            key = (depth_key, canonical)
            nodes = buckets.get(key)
            if not nodes:
                continue
            assigned.append(nodes.pop(0))
        reactants_map[rxn.node_id] = assigned
    return reactants_map


def _run_namerxn(rxn_smiles: Iterable[str], binary_path: Path) -> List[str]:
    """Run the NameRxn binary on reaction SMILES and return mapped reactions."""
    rxn_list = list(rxn_smiles)
    if not rxn_list:
        return []
    with tempfile.TemporaryDirectory() as tmpdir:
        input_path = Path(tmpdir) / 'rxn.smi'
        output_path = Path(tmpdir) / 'rxn_out.smi'
        input_path.write_text("\n".join(rxn_list) + "\n")
        subprocess.run([str(binary_path), str(input_path), str(output_path)], check=True)
        lines = [line.strip() for line in output_path.read_text().splitlines() if line.strip()]
    if len(lines) != len(rxn_list):
        raise RuntimeError('NameRxn returned unexpected number of reactions.')
    return lines


def compute_biomass_utilization(
    route_entry: dict,
    block_smiles: set[str],
    block_mass_lookup: Dict[str, float],
    namerxn_binary: Path,
) -> BiomassComputation:
    """
    Compute BUE for one route.

    Starting molecules that match `block_smiles` define the biomass input mass
    and receive atom flags. NameRxn atom maps are then used reaction-by-reaction
    to propagate those flags to the final target product.
    """
    if not block_smiles:
        return BiomassComputation(utilization=0.0, biomass_mass=0.0, block_input_mass=0.0)

    mol_nodes, rxn_nodes = _collect_nodes(route_entry)
    if not mol_nodes or not rxn_nodes:
        return BiomassComputation(utilization=0.0, biomass_mass=0.0, block_input_mass=0.0)

    product_of_reaction, produced_by = _map_products_to_reactions(mol_nodes, rxn_nodes)
    reactants_per_reaction = _map_reactants_to_reactions(mol_nodes, rxn_nodes)

    molecule_flags: Dict[str, tuple[bool, ...]] = {}
    block_input_mass = 0.0
    for node in mol_nodes.values():
        # Starting materials are the molecule nodes that are not produced by any
        # reaction in this route. Only these can contribute new biomass input.
        if node.node_id in produced_by:
            continue
        mol = Chem.MolFromSmiles(node.smiles)
        if mol is None:
            continue
        canonical = node.canonical
        if canonical and canonical in block_mass_lookup:
            block_input_mass += block_mass_lookup.get(canonical, 0.0)
            atom_flags = [True] * mol.GetNumAtoms()
        else:
            atom_flags = [False] * mol.GetNumAtoms()
        # Flags are stored in canonical atom order so they remain comparable
        # when the same molecule later appears in a differently ordered SMILES.
        molecule_flags[node.node_id] = _collapse_flags(mol, atom_flags)

    if block_input_mass <= 0.0:
        return BiomassComputation(utilization=0.0, biomass_mass=0.0, block_input_mass=0.0)

    # Process from leaf reactions toward the target product so that any
    # intermediate has already received biomass flags before it is consumed.
    rxn_order = sorted(rxn_nodes.values(), key=lambda r: r.depth, reverse=True)
    mapped_lines = _run_namerxn([rxn.smiles for rxn in rxn_order], namerxn_binary)
    mapped_by_reaction = {rxn.node_id: line for rxn, line in zip(rxn_order, mapped_lines)}

    for rxn in rxn_order:
        mapped_line = mapped_by_reaction.get(rxn.node_id)
        if not mapped_line:
            continue
        if ' ' in mapped_line:
            # NameRxn appends class metadata after the mapped reaction; keep
            # only the reaction SMILES before splitting reactants/products.
            mapped_rxn, _ = mapped_line.rsplit(' ', 1)
        else:
            mapped_rxn = mapped_line
        mapped_reactants, mapped_products = _parse_reaction_smiles(mapped_rxn)
        reactant_nodes = reactants_per_reaction.get(rxn.node_id, [])
        unmatched = reactant_nodes.copy()
        reactant_records: List[tuple[str, str]] = []
        for mapped_reactant in mapped_reactants:
            # Mapped reactants still need to be matched back to tracker nodes;
            # removing atom maps makes this a normal canonical-SMILES match.
            canonical = _canonicalize_without_maps(mapped_reactant)
            if canonical is None:
                continue
            node_id = None
            for candidate in unmatched:
                if mol_nodes[candidate].canonical == canonical:
                    node_id = candidate
                    break
            if node_id is None:
                continue
            unmatched.remove(node_id)
            reactant_records.append((node_id, mapped_reactant))

        product_node_id = product_of_reaction.get(rxn.node_id)
        if not product_node_id:
            continue
        canonical_product = mol_nodes[product_node_id].canonical
        mapped_product_str = None
        for prod in mapped_products:
            if _canonicalize_without_maps(prod) == canonical_product:
                mapped_product_str = prod
                break
        if mapped_product_str is None:
            continue

        biomass_mapnums: set[int] = set()
        for node_id, mapped_smiles in reactant_records:
            mol = Chem.MolFromSmiles(mapped_smiles)
            if mol is None:
                continue
            canonical_flags = molecule_flags.get(node_id)
            if canonical_flags is None:
                continue
            # Convert saved canonical flags into the atom order used by the
            # mapped NameRxn reactant, then collect the map numbers of biomass atoms.
            atom_flags = _flags_for_mol(mol, canonical_flags)
            for atom, flag in zip(mol.GetAtoms(), atom_flags):
                amap = atom.GetAtomMapNum()
                if flag and amap:
                    biomass_mapnums.add(amap)

        product_mol = Chem.MolFromSmiles(mapped_product_str)
        if product_mol is None:
            continue
        product_atom_flags: List[bool] = []
        for atom in product_mol.GetAtoms():
            amap = atom.GetAtomMapNum()
            product_atom_flags.append(bool(amap and amap in biomass_mapnums))
        # Save product flags in canonical order; this product may be an
        # intermediate reactant in the next, shallower reaction.
        molecule_flags[product_node_id] = _collapse_flags(product_mol, product_atom_flags)

    root_node = min(mol_nodes.values(), key=lambda node: node.depth)
    root_mol = Chem.MolFromSmiles(root_node.smiles)
    if root_mol is None:
        return BiomassComputation(utilization=0.0, biomass_mass=0.0, block_input_mass=block_input_mass)
    canonical_flags = molecule_flags.get(root_node.node_id, tuple(False for _ in range(root_mol.GetNumAtoms())))
    atom_flags = _flags_for_mol(root_mol, canonical_flags)
    periodic_table = Chem.GetPeriodicTable()
    biomass_mass = 0.0
    for atom, flag in zip(root_mol.GetAtoms(), atom_flags):
        if flag:
            # Side products are absent from the route JSON, so the numerator is
            # simply the mass of biomass-tagged atoms retained in the target.
            biomass_mass += periodic_table.GetAtomicWeight(atom.GetAtomicNum())

    utilization = biomass_mass / block_input_mass if block_input_mass > 0 else 0.0
    return BiomassComputation(utilization=utilization, biomass_mass=biomass_mass, block_input_mass=block_input_mass)
