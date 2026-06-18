import json
import networkx as nx
from os import environ
from dotenv import load_dotenv
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem, rdmolops, Descriptors


RDLogger.DisableLog('rdApp.*')


# === Reading MolView diagram JSON ===

bond_type = {
    'single': Chem.rdchem.BondType.SINGLE,
    'double': Chem.rdchem.BondType.DOUBLE,
    'triple': Chem.rdchem.BondType.TRIPLE,
    'single_up': Chem.rdchem.BondType.SINGLE,
    'single_down': Chem.rdchem.BondType.SINGLE
}


bond_dir = {
    'single_up': Chem.BondDir.BEGINWEDGE,
    'single_down': Chem.BondDir.BEGINDASH
}


def vec2(data) -> tuple[float, float]:
    if isinstance(data, list):
        return (data[0] / 100, data[1] / 100)
    else:
        return (data['x'] / 100, data['y'] / 100)


def lone_pairs(atom: dict):
    if 'non_bonded_ve' in atom:
        return atom['non_bonded_ve'] // 2
    else:
        return atom.get('lone_pairs', 0)


def unpaired_electrons(atom: dict):
    if 'non_bonded_ve' in atom:
        return atom['non_bonded_ve'] % 2
    else:
        return atom.get('unpaired_electrons', 0)


def bond_edge(bond):
    if 'atoms' in bond:
        return bond['atoms']
    else:
        return [bond['from'], bond['to']]


def attach_arrow_endpoint(g: nx.Graph, endpoint, anchor, bonds):
    if 'atom' in anchor:
        g.add_edge(anchor['atom'], endpoint, type='arrow')
    elif 'bond' in anchor:
        i, j = bond_edge(bonds[anchor['bond']])
        g.add_edge(i, endpoint, type='arrow')
        g.add_edge(j, endpoint, type='arrow')
        g.edges[(i,j)]['attached_to_arrow'] = True


def json_to_graph(diagram, add_arrows=False) -> nx.Graph:
    g = nx.Graph()

    for i, atom in enumerate(diagram['atoms']):
        g.add_node(atom.get('id') or i, **{
            'type': 'atom',
            'label': atom.get('label') or atom['symbol'],
            'position': vec2(atom['position']),
            'formal_charge': atom.get('formal_charge', 0),
            'lone_pairs': lone_pairs(atom),
            'unpaired_electrons': unpaired_electrons(atom)
        })

    for bond in diagram['bonds']:
        i, j = bond_edge(bond)
        g.add_edge(i, j, type='bond', bond_type=bond['type'])

    if 'arrows' in diagram and add_arrows:
        for i, arrow in enumerate(diagram['arrows']):
            start = f'arrow{i}start'
            end = f'arrow{i}end'
            attr = { 'type': 'arrow', 'arrow_type': arrow['type'] }
            g.add_node(start, **attr, point='start')
            g.add_node(end, **attr, point='end')
            g.add_edge(start, end, type='arrow')
            attach_arrow_endpoint(g, start, arrow['start'], diagram['bonds'])
            attach_arrow_endpoint(g, end, arrow['end'], diagram['bonds'])

    return g


def graph_to_mol(g: nx.Graph) -> Chem.rdchem.RWMol:
    mol = Chem.rdchem.RWMol()
    atoms = {}

    # 2D conformer.
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(False)

    for i, attr in g.nodes(data=True):
        if attr['type'] == 'atom':
            a = Chem.rdchem.Atom(attr['label'])
            a.SetFormalCharge(attr['formal_charge'])
            atom_index = mol.AddAtom(a)
            atoms[i] = atom_index

            # Update conformer. Note that we add a Z-coordinate with value 0.
            conf.SetAtomPosition(atom_index, attr['position'] + (0,))

    for i, j, attr in g.edges(data=True):
        if attr['type'] == 'bond':
            type = attr['bond_type']
            bond_count = mol.AddBond(atoms[i], atoms[j], bond_type[type])
            bond = mol.GetBondWithIdx(bond_count - 1)
            if type in bond_dir:
                bond.SetBondDir(bond_dir[type])

    # Assign conformer and compute stereochemistry.
    mol = mol.GetMol()
    mol.AddConformer(conf)

    # SanitizeMol applies some standardizations and sets some fields. This is
    # critical for subsequent stereochemistry operations. However, we don't want
    # to fail if non-standard valences are used (Example: ClF3).
    try:
        Chem.SanitizeMol(mol)
    except:
        pass

    Chem.DetectBondStereochemistry(mol)
    Chem.AssignChiralTypesFromBondDirs(mol)
    Chem.AssignStereochemistry(mol)

    return mol


def json_to_mol(data) -> Chem.rdchem.RWMol:
   return graph_to_mol(json_to_graph(data))


# === Utilities based on NetworkX ===


def match_props(x, y, props):
    return all(map(lambda p: x.get(p[0], p[1]) == y.get(p[0], p[1]), props))


def match_nodes(v1, v2):
    # We have atom nodes and arrow nodes.
    if v1['type'] == 'atom':
        return match_props(v1, v2, [
            ['type', ''],
            ['label', ''],
            ['formal_charge', 0],
            ['lone_pairs', 0],
            ['unpaired_electrons', 0]
        ])
    elif v1['type'] == 'arrow':
        return match_props(v1, v2, [
            ['type', ''],
            ['arrow_type', ''],
            ['point', '']
        ])
    else:
        return False


def match_edges(e1, e2):
    if e1['type'] != e2['type']:
        return False
    # We match bond types on bonds attached to arrows.
    if e1['type'] == 'bond' and 'attached_to_arrow' in e1:
        return e1['bond_type'] == e2['bond_type']
    else:
        return True


def graph_basic_hydrogens(g: nx.Graph) -> list:
    # Select hdrogen atoms _where all other attributes are empty_.
    return [v for v, attr in g.nodes(data=True) if
            attr['type'] == 'atom' and
            attr['label'] == 'H' and
            attr['formal_charge'] == 0 and
            attr['unpaired_electrons'] == 0]


def compare_component(g1: nx.Graph, g2: nx.Graph, match_stereo: bool) -> bool:
    g1_no_h = g1.copy()
    g2_no_h = g2.copy()
    g1_no_h.remove_nodes_from(graph_basic_hydrogens(g1_no_h))
    g2_no_h.remove_nodes_from(graph_basic_hydrogens(g2_no_h))

    # Don't compare bond types to avoid filtering resonance structures.
    if nx.is_isomorphic(g1_no_h, g2_no_h, match_nodes, match_edges):
        # This can fail because of non-standard valences, like in ClF3, which is
        # not supported by RDKit. Therefore, in the case of an error, we default
        # to accepting the bare isomorphism.
        try:
            smi1 = get_smiles(graph_to_mol(g1), match_stereo)
            smi2 = get_smiles(graph_to_mol(g2), match_stereo)
            return smi1 == smi2
        except ValueError:
            return True

    else:
        return False


# === Utilities based on RDKit ===


def get_smiles(mol: Chem.Mol, isomeric=True) -> str:
    mol = Chem.RemoveHs(mol)
    return Chem.MolToSmiles(mol, isomericSmiles=isomeric)


def get_fragments(mol):
    return rdmolops.GetMolFrags(mol, asMols=True)


# === Main operations ===


def validate(diagram) -> bool:
    try:
        json_to_mol(diagram)
        return True
    except Exception:
        return False


def compare(diagram1, diagram2, match_stereo=True) -> bool:
    g1 = json_to_graph(diagram1, True)
    g2 = json_to_graph(diagram2, True)
    cs1 = list(map(lambda c: g1.subgraph(c), nx.connected_components(g1)))
    cs2 = list(map(lambda c: g2.subgraph(c), nx.connected_components(g2)))

    # Check if every component in g1 matches a unique component in g2.
    for c1 in cs1:
        matched = False
        for i, c2 in enumerate(cs2):
            if compare_component(c1, c2, match_stereo):
                matched = True
                cs2.pop(i)
                break
        if not matched:
            return False

    return len(cs2) == 0


def evaluate(diagram, evaluator, params):
    mol = json_to_mol(diagram)
    match evaluator:
        case 'NumMols':
            return len(set(map(get_smiles, get_fragments(mol))))
        case _:
            return {'err': 'Unknown evaluator'}


# === AWS Lambda handler ===


def check_authorization(provided_key: str) -> bool:
    allowed_key = environ.get('API_KEY')
    if not allowed_key:
        return False
    return provided_key == allowed_key


def handler(event, context):
    try:
        load_dotenv()
        headers = event.get('headers', {})
        auth_header: str = headers.get('Authorization') or headers.get('authorization')
        if not auth_header:
            return {
                'body': json.dumps({'err': 'Authorization Required'}),
                'headers': {
                    'Content-Type': 'application/json'
                },
                'statusCode': 401,
            }
        auth_header = auth_header.replace('Bearer ', '')

        if not check_authorization(auth_header):
            return {
                'body': json.dumps({'err': 'Invalid authorization token'}),
                'headers': {
                    'Content-Type': 'application/json'
                },
                'statusCode': 401,
            }

        request_path_raw: str = event.get('requestContext', {}).get('path', '')
        request_path = request_path_raw.replace('/api/v1/', '')
        request_body = json.loads(event['body'])

        match request_path:
            case 'compare':
                response_body = {
                    'equal': compare(
                        request_body['reference_diagram'],
                        request_body['student_diagram'],
                        request_body['match_stereo'])
                }
            case 'validate':
                response_body = {
                    'valid': validate(
                        request_body['diagram'])
                }
            case 'evaluate':
                response_body = {
                    'result': evaluate(
                        request_body['diagram'],
                        request_body['evaluator'],
                        request_body['params'])
                }
            case _:
                response_body = {'err': 'Invalid request'}

        return {
            'body': json.dumps(response_body),
            'headers': {
                'Content-Type': 'application/json'
            },
            'statusCode': 200,
        }
    except Exception as e:
        print(e)
        return {
            'body': json.dumps({'err': 'Unable to parse diagrams or an internal error occurred.'}),
            'headers': {
                'Content-Type': 'application/json'
            },
            'statusCode': 400,
        }
