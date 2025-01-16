import json
import networkx as nx
from os import environ
from dotenv import load_dotenv
from rdkit import Chem
from rdkit import RDLogger


RDLogger.DisableLog('rdApp.*')


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


def check_authorization(provided_key: str) -> bool:
    allowed_key = environ.get('API_KEY')
    if not allowed_key:
        return False
    return provided_key == allowed_key


def read_vec2(data) -> tuple[float, float]:
    if isinstance(data, list):
        return (data[0] / 100, data[1] / 100)
    else:
        return (data['x'] / 100, data['y'] / 100)


def json_to_graph(data) -> nx.Graph:
    g = nx.Graph()

    for v, atom in enumerate(data['atoms']):
        g.add_node(v, **{
            'symbol': atom['symbol'],
            'position': read_vec2(atom['position']),
            'formal_charge': atom.get('formal_charge', 0),
            'non_bonded_ve': atom.get('non_bonded_ve', 0)
        })

    for bond in data['bonds']:
        g.add_edge(bond['from'], bond['to'], type=bond['type'])

    return g


def graph_to_mol(g: nx.Graph) -> Chem.rdchem.RWMol:
    mol = Chem.rdchem.RWMol()
    index = {}

    for v in g.nodes:
        a = Chem.rdchem.Atom(g.nodes[v]['symbol'])
        a.SetFormalCharge(g.nodes[v]['formal_charge'])
        index[v] = mol.AddAtom(a)

    for u, v in g.edges:
        type = g.edges[(u,v)]['type']
        bond = mol.AddBond(index[u], index[v], bond_type[type])
        if type in bond_dir:
            mol.GetBondWithIdx(bond - 1).SetBondDir(bond_dir[type])

    # Build 2D conformer.
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(False)
    for v in g.nodes:
        # Note that we add a Z-coordinate with value 0.
        conf.SetAtomPosition(index[v], g.nodes[v]['position'] + (0,))

    # Assign conformer and compute stereochemistry.
    mol = mol.GetMol()
    mol.AddConformer(conf)

    Chem.SanitizeMol(mol)
    Chem.DetectBondStereochemistry(mol)
    Chem.AssignChiralTypesFromBondDirs(mol)
    Chem.AssignStereochemistry(mol)

    return mol


def json_to_mol(data) -> Chem.rdchem.RWMol:
   return graph_to_mol(json_to_graph(data))


def _node_match(v1, v2):
    props = [
        ['symbol', ''],
        ['formal_charge', 0],
        ['non_bonded_ve', 0],
        ['mark', -1]
    ]
    return all(map(lambda p: v1.get(p[0], p[1]) == v2.get(p[0], p[1]), props))


def _edge_match(e1, e2):
    return e1['type'] == e2['type']


def _graph_to_smiles(g: nx.Graph) -> str:
    mol = graph_to_mol(g)
    return Chem.rdmolfiles.MolToSmiles(mol, isomericSmiles=True)


def compare_components(g1: nx.Graph, g2: nx.Graph, match_stereo: bool) -> bool:
    if nx.is_isomorphic(g1, g2, _node_match, _edge_match):
        if match_stereo:
            try:
                return _graph_to_smiles(g1) == _graph_to_smiles(g2)
            except ValueError:
                return False
        else:
            return True


def compare(diagram1, diagram2, match_stereo = True) -> bool:
    g1 = json_to_graph(diagram1)
    g2 = json_to_graph(diagram2)
    cs1 = list(map(lambda c: g1.subgraph(c), nx.connected_components(g1)))
    cs2 = list(map(lambda c: g2.subgraph(c), nx.connected_components(g2)))

    # Check if every component in g1 matches a unique component in g2.
    for c1 in cs1:
        match = (i for i, c2 in enumerate(cs2) if
                 compare_components(c1, c2, match_stereo))
        index = next(match, None)
        if index is None:
            return False
        else:
            # Every component can only be matched once.
            cs2.pop(index)

    return len(cs2) == 0


def validate(diagram) -> bool:
    try:
        json_to_mol(diagram)
        return True
    except ValueError:
        return False


def handle_validate(body):
    return {
        'valid': validate(body['diagram'])
    }


def handle_compare(body):
    return {
        'equal': compare(
            body['reference_diagram'],
            body['student_diagram'],
            body['match_stereo'])
    }


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
                response_body = handle_compare(request_body)
            case 'validate':
                response_body = handle_validate(request_body)
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
