from dotenv import load_dotenv
from os import environ
import json
from rdkit import Chem

bond_type = {
    'single': Chem.rdchem.BondType.SINGLE,
    'double': Chem.rdchem.BondType.DOUBLE,
    'triple': Chem.rdchem.BondType.TRIPLE,
    'single_up': Chem.rdchem.BondType.SINGLE,
    'single_down': Chem.rdchem.BondType.SINGLE
}


def check_authorization(provided_key: str) -> bool:
    allowed_key = environ.get('API_KEY')
    if not allowed_key:
        return False
    return provided_key == allowed_key


def read_mol(data) -> Chem.rdchem.RWMol:
    mol = Chem.rdchem.RWMol()

    for atom in data['atoms']:
        a = Chem.rdchem.Atom(atom['symbol'])
        if 'formal_charge' in atom:
            a.SetFormalCharge(atom['formal_charge'])
        mol.AddAtom(a)

    for bond in data['bonds']:
        mol.AddBond(
            bond['from'],
            bond['to'],
            bond_type[bond['type']])

    return mol


def compare(data1, data2) -> bool:
    mol1 = read_mol(data1)
    mol2 = read_mol(data2)
    smiles1 = Chem.rdmolfiles.MolToSmiles(mol1)
    smiles2 = Chem.rdmolfiles.MolToSmiles(mol2)
    return smiles1 == smiles2


def validate(diagram) -> bool:
    mol = read_mol(diagram)
    try:
        Chem.SanitizeMol(mol)
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
            body['student_diagram'])
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
