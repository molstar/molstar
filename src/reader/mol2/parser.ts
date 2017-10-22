//               NOTES                
//When want to created undefined string column, must use 
// undefStr = UndefinedColumn(molecule.num_atoms, ColumnType.str)
// but not 
// const undefPooledStr = UndefinedColumn(molecule.num_atoms, ColumnType.pooledStr);
// because latter actuall return a column of zeros
import Tokenizer from '../common/text/tokenizer'
import FixedColumn from '../common/text/column/fixed'
import { ColumnType, UndefinedColumn } from '../common/column'
import * as Schema from './schema'
import Result from '../result'
import Computation from '../../utils/computation' 

interface State {
    tokenizer: Tokenizer,
    molecule: Schema.Molecule,
    chunker: Computation.Chunker
}



function createEmptyMolecule(): Schema.Molecule {
    return {
        mol_name: '',
        num_atoms: 0,
        num_bonds: 0,
        num_subst: 0,
        num_feat: 0,
        num_sets: 0,
        mol_type: '',
        charge_type: '',
        status_bits:'',
        mol_comment: ''
    };
}




function State(tokenizer: Tokenizer, ctx: Computation.Context): State { 
    return {
        tokenizer,
        molecule: createEmptyMolecule(),
        chunker: Computation.chunker(ctx, 100000)
    };
}





function handleMolecule(state: State) {
    const { tokenizer, molecule } = state;
    Tokenizer.markLine(tokenizer); // skip the line '@<TRIPOS>MOLECULE'
    Tokenizer.markLine(tokenizer);
    let name = Tokenizer.getTokenString(tokenizer);
    molecule.mol_name = name;

    Tokenizer.markLine(tokenizer);
    const values = Tokenizer.getTokenString(tokenizer).trim().split(/\s+/g);
    molecule.num_atoms = parseInt(values[0]) ? parseInt(values[1]) : 0;
    molecule.num_bonds = parseInt(values[1]) ? parseInt(values[1]) : 0;
    molecule.num_subst = parseInt(values[2]) ? parseInt(values[1]) : 0;
    molecule.num_feat = parseInt(values[3]) ? parseInt(values[1]) : 0;
    molecule.num_sets = parseInt(values[4]) ? parseInt(values[1]) : 0;

    Tokenizer.markLine(tokenizer);
    molecule.mol_type = Tokenizer.getTokenString(tokenizer);

    Tokenizer.markLine(tokenizer);
    molecule.charge_type = Tokenizer.getTokenString(tokenizer);

    Tokenizer.markLine(tokenizer);
    if(Tokenizer.getTokenString(tokenizer) == ''){return}
    else{molecule.status_bits = Tokenizer.getTokenString(tokenizer)}


    Tokenizer.markLine(tokenizer);
    if(Tokenizer.getTokenString(tokenizer) == ''){return}
    else{molecule.mol_comment = Tokenizer.getTokenString(tokenizer)}
}





async function handleAtoms(state: State): Promise<Schema.Atoms> {
    const { tokenizer, molecule } = state;
    let hasSubst_id = false;
    let hasSubst_name = false;
    let hasCharge = false;
    let hasStatus_bit = false;

    // skip empty lines and '@<TRIPOS>ATOM'
    while(Tokenizer.getTokenString(tokenizer) != '@<TRIPOS>ATOM'){
        Tokenizer.markLine(tokenizer);
    }
    const lines =  await Tokenizer.readLinesAsync(tokenizer, molecule.num_atoms, state.chunker);
    const firstLine = tokenizer.data.substring(lines.indices[0], lines.indices[1]);
    const firstLineArray = firstLine.trim().split(/\s+/g)
    const length = firstLineArray.length;
    if(length == 9){
        hasSubst_id = true;
        hasSubst_name = true;
        hasCharge = true;
    }else if(length == 10){
        hasSubst_id = true;
        hasSubst_name = true;
        hasCharge = true;
        hasStatus_bit = true;
    }    

    const col = FixedColumn(lines);
    const undefInt = UndefinedColumn(molecule.num_atoms, ColumnType.int);
    const undefFloat = UndefinedColumn(molecule.num_atoms, ColumnType.float);
    //const undefPooledStr = UndefinedColumn(molecule.num_atoms, ColumnType.pooledStr);
    // created below column to pass unit tests
    const undefStr = UndefinedColumn(molecule.num_atoms, ColumnType.str);

    const ret = {
        count: molecule.num_atoms,
        atom_id: col(0, 7, ColumnType.int),
        atom_name: col(7, 9, ColumnType.pooledStr), 
        x: col(16, 10, ColumnType.float),
        y: col(26, 10, ColumnType.float),
        z: col(36, 10, ColumnType.float),
        atom_type: col(46, 6, ColumnType.pooledStr),
        // optional properties
        subst_id: hasSubst_id ? col(52, 6, ColumnType.int) : undefInt, 
        subst_name: hasSubst_name ? col(58, 8, ColumnType.pooledStr) : undefStr,
        charge: hasCharge ? col(66, 10, ColumnType.float) : undefFloat, 
        // undefPooledStr cannot pass unit tests because it create a column of zeros but not empty strings
        //status_bit: hasStatus_bit ? col(76, 100, ColumnType.pooledStr) : undefPooledStr, 

        // use undefStr instead to pass unit tests
        status_bit: hasStatus_bit ? col(76, 100, ColumnType.pooledStr) : undefStr, 

    };

    return ret;
}




async function handleBonds(state: State): Promise<Schema.Bonds> {
    const { tokenizer, molecule } = state;
    let hasStatus_bit = false;

    while(Tokenizer.getTokenString(tokenizer) != '@<TRIPOS>BOND'){
        Tokenizer.markLine(tokenizer);
    }
    const lines =  await Tokenizer.readLinesAsync(tokenizer, molecule.num_bonds, state.chunker);
    const firstLine = tokenizer.data.substring(lines.indices[0], lines.indices[1]);
    const length = firstLine.split(' ').length;
    if(length == 4){
        hasStatus_bit = true;
    }

    const col = FixedColumn(lines);
    //const undefPooledStr = UndefinedColumn(molecule.num_bonds, ColumnType.pooledStr);
    // created below column to pass unit tests
    const undefStr = UndefinedColumn(molecule.num_atoms, ColumnType.str);

    const ret = {
        count: molecule.num_bonds,
        bond_id: col(0, 6, ColumnType.int),
        origin_atom_id: col(6, 6, ColumnType.int), 
        target_atom_id: col(12, 6, ColumnType.int),
        bond_type: col(18, 5, ColumnType.pooledStr), 
        // optional properties
        // undefPooledStr cannot pass unit tests because it create a column of zeros but not empty strings
        //status_bits: hasStatus_bit ? col(23, 50, ColumnType.pooledStr) : undefPooledStr, 
        status_bits: hasStatus_bit ? col(23, 50, ColumnType.pooledStr) : undefStr, 
    };

    return ret;
}




async function parseInternal(data: string, ctx: Computation.Context): Promise<Result<Schema.File>> {
    const tokenizer = Tokenizer(data);

    ctx.update({ message: 'Parsing...', current: 0, max: data.length });
    const structures: Schema.Structure[] = [];
    while (tokenizer.position < data.length) {
        const state = State(tokenizer, ctx);
        handleMolecule(state);
        const atoms = await handleAtoms(state);
        const bonds = await handleBonds(state);
        structures.push({ molecule: state.molecule, atoms, bonds });
    }

    const result: Schema.File = { structures };
    return Result.success(result);
}





export function parse(data: string) {
    return Computation.create<Result<Schema.File>>(async ctx => {
        return await parseInternal(data, ctx);
    });
}

export default parse;