import Tokenizer from '../common/text/tokenizer'
import FixedColumn from '../common/text/column/fixed'
import { ColumnType, UndefinedColumn } from '../common/column'
import * as Schema from './schema'
import Result from '../result'
import Computation from '../../utils/computation' ////////// not using this



interface State {
    tokenizer: Tokenizer,
    molecule: Schema.Molecule,
    ///////////// not including Computation.chunker /////////////
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




function State(tokenizer: Tokenizer): State { //////////// not having ctx: Computation.Context as a parameter //////////////
    return {
        tokenizer,
        molecule: createEmptyMolecule(),
        //////////// not having chunker: Computation.chunker(ctx, 100000) ///////////
    };
}





function handleMolecule(state: State) {
    const { tokenizer, molecule } = state;
    
    Tokenizer.markLine(tokenizer);
    let name = Tokenizer.getTokenString(tokenizer);
    molecule.mol_name = name;

    Tokenizer.markLine(tokenizer);
    const values = Tokenizer.getTokenString(tokenizer).trim().split(/\s+/g);
    molecule.num_atoms = parseInt(values[0]);
    molecule.num_bonds = parseInt(values[1]);
    molecule.num_subst = parseInt(values[2]);
    molecule.num_feat = parseInt(values[3]);
    molecule.num_sets = parseInt(values[4]);

    Tokenizer.markLine(tokenizer);
    molecule.mol_type = Tokenizer.getTokenString(tokenizer);

    Tokenizer.markLine(tokenizer);
    molecule.charge_type = Tokenizer.getTokenString(tokenizer);

    // skip the empty line
    Tokenizer.markLine(tokenizer)

}





function handleAtoms(state: State): Schema.Atoms {
    const { tokenizer, molecule } = state;

    ////////// not using readLinesAsync /////////
    const lines =  Tokenizer.readLines(tokenizer, molecule.num_atoms);

    // default all false
    const hasSubst_id = false;
    const hasSubst_name = false;
    const hasCharge = false;
    const hasStatus_bit = false;

    /*
    const pO = 20;
    const pW = state.header.precision.position + 5;
    const vO = pO + 3 * pW;
    const vW = state.header.precision.velocity + 4;
    */

    const col = FixedColumn(lines);
    const undefInt = UndefinedColumn(molecule.num_atoms, ColumnType.int);
    const undefFloat = UndefinedColumn(molecule.num_atoms, ColumnType.float);
    const undefStr = UndefinedColumn(molecule.num_atoms, ColumnType.str);
    const undefPooledStr = UndefinedColumn(molecule.num_atoms, ColumnType.pooledStr);

    /////// wanted to have const undef = UndefinedColumn(molecule.num_atoms) like col, but failed

    /////// some unclear about the formatting, like the field sizes
    const ret = {
        count: molecule.num_atoms,
        atom_id: col(0, 7, ColumnType.int),
        atom_name: col(7, 9, ColumnType.str), ////// don't know use str or pooledStr
        x: col(16, 10, ColumnType.float),
        y: col(26, 10, ColumnType.float),
        z: col(36, 10, ColumnType.float),
        atom_type: col(46, 0, ColumnType.str), ////// don't know use str or pooledStr //////// don't know which is the atom_type
        // optional properties
        subst_id: hasSubst_id ? col(0, 0, ColumnType.int) : undefInt, 
        subst_name: hasSubst_name ? col(0, 0, ColumnType.str) : undefStr,///////// don't know use str or pooledStr
        charge: hasCharge ? col(0, 0, ColumnType.float) : undefFloat, //////// don't know use int or float
        status_bit: hasStatus_bit ? col(0, 0, ColumnType.pooledStr) : undefPooledStr, ////////// don't know use str or pooledStr
    };

    return ret;
}




function handleBonds(state: State): Schema.Bonds {
    const { tokenizer, molecule } = state;

    ////////// not using readLinesAsync /////////
    const lines =  Tokenizer.readLines(tokenizer, molecule.num_bonds);

    // default all false
    const hasStatus_bit = false;

    /*
    const pO = 20;
    const pW = state.header.precision.position + 5;
    const vO = pO + 3 * pW;
    const vW = state.header.precision.velocity + 4;
    */

    const col = FixedColumn(lines);
    const undefInt = UndefinedColumn(molecule.num_bonds, ColumnType.int);
    const undefFloat = UndefinedColumn(molecule.num_bonds, ColumnType.float);
    const undefStr = UndefinedColumn(molecule.num_bonds, ColumnType.str);
    const undefPooledStr = UndefinedColumn(molecule.num_bonds, ColumnType.pooledStr);

    /////// wanted to have const undef = UndefinedColumn(molecule.num_atoms) like col, but failed

    /////// some unclear about the formatting, like the field sizes
    const ret = {
        count: molecule.num_bonds,
        bond_id: col(0, 6, ColumnType.int),
        origin_atom_id: col(6, 6, ColumnType.int), 
        target_atom_id: col(12, 6, ColumnType.int),
        bond_type: col(18, 5, ColumnType.str), ///////// don't know use str or pooledStr
        // optional properties
        status_bits: hasStatus_bit ? col(0, 0, ColumnType.str) : undefStr, ///////// don't know use str or pooledStr
    };

    return ret;
}




//////// not using async here
function parseInternal(data: string): Result<Schema.File> { /////// not having ctx as a parameter, and not returning Promise
    const tokenizer = Tokenizer(data);

    const structures: Schema.Structure[] = [];
    while (tokenizer.position < data.length) {
        const state = State(tokenizer);//////////different
        handleMolecule(state);
        const atoms = handleAtoms(state);
        const bonds = handleBonds(state);
        structures.push({ molecule: state.molecule, atoms, bonds });
    }

    const result: Schema.File = { structures };
    return Result.success(result);
}





///////// diffrent than gro parser
export function parse(data: string) {
    return parseInternal(data);
}

export default parse;