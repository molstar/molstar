import Tokenizer from '../common/text/tokenizer'
import FixedColumn from '../common/text/column/fixed'
import { ColumnType, UndefinedColumn } from '../common/column'
import * as Schema from './schema'
import Result from '../result'
//import Computation from '../../utils/computation'

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

/**
 * title string (free format string, optional time in ps after 't=')
 */
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

/**
 * This format is fixed, ie. all columns are in a fixed position.
 * Optionally (for now only yet with trjconv) you can write gro files
 * with any number of decimal places, the format will then be n+5
 * positions with n decimal places (n+1 for velocities) in stead
 * of 8 with 3 (with 4 for velocities). Upon reading, the precision
 * will be inferred from the distance between the decimal points
 * (which will be n+5). Columns contain the following information
 * (from left to right):
 *     residue number (5 positions, integer)
 *     residue name (5 characters)
 *     atom name (5 characters)
 *     atom number (5 positions, integer)
 *     position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
 *     velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)
 */
function handleAtoms(state: State): Schema.Atoms {
    const { tokenizer, molecule } = state;

    ////////// not using readLinesAsync /////////
    const lines =  Tokenizer.readLines(tokenizer, molecule.num_atoms);

    const pO = 20;
    const pW = state.header.precision.position + 5;
    const vO = pO + 3 * pW;
    const vW = state.header.precision.velocity + 4;

    const col = FixedColumn({ data: tokenizer.data, lines, rowCount: state.numberOfAtoms });

    const ret = {
        count: molecule.num_atoms,
        atom_id: col(0, 0, ColumnType.int),
        atom_name: col(0, 0, ColumnType.str),
        x: col(0, 0, ColumnType.float),
    };

    return ret;
}

/**
 * box vectors (free format, space separated reals), values:
 * v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y),
 * the last 6 values may be omitted (they will be set to zero).
 * Gromacs only supports boxes with v1(y)=v1(z)=v2(z)=0.
 */
function handleBoxVectors(state: State) {
    const { tokenizer } = state;
    markLine(tokenizer);
    const values = getTokenString(tokenizer).trim().split(/\s+/g);
    state.header.box = [+values[0], +values[1], +values[2]];
}

function parseInternal(data: string): Result<Schema.File> {
    const tokenizer = TokenizerState(data);

    const structures: Schema.Structure[] = [];
    while (tokenizer.position < data.length) {
        const state = createState(tokenizer);
        handleMolecule(state);
        const atoms = handleAtoms(state);
        handleBoxVectors(state);
        structures.push({ header: state.header, atoms });
    }

    const result: Schema.File = { structures };
    return Result.success(result);
}

export function parse(data: string) {
    return parseInternal(data);
}

export default parse;