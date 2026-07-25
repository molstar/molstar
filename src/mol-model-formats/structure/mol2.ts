/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column, Table } from '../../mol-data/db';
import { Model, Symmetry } from '../../mol-model/structure/model';
import { BondType, getMoleculeType, MoleculeType } from '../../mol-model/structure/model/types';
import { RuntimeContext, Task } from '../../mol-task';
import { createModels } from './basic/parser';
import { BasicSchema, createBasic } from './basic/schema';
import { ComponentBuilder } from './common/component';
import { EntityBuilder } from './common/entity';
import { ModelFormat } from '../format';
import { IndexPairBonds } from './property/bonds/index-pair';
import { Mol2Crysin, Mol2File } from '../../mol-io/reader/mol2/schema';
import { AtomPartialCharge } from './property/partial-charge';
import { Trajectory, ArrayTrajectory } from '../../mol-model/structure';
import { guessElementSymbolString } from './util';
import { ModelSymmetry } from './property/symmetry';
import { Vec3 } from '../../mol-math/linear-algebra';
import { getChainId } from './common/util';
import { Spacegroup, SpacegroupCell } from '../../mol-math/geometry/spacegroup/construction';
import { AxisPermutations, TetragonalCenteringBasisop, transformOperators, transformTetragonalCentering } from '../../mol-math/geometry/spacegroup/common';
import { Mat4 } from '../../mol-math/linear-algebra/3d/mat4';

type Subst = { chain: string, subType: string }

async function getModels(mol2: Mol2File, ctx: RuntimeContext) {
    const models: Model[] = [];

    for (let i = 0, il = mol2.structures.length; i < il; ++i) {
        const { molecule, atoms, bonds, substructures, crysin } = mol2.structures[i];

        const entityIds = new Array<string>(atoms.count);
        const asymIds = new Array<string>(atoms.count);
        const typeSymbols = new Array<string>(atoms.count);
        const compIds = new Array<string>(atoms.count);
        const seqIds = new Array<number>(atoms.count);

        let hasAtomType = false;
        for (let i = 0; i < atoms.count; ++i) {
            if (atoms.atom_type.value(i).includes('.')) {
                hasAtomType = true;
                break;
            }
        }

        const substByRoot = new Map<number, Subst>();
        if (substructures) {
            for (let i = 0; i < substructures.count; ++i) {
                const rootAtom = substructures.root_atom.value(i);
                const chain = substructures.chain.value(i);
                const subType = substructures.sub_type.value(i);
                substByRoot.set(rootAtom, { chain, subType });
            }
        }

        const entityBuilder = new EntityBuilder();
        const componentBuilder = new ComponentBuilder(atoms.subst_id, atoms.atom_name);

        let currentSubst: Subst | undefined = undefined;
        let currentEntityId = '';
        let currentAsymIndex = 0;
        let currentAsymId = '';
        let prevMoleculeType = MoleculeType.Unknown;
        let prevSubstId = -1;

        for (let i = 0; i < atoms.count; ++i) {
            const substName = atoms.subst_name.value(i);
            const substId = atoms.subst_id.value(i);

            let compId = substName.replace(/\d+$/, '') || 'MOL';

            const seqIdMatch = substName.match(/(\d+)$/);
            seqIds[i] = seqIdMatch ? parseInt(seqIdMatch[1]) : substId;

            const subst = substByRoot.get(atoms.atom_id.value(i));
            if (subst) currentSubst = subst;

            if (currentSubst) {
                if (currentSubst.subType && currentSubst.subType !== 'UNK') {
                    compId = currentSubst.subType;
                }
            }

            typeSymbols[i] = hasAtomType
                ? atoms.atom_type.value(i).split('.')[0].toUpperCase()
                : guessElementSymbolString(atoms.atom_name.value(i), compId);

            compIds[i] = compId;

            if (substId !== prevSubstId) {
                const moleculeType = getMoleculeType(componentBuilder.add(compId, i).type, compId);

                if (currentSubst?.chain) {
                    currentAsymId = currentSubst.chain;
                } else {
                    if (moleculeType !== prevMoleculeType || substId !== prevSubstId + 1) {
                        currentAsymId = getChainId(currentAsymIndex);
                        currentAsymIndex += 1;
                    }
                }

                currentEntityId = entityBuilder.getEntityId(compId, moleculeType, currentAsymId);

                prevSubstId = substId;
                prevMoleculeType = moleculeType;
            }

            entityIds[i] = currentEntityId;
            asymIds[i] = currentAsymId;
        }

        const type_symbol = Column.ofStringArray(typeSymbols);
        const comp_id = Column.ofStringArray(compIds);
        const entity_id = Column.ofStringArray(entityIds);
        const asym_id = Column.ofStringArray(asymIds);
        const atom_id = Column.asArrayColumn(atoms.atom_name);
        const seq_id = Column.ofIntArray(seqIds);

        const atom_site = Table.ofPartialColumns(BasicSchema.atom_site, {
            auth_asym_id: asym_id,
            auth_atom_id: atom_id,
            auth_comp_id: comp_id,
            auth_seq_id: seq_id,
            Cartn_x: Column.asArrayColumn(atoms.x, Float32Array),
            Cartn_y: Column.asArrayColumn(atoms.y, Float32Array),
            Cartn_z: Column.asArrayColumn(atoms.z, Float32Array),
            id: Column.asArrayColumn(atoms.atom_id),

            label_asym_id: asym_id,
            label_atom_id: atom_id,
            label_comp_id: comp_id,
            label_seq_id: atoms.subst_id,
            label_entity_id: entity_id,

            occupancy: Column.ofConst(1, atoms.count, Column.Schema.float),
            type_symbol,

            pdbx_PDB_model_num: Column.ofConst(i, atoms.count, Column.Schema.int),
        }, atoms.count);

        const basic = createBasic({
            entity: entityBuilder.getEntityTable(),
            chem_comp: componentBuilder.getChemCompTable(),
            atom_site
        });

        const _models = await createModels(basic, Mol2Format.create(mol2), ctx);

        if (_models.frameCount > 0) {
            const indexA = Column.ofIntArray(Column.mapToArray(bonds.origin_atom_id, x => x - 1, Int32Array));
            const indexB = Column.ofIntArray(Column.mapToArray(bonds.target_atom_id, x => x - 1, Int32Array));
            const key = bonds.bond_id;
            const order = Column.ofIntArray(Column.mapToArray(bonds.bond_type, x => {
                switch (x) {
                    case 'ar': // aromatic
                    case 'am': // amide
                    case 'un': // unknown
                        return 1;
                    case 'du': // dummy
                    case 'nc': // not connected
                        return 0;
                    default:
                        return parseInt(x);
                }
            }, Int8Array));
            const flag = Column.ofIntArray(Column.mapToArray(bonds.bond_type, x => {
                switch (x) {
                    case 'ar': // aromatic
                    case 'am': // amide
                        return BondType.Flag.Aromatic | BondType.Flag.Covalent;
                    case 'du': // dummy
                    case 'nc': // not connected
                        return BondType.Flag.None;
                    case 'un': // unknown
                    default:
                        return BondType.Flag.Covalent;
                }
            }, Int8Array));
            const pairBonds = IndexPairBonds.fromData(
                { pairs: { key, indexA, indexB, order, flag }, count: atoms.count },
                { maxDistance: crysin ? -1 : Infinity }
            );

            const first = _models.representative;
            IndexPairBonds.Provider.set(first, pairBonds);

            AtomPartialCharge.Provider.set(first, {
                data: atoms.charge,
                type: molecule.charge_type
            });

            if (crysin) {
                const symmetry = getSymmetry(crysin);
                if (symmetry) ModelSymmetry.Provider.set(first, symmetry);
            }

            models.push(first);
        }
    }

    return new ArrayTrajectory(models);
}

function getSymmetry(crysin: Mol2Crysin): Symmetry | undefined {
    const spaceCell = SpacegroupCell.create(
        crysin.spaceGroup,
        Vec3.create(crysin.a, crysin.b, crysin.c),
        Vec3.scale(Vec3(), Vec3.create(crysin.alpha, crysin.beta, crysin.gamma), Math.PI / 180)
    );

    return {
        spacegroup: getMol2CrysinSpacegroup(spaceCell, crysin.spaceGroup, crysin.setting),
        assemblies: [],
        isNonStandardCrystalFrame: false,
        ncsOperators: []
    };
}

//

/**
 * A non-default CRYSIN setting: its Hermann-Mauguin name and the recipe for its
 * symmetry operators. `axisPermutation >= 0` indexes `AxisPermutations` (a pure
 * axis relabelling of the canonical monoclinic/orthorhombic group); the sentinel
 * `-1` marks the tetragonal C/F "setting 2" descriptions, which are a genuine
 * cell doubling handled by `transformTetragonalCentering` (below) - this is a
 * CRYSIN/SYBYL-manual-specific convention, not something International Tables
 * catalogs as a numbered per-group setting.
 */
type SettingEntry = readonly [name: string, axisPermutation: number];

/**
 * The single source of non-default CRYSIN setting (2-6) data, keyed by ITA
 * space-group number. Each row is indexed by `[setting - 2]`; `null` marks a
 * setting with no alternate and trailing empties are omitted.
 */
export const SettingsByNumber: { [spacegroupNumber: number]: (SettingEntry | null)[] } = {
    3: [null, ['P 1 1 2', 2]],
    4: [null, ['P 1 1 21', 2]],
    5: [['A 2', 3], ['B 1 1 2', 5], ['A 1 1 2', 2]],
    6: [null, ['P 1 1 m', 2]],
    7: [['P a', 3], ['P 1 1 b', 5], ['P 1 1 a', 2]],
    8: [['A m', 3], ['B 1 1 m', 5], ['A 1 1 m', 2]],
    9: [['A a', 3], ['B 1 1 b', 5], ['A 1 1 a', 2]],
    10: [null, ['P 1 1 2/m', 2]],
    11: [null, ['P 1 1 21/m', 2]],
    12: [['A 2/m', 3], ['B 1 1 2/m', 5], ['A 1 1 2/m', 2]],
    13: [['P 2/a', 3], ['P 1 1 2/b', 5], ['P 1 1 2/a', 2]],
    14: [['P 21/a', 3], ['P 1 1 21/b', 5], ['P 1 1 21/a', 2]],
    15: [['A 2/a', 3], ['B 1 1 2/b', 5], ['A 1 1 2/a', 2]],
    17: [['P 21 2 2', 2], ['P 2 21 2', 4]],
    18: [['P 2 21 21', 2], ['P 21 2 21', 4]],
    20: [['A 21 2 2', 2], ['B 2 21 2', 4]],
    21: [['A 2 2 2', 2], ['B 2 2 2', 4]],
    25: [['P 2 m m', 2], ['P m 2 m', 4]],
    26: [['P 21 m a', 2], ['P b 21 m', 4], ['P m 21 b', 5], ['P c m 21', 1], ['P 21 a m', 3]],
    27: [['P 2 a a', 2], ['P b 2 b', 4]],
    28: [['P 2 m b', 2], ['P c 2 m', 4], ['P m 2 a', 5], ['P b m 2', 1], ['P 2 c m', 3]],
    29: [['P 21 a b', 2], ['P c 21 b', 4], ['P b 21 a', 5], ['P b c 21', 1], ['P 21 c a', 3]],
    30: [['P 2 n a', 2], ['P b 2 n', 4], ['P n 2 b', 5], ['P c n 2', 1], ['P 2 a n', 3]],
    31: [['P 21 m n', 2], ['P n 21 m', 4], ['P m 21 n', 5], ['P n m 21', 1], ['P 21 n m', 3]],
    32: [['P 2 c b', 2], ['P c 2 a', 4]],
    33: [['P 21 n b', 2], ['P c 21 n', 4], ['P n 21 a', 5], ['P b n 21', 1], ['P 21 c n', 3]],
    34: [['P 2 n n', 2], ['P n 2 n', 4]],
    35: [['A 2 m m', 2], ['B m 2 m', 4]],
    36: [['A 21 m a', 2], ['B b 21 m', 4], ['B m 21 b', 5], ['C c m 21', 1], ['A 21 a m', 3]],
    37: [['A 2 a a', 2], ['B b 2 b', 4]],
    38: [['B 2 m m', 2], ['C m 2 m', 4], ['A m 2 m', 5], ['B m m 2', 1], ['C 2 m m', 3]],
    39: [['B 2 c m', 2], ['C m 2 a', 4], ['A c 2 m', 5], ['B m a 2', 1], ['C 2 m b', 3]],
    40: [['B 2 m b', 2], ['C c 2 m', 4], ['A m 2 a', 5], ['B b m 2', 1], ['C 2 c m', 3]],
    41: [['B 2 c b', 2], ['C c 2 a', 4], ['A c 2 a', 5], ['B b a 2', 1], ['C 2 c b', 3]],
    42: [['F 2 m m', 2], ['F m 2 m', 4]],
    43: [['F 2 d d', 2], ['F d 2 d', 4]],
    44: [['I 2 m m', 2], ['I m 2 m', 4]],
    45: [['I 2 c b', 2], ['I c 2 a', 4]],
    46: [['I 2 m b', 2], ['I c 2 m', 4], ['I m 2 a', 5], ['I b m 2', 1], ['I 2 c m', 3]],
    49: [['P m a a', 2], ['P b m b', 4]],
    50: [['P n c b', 2], ['P c n a', 4]],
    51: [['P b m m', 2], ['P m c m', 4], ['P m a m', 5], ['P m m b', 1], ['P c m m', 3]],
    52: [['P b n n', 2], ['P n c n', 4], ['P n a n', 5], ['P n n b', 1], ['P c n n', 3]],
    53: [['P b m n', 2], ['P n c m', 4], ['P m a n', 5], ['P n m b', 1], ['P c n m', 3]],
    54: [['P b a a', 2], ['P b c b', 4], ['P b a b', 5], ['P c c b', 1], ['P c a a', 3]],
    55: [['P m c b', 2], ['P c m a', 4]],
    56: [['P n a a', 2], ['P b n b', 4]],
    57: [['P m c a', 2], ['P b m a', 4], ['P c m b', 5], ['P c a m', 1], ['P m a b', 3]],
    58: [['P m n n', 2], ['P n m n', 4]],
    59: [['P n m m', 2], ['P m n m', 4]],
    60: [['P n c a', 2], ['P b n a', 4], ['P c n b', 5], ['P c a n', 1], ['P n a b', 3]],
    61: [null, ['P c a b', 1]],
    62: [['P b n m', 2], ['P m c n', 4], ['P n a m', 5], ['P m n b', 1], ['P c m n', 3]],
    63: [['A m m a', 2], ['B b m m', 4], ['B m m b', 5], ['C c m m', 1], ['A m a m', 3]],
    64: [['A b m a', 2], ['B b c m', 4], ['B m a b', 5], ['C c m b', 1], ['A c a m', 3]],
    65: [['A m m m', 2], ['B m m m', 4]],
    66: [['A m a a', 2], ['B b m b', 4]],
    67: [['A b m m', 2], ['B m c m', 4], ['B m a m', 5], ['C m m b', 1], ['A c m m', 3]],
    68: [['A b a a', 2], ['B b c b', 4], ['B b a b', 4], ['C c c b', 0], ['A c a a', 2]],
    72: [['I m c b', 2], ['I c m a', 4]],
    73: [null, ['I c a b', 1]],
    74: [['I b m m', 2], ['I m c m', 4], ['I m a m', 5], ['I m m b', 1], ['I c m m', 3]],
    75: [['C 4', -1]],
    76: [['C 41', -1]],
    77: [['C 42', -1]],
    78: [['C 43', -1]],
    79: [['F 4', -1]],
    80: [['F 41', -1]],
    81: [['C -4', -1]],
    82: [['F -4', -1]],
    83: [['C 4/m', -1]],
    84: [['C 42/m', -1]],
    85: [['C 4/a', -1]],
    86: [['C 42/a', -1]],
    87: [['F 4/m', -1]],
    88: [['F 41/d', -1]],
    89: [['C 4 2 2', -1]],
    90: [['C 4 2 21', -1]],
    91: [['C 41 2 2', -1]],
    92: [['C 41 2 21', -1]],
    93: [['C 42 2 2', -1]],
    94: [['C 42 2 21', -1]],
    95: [['C 43 2 2', -1]],
    96: [['C 43 2 21', -1]],
    97: [['F 4 2 2', -1]],
    98: [['F 41 2 2', -1]],
    99: [['C 4 m m', -1]],
    100: [['C 4 m b', -1]],
    101: [['C 42 m c', -1]],
    102: [['C 42 m n', -1]],
    103: [['C 4 c c', -1]],
    104: [['C 4 c n', -1]],
    105: [['C 42 c m', -1]],
    106: [['C 42 c b', -1]],
    107: [['F 4 m m', -1]],
    108: [['F 4 m c', -1]],
    109: [['F 41 d m', -1]],
    110: [['F 41 d c', -1]],
    111: [['C -4 m 2', -1]],
    112: [['C -4 c 2', -1]],
    113: [['C -4 m 21', -1]],
    114: [['C -4 c 21', -1]],
    115: [['C -4 2 m', -1]],
    116: [['C -4 2 c', -1]],
    117: [['C -4 2 b', -1]],
    118: [['C -4 2 n', -1]],
    119: [['F -4 2 m', -1]],
    120: [['F -4 2 c', -1]],
    121: [['F -4 m 2', -1]],
    122: [['F -4 d 2', -1]],
    123: [['C 4/m m m', -1]],
    124: [['C 4/m c c', -1]],
    125: [['C 4/a m b', -1]],
    126: [['C 4/a c n', -1]],
    127: [['C 4/m m b', -1]],
    128: [['C 4/m c n', -1]],
    129: [['C 4/a m m', -1]],
    130: [['C 4/a c c', -1]],
    131: [['C 42/m c m', -1]],
    132: [['C 42/m m c', -1]],
    133: [['C 42/a c b', -1]],
    134: [['C 42/a m n', -1]],
    135: [['C 42/m c b', -1]],
    136: [['C 42/m m n', -1]],
    137: [['C 42/a c m', -1]],
    138: [['C 42/a m c', -1]],
    139: [['F 4/m m m', -1]],
    140: [['F 4/m m c', -1]],
    141: [['F 41/d d m', -1]],
    142: [['F 41/d d c', -1]],
};

/**
 * Computes the alternate-setting name, ITA change-of-basis symbol, and
 * symmetry operators for a MOL2 CRYSIN-style (spacegroupNumber, setting)
 * pair. Monoclinic and orthorhombic settings are a pure axis relabelling of
 * the canonical operators (via the entry's `AxisPermutations` index and
 * `transformOperators`); the tetragonal "setting 2" C/F-centered descriptions
 * are a genuine cell change derived from `baseOperators` (the canonical
 * setting-1 operators). Returns `undefined` when `setting === 1` or the
 * requested setting is unsupported - callers should fall back to the
 * canonical description in that case.
 */
export function getOperatorsForSetting(spacegroupNumber: number, setting: number, baseOperators: ReadonlyArray<Mat4>): { name: string, basisop: string, operators: Mat4[] } | undefined {
    if (setting === 1) return undefined;

    const entry = SettingsByNumber[spacegroupNumber]?.[setting - 2];
    if (!entry) return undefined;

    const [name, axisPermutation] = entry;
    if (axisPermutation >= 0) {
        const { basisop, matrix } = AxisPermutations[axisPermutation];
        return { name, basisop, operators: transformOperators(baseOperators, matrix) };
    }
    // tetragonal C/F setting 2: a genuine cell change, not an axis permutation
    return { name, basisop: TetragonalCenteringBasisop, operators: transformTetragonalCentering(baseOperators) };
}

/**
 * Resolves a MOL2 `CRYSIN` record's `(spaceGroup, setting)` pair (see
 * `Mol2Crysin` in `mol-io/reader/mol2/schema`) into a full `Spacegroup`.
 * Falls back to the canonical `Spacegroup.create(cell)` result when `setting`
 * is 1, or when the requested (spacegroupNumber, setting) combination is
 * unsupported (blank table entry, or no applicable transform) - this avoids
 * ever applying incorrect operators to a cell measured under a different
 * setting.
 */
export function getMol2CrysinSpacegroup(cell: SpacegroupCell, spacegroupNumber: number, setting: number): Spacegroup {
    const canonical = Spacegroup.create(cell);
    if (setting === 1) return canonical;

    const result = getOperatorsForSetting(spacegroupNumber, setting, canonical.operators);
    if (!result) {
        console.warn(`Unsupported spacegroup setting '${setting}' for spacegroup number ${spacegroupNumber}, using the canonical setting instead`);
        return canonical;
    }

    return Spacegroup.createWithSetting(cell, result.name, result.basisop, result.operators);
}

//

export { Mol2Format };

type Mol2Format = ModelFormat<Mol2File>

namespace Mol2Format {
    export function is(x?: ModelFormat): x is Mol2Format {
        return x?.kind === 'mol2';
    }

    export function create(mol2: Mol2File): Mol2Format {
        return { kind: 'mol2', name: mol2.name, data: mol2 };
    }
}

export function trajectoryFromMol2(mol2: Mol2File): Task<Trajectory> {
    return Task.create('Parse MOL2', ctx => getModels(mol2, ctx));
}
