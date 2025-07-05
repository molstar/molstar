/**
 * Copyright (c) 2020-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { Spacegroup, SpacegroupCell } from '../../mol-math/geometry';
import { Vec3 } from '../../mol-math/linear-algebra';
import { getChainId } from './common/util';

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
    // TODO handle `crysin.setting`
    if (crysin.setting !== 1) return;

    const spaceCell = SpacegroupCell.create(
        crysin.spaceGroup,
        Vec3.create(crysin.a, crysin.b, crysin.c),
        Vec3.scale(Vec3(), Vec3.create(crysin.alpha, crysin.beta, crysin.gamma), Math.PI / 180)
    );

    return {
        spacegroup: Spacegroup.create(spaceCell),
        assemblies: [],
        isNonStandardCrystalFrame: false,
        ncsOperators: []
    };
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
