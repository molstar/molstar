/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Augment Agent
 */

import { Column, Table } from '../../mol-data/db';
import { Spacegroup, SpacegroupCell } from '../../mol-math/geometry';
import { Vec3 } from '../../mol-math/linear-algebra';
import { ZmlFile } from '../../mol-io/reader/zml/schema';
import { Trajectory } from '../../mol-model/structure';
import { BondType, getElementFromAtomicNumber, getMoleculeType } from '../../mol-model/structure/model/types';
import { RuntimeContext, Task } from '../../mol-task';
import { ModelFormat } from '../format';
import { createModels } from './basic/parser';
import { BasicSchema, createBasic } from './basic/schema';
import { ComponentBuilder } from './common/component';
import { EntityBuilder } from './common/entity';
import { getChainId } from './common/util';
import { IndexPairBonds } from './property/bonds/index-pair';
import { ModelSymmetry } from './property/symmetry';

async function getModels(zml: ZmlFile, ctx: RuntimeContext): Promise<Trajectory> {
    const { molsys, positions, atomCount } = zml;

    const type_symbols = new Array<string>(atomCount);
    const atom_ids = new Int32Array(atomCount);
    const atom_names = new Array<string>(atomCount);
    const comp_ids = new Array<string>(atomCount);
    const seq_ids = new Int32Array(atomCount);
    const ins_codes = new Array<string>(atomCount);
    const asym_ids = new Array<string>(atomCount);
    const entity_ids = new Array<string>(atomCount);
    const formal_charges = new Int32Array(atomCount);
    const x = new Float32Array(atomCount);
    const y = new Float32Array(atomCount);
    const z = new Float32Array(atomCount);

    const entityBuilder = new EntityBuilder();
    const seqIdColumn = Column.ofIntArray(seq_ids);
    const atomNameColumn = Column.ofStringArray(atom_names);
    const componentBuilder = new ComponentBuilder(seqIdColumn, atomNameColumn);

    let atomIdx = 0;
    for (let m = 0; m < molsys.molecules.length; m++) {
        const mol = molsys.molecules[m];
        const asymId = mol.chain_id || getChainId(m);
        for (const res of mol.residues) {
            const compId = res.name || 'UNK';
            for (const atom of res.atoms) {
                atom_ids[atomIdx] = atom.id;
                atom_names[atomIdx] = atom.name;
                type_symbols[atomIdx] = getElementFromAtomicNumber(atom.atomic_number);
                comp_ids[atomIdx] = compId;
                seq_ids[atomIdx] = res.id;
                ins_codes[atomIdx] = res.insertion_code || '';
                asym_ids[atomIdx] = asymId;
                formal_charges[atomIdx] = atom.formal_charge;
                x[atomIdx] = positions[atomIdx * 3 + 0];
                y[atomIdx] = positions[atomIdx * 3 + 1];
                z[atomIdx] = positions[atomIdx * 3 + 2];
                atomIdx++;
            }
        }
    }

    atomIdx = 0;
    for (const mol of molsys.molecules) {
        for (const res of mol.residues) {
            const compId = res.name || 'UNK';
            const moleculeType = getMoleculeType(componentBuilder.add(compId, atomIdx).type, compId);
            const entityId = entityBuilder.getEntityId(compId, moleculeType, asym_ids[atomIdx]);
            for (let i = 0; i < res.atoms.length; i++) {
                entity_ids[atomIdx + i] = entityId;
            }
            atomIdx += res.atoms.length;
        }
    }

    const type_symbol = Column.ofStringArray(type_symbols);
    const asym_id = Column.ofStringArray(asym_ids);
    const comp_id = Column.ofStringArray(comp_ids);
    const entity_id = Column.ofStringArray(entity_ids);
    const ins_code = Column.ofStringArray(ins_codes);

    const atom_site = Table.ofPartialColumns(BasicSchema.atom_site, {
        auth_asym_id: asym_id,
        auth_atom_id: atomNameColumn,
        auth_comp_id: comp_id,
        auth_seq_id: seqIdColumn,
        Cartn_x: Column.ofFloatArray(x),
        Cartn_y: Column.ofFloatArray(y),
        Cartn_z: Column.ofFloatArray(z),
        id: Column.ofIntArray(atom_ids),
        pdbx_formal_charge: Column.ofIntArray(formal_charges),
        pdbx_PDB_ins_code: ins_code,

        label_asym_id: asym_id,
        label_atom_id: atomNameColumn,
        label_comp_id: comp_id,
        label_seq_id: seqIdColumn,
        label_entity_id: entity_id,

        occupancy: Column.ofConst(1, atomCount, Column.Schema.float),
        type_symbol,

        pdbx_PDB_model_num: Column.ofConst(1, atomCount, Column.Schema.int),
    }, atomCount);

    const basic = createBasic({
        entity: entityBuilder.getEntityTable(),
        chem_comp: componentBuilder.getChemCompTable(),
        atom_site,
    });

    const models = await createModels(basic, ZmlFormat.create(zml), ctx);
    if (models.frameCount === 0) return models;

    const bondCount = molsys.bonds.length;
    if (bondCount > 0) {
        const indexA = new Int32Array(bondCount);
        const indexB = new Int32Array(bondCount);
        const order = new Int8Array(bondCount);
        const flag = new Int8Array(bondCount);
        for (let i = 0; i < bondCount; i++) {
            const [a, b, o] = molsys.bonds[i];
            indexA[i] = a;
            indexB[i] = b;
            order[i] = o;
            flag[i] = BondType.Flag.Covalent;
        }
        const pairBonds = IndexPairBonds.fromData(
            { pairs: { indexA: Column.ofIntArray(indexA), indexB: Column.ofIntArray(indexB), order: Column.ofIntArray(order), flag: Column.ofIntArray(flag) }, count: bondCount },
            { maxDistance: molsys.box ? -1 : Infinity }
        );
        IndexPairBonds.Provider.set(models.representative, pairBonds);
    }

    if (molsys.box) {
        const symmetry = getSymmetry(molsys.box);
        if (symmetry) ModelSymmetry.Provider.set(models.representative, symmetry);
    }

    return models;
}

function getSymmetry(box: ReadonlyArray<ReadonlyArray<number>>) {
    if (box.length !== 3) return undefined;
    const ax = Vec3.fromArray(Vec3(), box[0] as ArrayLike<number> as number[], 0);
    const bx = Vec3.fromArray(Vec3(), box[1] as ArrayLike<number> as number[], 0);
    const cx = Vec3.fromArray(Vec3(), box[2] as ArrayLike<number> as number[], 0);
    const a = Vec3.magnitude(ax), b = Vec3.magnitude(bx), c = Vec3.magnitude(cx);
    if (a <= 0 || b <= 0 || c <= 0) return undefined;
    const alpha = Math.acos(Math.max(-1, Math.min(1, Vec3.dot(bx, cx) / (b * c))));
    const beta = Math.acos(Math.max(-1, Math.min(1, Vec3.dot(ax, cx) / (a * c))));
    const gamma = Math.acos(Math.max(-1, Math.min(1, Vec3.dot(ax, bx) / (a * b))));
    const cell = SpacegroupCell.create('P 1', Vec3.create(a, b, c), Vec3.create(alpha, beta, gamma));
    return { spacegroup: Spacegroup.create(cell), assemblies: [], isNonStandardCrystalFrame: false, ncsOperators: [] };
}

export { ZmlFormat };
type ZmlFormat = ModelFormat<ZmlFile>
namespace ZmlFormat {
    export function is(x?: ModelFormat): x is ZmlFormat { return x?.kind === 'zml'; }
    export function create(zml: ZmlFile): ZmlFormat { return { kind: 'zml', name: zml.name, data: zml }; }
}

export function trajectoryFromZml(zml: ZmlFile): Task<Trajectory> {
    return Task.create('Parse ZML', ctx => getModels(zml, ctx));
}
