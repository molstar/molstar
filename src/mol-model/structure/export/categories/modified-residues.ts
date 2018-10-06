/**
 * Copyright (c) 2017-2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Segmentation } from 'mol-data/int';
import { CifWriter } from 'mol-io/writer/cif';
import { StructureElement, StructureProperties as P, Unit } from '../../../structure';
import { CifExportContext } from '../mmcif';

import CifField = CifWriter.Field
import CifCategory = CifWriter.Category

const pdbx_struct_mod_residue_fields: CifField<number, StructureElement[]>[] = [
    CifField.index('id'),
    CifField.str(`label_comp_id`, (i, xs) => P.residue.label_comp_id(xs[i])),
    CifField.int(`label_seq_id`, (i, xs) => P.residue.label_seq_id(xs[i])),
    CifField.str(`pdbx_PDB_ins_code`, (i, xs) => P.residue.pdbx_PDB_ins_code(xs[i])),
    CifField.str(`label_asym_id`, (i, xs) => P.chain.label_asym_id(xs[i])),
    CifField.str(`label_entity_id`, (i, xs) => P.chain.label_entity_id(xs[i])),
    CifField.str(`auth_comp_id`, (i, xs) => P.residue.auth_comp_id(xs[i])),
    CifField.int(`auth_seq_id`, (i, xs) => P.residue.auth_seq_id(xs[i])),
    CifField.str(`auth_asym_id`, (i, xs) => P.chain.auth_asym_id(xs[i])),
    CifField.str<number, StructureElement[]>('parent_comp_id', (i, xs) => xs[i].unit.model.properties.modifiedResidues.parentId.get(P.residue.label_comp_id(xs[i]))!),
    CifField.str('details', (i, xs) => xs[i].unit.model.properties.modifiedResidues.details.get(P.residue.label_comp_id(xs[i]))!)
];

function getModifiedResidues({ structures }: CifExportContext): StructureElement[] {
    // TODO: can different models have differnt modified residues?
    const structure = structures[0], model = structure.model;
    const map = model.properties.modifiedResidues.parentId;
    if (!map.size) return [];

    const ret = [];
    const prop = P.residue.label_comp_id;
    const loc = StructureElement.create();
    for (const unit of structure.units) {
        if (!Unit.isAtomic(unit) || !unit.conformation.operator.isIdentity) continue;
        const residues = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, unit.elements);
        loc.unit = unit;
        while (residues.hasNext) {
            const seg = residues.move();
            loc.element = unit.elements[seg.start];
            const name = prop(loc);
            if (map.has(name)) {
                ret[ret.length] = StructureElement.create(loc.unit, loc.element);
            }
        }
    }
    return ret;
}

export const _pdbx_struct_mod_residue: CifCategory<CifExportContext> = {
    name: 'pdbx_struct_mod_residue',
    instance(ctx) {
        const residues = getModifiedResidues(ctx);
        return {
            fields: pdbx_struct_mod_residue_fields,
            source: [{ data: residues, rowCount: residues.length }]
        };
    }
}