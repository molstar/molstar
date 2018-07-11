/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Segmentation } from 'mol-data/int';
import { CifWriter } from 'mol-io/writer/cif';
import { SecondaryStructure } from '../../model/properties/seconday-structure';
import { StructureElement, Unit, StructureProperties as P } from '../../structure';
import { CifExportContext } from '../mmcif';
import CifField = CifWriter.Field
import CifCategory = CifWriter.Category
import { Column } from 'mol-data/db';

export function _struct_conf(ctx: CifExportContext): CifCategory {
    const elements = findElements(ctx, 'helix');
    return {
        data: elements,
        name: 'struct_conf',
        fields: struct_conf_fields,
        rowCount: elements.length
    };
}

export function _struct_sheet_range(ctx: CifExportContext): CifCategory {
    const elements = (findElements(ctx, 'sheet') as SSElement<SecondaryStructure.Sheet>[]).sort(compare_ssr);
    return {
        data: elements,
        name: 'struct_sheet_range',
        fields: struct_sheet_range_fields,
        rowCount: elements.length
    };
}

function compare_ssr(x: SSElement<SecondaryStructure.Sheet>, y: SSElement<SecondaryStructure.Sheet>) {
    const a = x.element, b = y.element;
    return a.sheet_id < b.sheet_id ? -1 : a.sheet_id === b.sheet_id ? x.start.element - y.start.element : 1
};

const struct_conf_fields: CifField[] = [
    CifField.str<number, SSElement<SecondaryStructure.Helix>[]>('conf_type_id', (i, data) => data[i].element.type_id),
    CifField.str<number, SSElement<SecondaryStructure.Helix>[]>('id', (i, data, idx) => `${data[i].element.type_id}${idx + 1}`),
    ...residueIdFields('beg_', e => e.start),
    ...residueIdFields('end_', e => e.end),
    CifField.str<number, SSElement<SecondaryStructure.Helix>[]>('pdbx_PDB_helix_class', (i, data) => data[i].element.helix_class),
    CifField.str<number, SSElement<SecondaryStructure.Helix>[]>('details', (i, data) => data[i].element.details || '', {
        valueKind: (i, d) => !!d[i].element.details ? Column.ValueKind.Present : Column.ValueKind.Unknown
    }),
    CifField.int<number, SSElement<SecondaryStructure.Helix>[]>('pdbx_PDB_helix_class', (i, data) => data[i].length)
];

const struct_sheet_range_fields: CifField[] = [
    CifField.str<number, SSElement<SecondaryStructure.Sheet>[]>('sheet_id', (i, data) => data[i].element.sheet_id),
    CifField.index('id'),
    ...residueIdFields('beg_', e => e.start),
    ...residueIdFields('end_', e => e.end),
    CifField.str('symmetry', (i, data) => '', { valueKind: (i, d) => Column.ValueKind.Unknown })
];

function residueIdFields(prefix: string, loc: (e: SSElement<any>) => StructureElement): CifField<number, SSElement<SecondaryStructure.Helix>[]>[] {
    return [
        CifField.str(`${prefix}label_comp_id`, (i, d) => P.residue.label_comp_id(loc(d[i]))),
        CifField.int(`${prefix}label_seq_id`, (i, d) => P.residue.label_seq_id(loc(d[i]))),
        CifField.str(`pdbx_${prefix}PDB_ins_code`, (i, d) => P.residue.pdbx_PDB_ins_code(loc(d[i]))),
        CifField.str(`${prefix}label_asym_id`, (i, d) => P.chain.label_asym_id(loc(d[i]))),
        CifField.str(`${prefix}label_entity_id`, (i, d) => P.chain.label_entity_id(loc(d[i]))),
        CifField.str(`${prefix}auth_comp_id`, (i, d) => P.residue.auth_comp_id(loc(d[i]))),
        CifField.int(`${prefix}auth_seq_id`, (i, d) => P.residue.auth_seq_id(loc(d[i]))),
        CifField.str(`${prefix}auth_asym_id`, (i, d) => P.chain.auth_asym_id(loc(d[i])))
    ];
}

interface SSElement<T extends SecondaryStructure.Element> {
    start: StructureElement,
    end: StructureElement,
    length: number,
    element: T
}

function findElements<T extends SecondaryStructure.Element>(ctx: CifExportContext, kind: SecondaryStructure.Element['kind']) {
    const { key, elements } = ctx.model.properties.secondaryStructure;

    const ssElements: SSElement<any>[] = [];

    for (const unit of ctx.structure.units) {
        // currently can only support this for "identity" operators.
        if (!Unit.isAtomic(unit) || !unit.conformation.operator.isIdentity) continue;

        const segs = unit.model.atomicHierarchy.residueAtomSegments;
        const residues = Segmentation.transientSegments(segs, unit.elements);

        let current: Segmentation.Segment, move = true;
        while (residues.hasNext) {
            if (move) current = residues.move();

            const start = current!.index;
            const startIdx = key[start];
            const element = elements[startIdx];
            if (element.kind !== kind) {
                move = true;
                continue;
            }

            let prev = start;
            while (residues.hasNext) {
                prev = current!.index;
                current = residues.move();
                if (startIdx !== key[current.index]) {
                    move = false;
                    ssElements[ssElements.length] = {
                        start: StructureElement.create(unit, segs.offsets[start]),
                        end: StructureElement.create(unit, segs.offsets[prev]),
                        length: prev - start + 1,
                        element
                    }
                    break;
                }
            }
        }
    }

    return ssElements as SSElement<T>[];
}