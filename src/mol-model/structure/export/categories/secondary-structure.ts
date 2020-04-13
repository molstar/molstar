/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Segmentation } from '../../../../mol-data/int';
import { CifWriter } from '../../../../mol-io/writer/cif';
import { SecondaryStructure } from '../../model/properties/seconday-structure';
import { StructureElement, Unit } from '../../structure';
import { CifExportContext } from '../mmcif';
import CifField = CifWriter.Field
import CifCategory = CifWriter.Category
import { Column } from '../../../../mol-data/db';
import { residueIdFields } from './atom_site';
import { ModelSecondaryStructure } from '../../../../mol-model-formats/structure/property/secondary-structure';

export const _struct_conf: CifCategory<CifExportContext> = {
    name: 'struct_conf',
    instance(ctx) {
        const elements = findElements(ctx, 'helix');
        return {
            fields: struct_conf_fields,
            source: [{ data: elements, rowCount: elements.length }]
        };
    }
};

export const _struct_sheet_range: CifCategory<CifExportContext> = {
    name: 'struct_sheet_range',
    instance(ctx) {
        const elements = (findElements(ctx, 'sheet') as SSElement<SecondaryStructure.Sheet>[]).sort(compare_ssr);
        return {
            fields: struct_sheet_range_fields,
            source: [{ data: elements, rowCount: elements.length }]
        };
    }
};

function compare_ssr(x: SSElement<SecondaryStructure.Sheet>, y: SSElement<SecondaryStructure.Sheet>) {
    const a = x.element, b = y.element;
    return a.sheet_id < b.sheet_id ? -1 : a.sheet_id === b.sheet_id ? x.start.element - y.start.element : 1;
};

const struct_conf_fields: CifField[] = [
    CifField.str<number, SSElement<SecondaryStructure.Helix>[]>('conf_type_id', (i, data) => data[i].element.type_id),
    CifField.str<number, SSElement<SecondaryStructure.Helix>[]>('id', (i, data, idx) => `${data[i].element.type_id}${idx + 1}`),
    ...residueIdFields<number, SSElement<SecondaryStructure.Helix>[]>((i, e) => e[i].start, { prefix: 'beg' }),
    ...residueIdFields<number, SSElement<SecondaryStructure.Helix>[]>((i, e) => e[i].end, { prefix: 'end' }),
    CifField.str<number, SSElement<SecondaryStructure.Helix>[]>('pdbx_PDB_helix_class', (i, data) => data[i].element.helix_class),
    CifField.str<number, SSElement<SecondaryStructure.Helix>[]>('details', (i, data) => data[i].element.details || '', {
        valueKind: (i, d) => !!d[i].element.details ? Column.ValueKind.Present : Column.ValueKind.Unknown
    }),
    CifField.int<number, SSElement<SecondaryStructure.Helix>[]>('pdbx_PDB_helix_length', (i, data) => data[i].length)
];

const struct_sheet_range_fields: CifField[] = [
    CifField.str<number, SSElement<SecondaryStructure.Sheet>[]>('sheet_id', (i, data) => data[i].element.sheet_id),
    CifField.index('id'),
    ...residueIdFields<number, SSElement<SecondaryStructure.Sheet>[]>((i, e) => e[i].start, { prefix: 'beg' }),
    ...residueIdFields<number, SSElement<SecondaryStructure.Sheet>[]>((i, e) => e[i].end, { prefix: 'end' }),
    CifField.str('symmetry', (i, data) => '', { valueKind: (i, d) => Column.ValueKind.Unknown })
];

interface SSElement<T extends SecondaryStructure.Element> {
    start: StructureElement.Location,
    end: StructureElement.Location,
    length: number,
    element: T
}

function findElements<T extends SecondaryStructure.Element>(ctx: CifExportContext, kind: SecondaryStructure.Element['kind']) {
    // TODO: encode secondary structure for different models?
    const secondaryStructure = ModelSecondaryStructure.Provider.get(ctx.firstModel);
    if (!secondaryStructure) return [] as SSElement<T>[];

    const { key, elements } = secondaryStructure;
    const ssElements: SSElement<any>[] = [];

    const structure = ctx.structures[0];
    for (const unit of structure.units) {
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
                        start: StructureElement.Location.create(structure, unit, segs.offsets[start]),
                        end: StructureElement.Location.create(structure, unit, segs.offsets[prev]),
                        length: prev - start + 1,
                        element
                    };
                    break;
                }
            }
        }
    }

    return ssElements as SSElement<T>[];
}