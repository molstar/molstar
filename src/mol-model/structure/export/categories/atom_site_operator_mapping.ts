/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SymmetryOperator } from '../../../../mol-math/geometry';
import { CifExportContext } from '../mmcif';
import { StructureElement, StructureProperties as P, CifExportCategoryInfo } from '../../structure';
import Unit from '../../structure/unit';
import { Segmentation } from '../../../../mol-data/int';
import { CifWriter } from '../../../../mol-io/writer/cif';
import { Column } from '../../../../mol-data/db';


export function atom_site_operator_mapping(ctx: CifExportContext): CifExportCategoryInfo | undefined {
    const entries = getEntries(ctx);
    if (entries.length === 0) return;
    return [Category, entries, { ignoreFilter: true }];
}

export const AtomSiteOperatorMappingCategoryName = 'molstar_atom_site_operator_mapping';

export const AtomSiteOperatorMappingSchema = {
    molstar_atom_site_operator_mapping: {
        label_asym_id: Column.Schema.Str(),
        auth_asym_id: Column.Schema.Str(),
        operator_name: Column.Schema.Str(),
        suffix: Column.Schema.Str(),

        // assembly
        assembly_id: Column.Schema.Str(),
        assembly_operator_id: Column.Schema.Int(),

        // symmetry
        symmetry_operator_index: Column.Schema.Int(),
        symmetry_hkl: Column.Schema.Vector(3),

        // NCS
        ncs_id: Column.Schema.Str(),
    }
};

const asmValueKind = (i: number, xs: Entry[]) => typeof xs[i].operator.assembly === 'undefined' ? Column.ValueKind.NotPresent : Column.ValueKind.Present;
const symmetryValueKind = (i: number, xs: Entry[]) => xs[i].operator.spgrOp === -1 ? Column.ValueKind.NotPresent : Column.ValueKind.Present;

const Fields = CifWriter.fields<number, Entry[], keyof (typeof AtomSiteOperatorMappingSchema)['molstar_atom_site_operator_mapping']>()
    .str('label_asym_id', (i, xs) => xs[i].label_asym_id)
    .str('auth_asym_id', (i, xs) => xs[i].auth_asym_id)
    .str('operator_name', (i, xs) => xs[i].operator.name)
    .str('suffix', (i, xs) => xs[i].operator.suffix)
    // assembly
    // TODO: include oper list as well?
    .str('assembly_id', (i, xs) => xs[i].operator.assembly?.id || '', { valueKind: asmValueKind })
    .int('assembly_operator_id', (i, xs) => xs[i].operator.assembly?.operId || 0, { valueKind: asmValueKind })
    // symmetry
    .int('symmetry_operator_index', (i, xs) => xs[i].operator.spgrOp, { valueKind: symmetryValueKind })
    .vec('symmetry_hkl', [(i, xs) => xs[i].operator.hkl[0], (i, xs) => xs[i].operator.hkl[1], (i, xs) => xs[i].operator.hkl[2]], { valueKind: symmetryValueKind })
    // NCS
    .str('ncs_id', (i, xs) => xs[i].operator.ncsId || '', { valueKind: (i, xs) => !xs[i].operator.ncsId ? Column.ValueKind.NotPresent : Column.ValueKind.Present })
    .getFields();

const Category: CifWriter.Category<Entry[]> = {
    name: 'molstar_atom_site_operator_mapping',
    instance(entries: Entry[]) {
        return { fields: Fields, source: [{ data: entries, rowCount: entries.length }] };
    }
};

interface Entry {
    label_asym_id: string,
    auth_asym_id: string,
    operator: SymmetryOperator
}

function getEntries(ctx: CifExportContext) {
    const existing = new Set<string>();
    const entries: Entry[] = [];

    for (const s of ctx.structures) {
        const l = StructureElement.Location.create(s);
        for (const unit of s.units) {
            const operator = unit.conformation.operator;
            if (!operator.suffix || unit.kind !== Unit.Kind.Atomic) continue;

            l.unit = unit;

            const { elements } = unit;
            const chainsIt = Segmentation.transientSegments(unit.model.atomicHierarchy.chainAtomSegments, elements);
            while (chainsIt.hasNext) {
                const chainSegment = chainsIt.move();
                l.element = elements[chainSegment.start];

                const label_asym_id = P.chain.label_asym_id(l);
                const key = `${label_asym_id}${operator.suffix}`;

                if (existing.has(key)) continue;
                existing.add(key);

                const auth_asym_id = P.chain.label_asym_id(l);
                entries.push({ label_asym_id, auth_asym_id, operator });
            }
        }
    }

    return entries;
}
