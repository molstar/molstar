/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { mmCIF_Database, mmCIF_Schema } from '@molstar/io/reader/cif/schema/mmcif';
import { SetUtils } from '@molstar/core/util/set';
import type { Model } from '../../model.js';
import type { Structure } from '../../structure.js';
import type { EntityIndex } from '../../model/indexing.js';
import { UniqueArray } from '@molstar/core/data/generic';
import { sortArray } from '@molstar/core/data/util';
import { CifWriter } from '@molstar/io/writer/cif';
import type { CifExportContext } from '../mmcif.js';
import { MmcifFormat } from '@molstar/model/formats/structure/mmcif';
import { type CifCategory, type CifField, getCifFieldType } from '@molstar/io/reader/cif';

export function getModelMmCifCategory<K extends keyof mmCIF_Schema>(model: Model, name: K): mmCIF_Database[K] | undefined {
    if (!MmcifFormat.is(model.sourceData)) return;
    return model.sourceData.data.db[name];
}

export function getUniqueResidueNamesFromStructures(structures: Structure[]) {
    return SetUtils.unionMany(...structures.map(s => s.uniqueResidueNames));
}

export function getUniqueEntityIdsFromStructures(structures: Structure[]): Set<string> {
    if (structures.length === 0) return new Set();

    const names = structures[0].model.entities.data.id;
    return new Set(getUniqueEntityIndicesFromStructures(structures).map(i => names.value(i)));
}

export function getUniqueEntityIndicesFromStructures(structures: Structure[]): ReadonlyArray<EntityIndex> {
    if (structures.length === 0) return [];
    if (structures.length === 1) return structures[0].entityIndices;
    const ret = UniqueArray.create<EntityIndex, EntityIndex>();
    for (const s of structures) {
        for (const e of s.entityIndices) {
            UniqueArray.add(ret, e, e);
        }
    }
    sortArray(ret.array);
    return ret.array;
}

export function copy_mmCif_category(name: keyof mmCIF_Schema, condition?: (structure: Structure) => boolean): CifWriter.Category<CifExportContext> {
    return {
        name,
        instance({ structures }) {
            if (condition && !condition(structures[0])) return CifWriter.Category.Empty;

            const model = structures[0].model;
            if (!MmcifFormat.is(model.sourceData)) return CifWriter.Category.Empty;

            const table = model.sourceData.data.db[name];
            if (!table || !table._rowCount) return CifWriter.Category.Empty;
            return CifWriter.Category.ofTable(table);
        }
    };
}

export function copy_source_mmCifCategory(encoder: CifWriter.Encoder, ctx: CifExportContext, category: CifCategory): CifWriter.Category<CifExportContext> | undefined {
    if (!MmcifFormat.is(ctx.firstModel.sourceData)) return;

    const fs = CifWriter.fields<number, undefined>();
    if (encoder.isBinary) {
        for (const f of category.fieldNames) {
            // TODO: this could be optimized
            const field = classifyField(f, category.getField(f)!);
            fs.add(field);
        }
    } else {
        for (const f of category.fieldNames) {
            const field = category.getField(f)!;
            fs.str(f, row => field.str(row));
        }
    }

    const fields = fs.getFields();
    return {
        name: category.name,
        instance() {
            return { fields, source: [{ data: void 0, rowCount: category.rowCount }] };
        }
    };
}

function classifyField(name: string, field: CifField): CifWriter.Field {
    const type = getCifFieldType(field);
    if (type['@type'] === 'str') {
        return { name, type: CifWriter.Field.Type.Str, value: field.str, valueKind: field.valueKind };
    } else if (type['@type'] === 'float') {
        return CifWriter.Field.float(name, field.float, { valueKind: field.valueKind, typedArray: Float64Array });
    } else {
        return CifWriter.Field.int(name, field.int, { valueKind: field.valueKind, typedArray: Int32Array });
    }
}