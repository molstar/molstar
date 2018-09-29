/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mmCIF_Database, mmCIF_Schema } from 'mol-io/reader/cif/schema/mmcif';
import { unionMany } from 'mol-util/set';
import { Model } from '../../model';
import { Structure } from '../../structure';
import { EntityIndex } from '../../model/indexing';
import { UniqueArray } from 'mol-data/generic';
import { sortArray } from 'mol-data/util';

export function getModelMmCifCategory<K extends keyof mmCIF_Schema>(model: Model, name: K): mmCIF_Database[K] | undefined {
    if (model.sourceData.kind !== 'mmCIF') return;
    return model.sourceData.data[name];
}

export function getUniqueResidueNamesFromStructures(structures: Structure[]) {
    return unionMany(structures.map(s => s.uniqueResidueNames));
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