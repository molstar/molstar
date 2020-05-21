/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column, Table } from '../../../mol-data/db';
import { RuntimeContext } from '../../../mol-task';
import UUID from '../../../mol-util/uuid';
import { Model } from '../../../mol-model/structure/model/model';
import { Entities } from '../../../mol-model/structure/model/properties/common';
import { CustomProperties } from '../../../mol-model/custom-property';
import { getAtomicHierarchyAndConformation } from './atomic';
import { getCoarse, EmptyCoarse, CoarseData } from './coarse';
import { getSequence } from './sequence';
import { sortAtomSite } from './sort';
import { ModelFormat } from '../../format';
import { getAtomicRanges } from '../../../mol-model/structure/model/properties/utils/atomic-ranges';
import { AtomSite, BasicData } from './schema';
import { getProperties } from './properties';
import { getEntities } from './entities';
import { getModelGroupName } from './util';

export async function createModels(data: BasicData, format: ModelFormat, ctx: RuntimeContext) {
    const properties = getProperties(data);
    const models = data.ihm_model_list._rowCount > 0
        ? await readIntegrative(ctx, data, properties, format)
        : await readStandard(ctx, data, properties, format);

    for (let i = 0; i < models.length; i++) {
        Model.TrajectoryInfo.set(models[i], { index: i, size: models.length });
    }

    return models;
}

/** Standard atomic model */
function createStandardModel(data: BasicData, atom_site: AtomSite, sourceIndex: Column<number>, entities: Entities, properties: Model['properties'], format: ModelFormat, previous?: Model): Model {

    const atomic = getAtomicHierarchyAndConformation(atom_site, sourceIndex, entities, properties.chemicalComponentMap, format, previous);
    const modelNum = atom_site.pdbx_PDB_model_num.value(0);
    if (previous && atomic.sameAsPrevious) {
        return {
            ...previous,
            id: UUID.create22(),
            modelNum,
            atomicConformation: atomic.conformation,
            _dynamicPropertyData: Object.create(null)
        };
    }

    const coarse = EmptyCoarse;
    const sequence = getSequence(data, entities, atomic.hierarchy, coarse.hierarchy);
    const atomicRanges = getAtomicRanges(atomic.hierarchy, entities, atomic.conformation, sequence);

    const entry = data.entry.id.valueKind(0) === Column.ValueKind.Present
        ? data.entry.id.value(0) : format.name;

    const label: string[] = [];
    if (entry) label.push(entry);
    if (data.struct.title.valueKind(0) === Column.ValueKind.Present) label.push(data.struct.title.value(0));

    return {
        id: UUID.create22(),
        entryId: entry,
        label: label.join(' | '),
        entry,
        sourceData: format,
        modelNum,
        entities,
        sequence,
        atomicHierarchy: atomic.hierarchy,
        atomicConformation: atomic.conformation,
        atomicRanges,
        atomicChainOperatorMappinng: atomic.chainOperatorMapping,
        coarseHierarchy: coarse.hierarchy,
        coarseConformation: coarse.conformation,
        properties,
        customProperties: new CustomProperties(),
        _staticPropertyData: Object.create(null),
        _dynamicPropertyData: Object.create(null)
    };
}

/** Integrative model with atomic/coarse parts */
function createIntegrativeModel(data: BasicData, ihm: CoarseData, properties: Model['properties'], format: ModelFormat): Model {
    const atomic = getAtomicHierarchyAndConformation(ihm.atom_site, ihm.atom_site_sourceIndex, ihm.entities, properties.chemicalComponentMap, format);
    const coarse = getCoarse(ihm, properties);
    const sequence = getSequence(data, ihm.entities, atomic.hierarchy, coarse.hierarchy);
    const atomicRanges = getAtomicRanges(atomic.hierarchy, ihm.entities, atomic.conformation, sequence);

    const entry = data.entry.id.valueKind(0) === Column.ValueKind.Present
        ? data.entry.id.value(0) : format.name;

    const label: string[] = [];
    if (entry) label.push(entry);
    if (data.struct.title.valueKind(0) === Column.ValueKind.Present) label.push(data.struct.title.value(0));
    if (ihm.model_name) label.push(ihm.model_name);
    if (ihm.model_group_name) label.push(ihm.model_group_name);

    return {
        id: UUID.create22(),
        entryId: entry,
        label: label.join(' | '),
        entry,
        sourceData: format,
        modelNum: ihm.model_id,
        entities: ihm.entities,
        sequence,
        atomicHierarchy: atomic.hierarchy,
        atomicConformation: atomic.conformation,
        atomicRanges,
        atomicChainOperatorMappinng: atomic.chainOperatorMapping,
        coarseHierarchy: coarse.hierarchy,
        coarseConformation: coarse.conformation,
        properties,
        customProperties: new CustomProperties(),
        _staticPropertyData: Object.create(null),
        _dynamicPropertyData: Object.create(null)
    };
}

function findModelEnd(num: Column<number>, startIndex: number) {
    const rowCount = num.rowCount;
    if (!num.isDefined) return rowCount;
    let endIndex = startIndex + 1;
    while (endIndex < rowCount && num.areValuesEqual(startIndex, endIndex)) endIndex++;
    return endIndex;
}

async function readStandard(ctx: RuntimeContext, data: BasicData, properties: Model['properties'], format: ModelFormat) {
    const models: Model[] = [];

    if (data.atom_site) {
        const atomCount = data.atom_site.id.rowCount;
        const entities = getEntities(data, properties);

        let modelStart = 0;
        while (modelStart < atomCount) {
            const modelEnd = findModelEnd(data.atom_site.pdbx_PDB_model_num, modelStart);
            const { atom_site, sourceIndex } = await sortAtomSite(ctx, data.atom_site, modelStart, modelEnd);
            const model = createStandardModel(data, atom_site, sourceIndex, entities, properties, format, models.length > 0 ? models[models.length - 1] : void 0);
            models.push(model);
            modelStart = modelEnd;
        }
    }
    return models;
}

function splitTable<T extends Table<any>>(table: T, col: Column<number>) {
    const ret = new Map<number, { table: T, start: number, end: number }>();
    const rowCount = table._rowCount;
    let modelStart = 0;
    while (modelStart < rowCount) {
        const modelEnd = findModelEnd(col, modelStart);
        const id = col.value(modelStart);
        ret.set(id, {
            table: Table.window(table, table._schema, modelStart, modelEnd) as T,
            start: modelStart,
            end: modelEnd
        });
        modelStart = modelEnd;
    }
    return ret;
}



async function readIntegrative(ctx: RuntimeContext, data: BasicData, properties: Model['properties'], format: ModelFormat) {
    const entities = getEntities(data, properties);
    // when `atom_site.ihm_model_id` is undefined fall back to `atom_site.pdbx_PDB_model_num`
    const atom_sites_modelColumn = data.atom_site.ihm_model_id.isDefined
        ? data.atom_site.ihm_model_id : data.atom_site.pdbx_PDB_model_num;
    const atom_sites = splitTable(data.atom_site, atom_sites_modelColumn);

    // TODO: will coarse IHM records require sorting or will we trust it?
    // ==> Probably implement a sort as as well and store the sourceIndex same as with atomSite
    // If the sorting is implemented, updated mol-model/structure/properties: atom.sourceIndex
    const sphere_sites = splitTable(data.ihm_sphere_obj_site, data.ihm_sphere_obj_site.model_id);
    const gauss_sites = splitTable(data.ihm_gaussian_obj_site, data.ihm_gaussian_obj_site.model_id);

    const models: Model[] = [];

    if (data.ihm_model_list) {
        const { model_id, model_name } = data.ihm_model_list;
        for (let i = 0; i < data.ihm_model_list._rowCount; i++) {
            const id = model_id.value(i);

            let atom_site, atom_site_sourceIndex;
            if (atom_sites.has(id)) {
                const e = atom_sites.get(id)!;
                // need to sort `data.atom_site` as `e.start` and `e.end` are indices into that
                const { atom_site: sorted, sourceIndex } = await sortAtomSite(ctx, data.atom_site, e.start, e.end);
                atom_site = sorted;
                atom_site_sourceIndex = sourceIndex;
            } else {
                atom_site = Table.window(data.atom_site, data.atom_site._schema, 0, 0);
                atom_site_sourceIndex = Column.ofIntArray([]);
            }

            const ihm: CoarseData = {
                model_id: id,
                model_name: model_name.value(i),
                model_group_name: getModelGroupName(id, data),
                entities: entities,
                atom_site,
                atom_site_sourceIndex,
                ihm_sphere_obj_site: sphere_sites.has(id) ? sphere_sites.get(id)!.table : Table.window(data.ihm_sphere_obj_site, data.ihm_sphere_obj_site._schema, 0, 0),
                ihm_gaussian_obj_site: gauss_sites.has(id) ? gauss_sites.get(id)!.table : Table.window(data.ihm_gaussian_obj_site, data.ihm_gaussian_obj_site._schema, 0, 0)
            };
            const model = createIntegrativeModel(data, ihm, properties, format);
            models.push(model);
        }
    }

    return models;
}