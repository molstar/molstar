/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from '../../mol-task';
import { Column, Table } from '../../mol-data/db';
import { MoleculeType } from '../../mol-model/structure/model/types';
import { EntityBuilder } from '../../mol-model-formats/structure/common/entity';
import { BasicSchema, createBasic } from '../../mol-model-formats/structure/basic/schema';
import { createModels } from '../../mol-model-formats/structure/basic/parser';
import { G3dDataBlock } from './data';
import { objectForEach } from '../../mol-util/object';
import { Trajectory } from '../../mol-model/structure';

interface Columns {
    entity_id: string[],
    chromosome: string[],
    start: Int32Array,
    end: Int32Array,
    x: Float32Array
    y: Float32Array
    z: Float32Array
    haplotype: string[]
}

function getColumns(block: G3dDataBlock) {
    const { data } = block;
    let size = 0;

    objectForEach(data, h => objectForEach(h, g => size += g.start.length));

    const columns: Columns = {
        entity_id: new Array(size),
        chromosome: new Array(size),
        start: new Int32Array(size),
        end: new Int32Array(size),
        x: new Float32Array(size),
        y: new Float32Array(size),
        z: new Float32Array(size),
        haplotype: new Array(size)
    };

    let o = 0;
    objectForEach(data, (hs, h) => {
        objectForEach(hs, (chs, ch) => {
            const entity_id = `${ch}-${h}`;
            for (let i = 0, _i = chs.start.length; i < _i; i++) {
                columns.entity_id[o] = entity_id;
                columns.chromosome[o] = ch;
                columns.start[o] = o + 1;
                columns.end[o] = o + 2;
                columns.x[o] = 10 * chs.x[i];
                columns.y[o] = 10 * chs.y[i];
                columns.z[o] = 10 * chs.z[i];
                columns.haplotype[o] = h;
                o++;
            }
        });
    });

    return columns;
}

function getBasic(data: G3dDataBlock) {
    const columns = getColumns(data);

    const rowCount = columns.start.length;
    const entityIds = new Array<string>(rowCount);
    const entityBuilder = new EntityBuilder();

    const stride = columns.start[1] - columns.start[0];

    const objectRadius = stride / 3500;

    for (let i = 0; i < rowCount; ++i) {
        const e = columns.entity_id[i];
        const entityId = entityBuilder.getEntityId(e, MoleculeType.DNA, e);
        entityIds[i] = entityId;
    }

    const ihm_sphere_obj_site = Table.ofPartialColumns(BasicSchema.ihm_sphere_obj_site, {
        id: Column.range(0, rowCount),
        entity_id: Column.ofStringArray(entityIds),
        seq_id_begin: Column.ofIntArray(columns.start),
        seq_id_end: Column.ofIntArray(columns.end),
        asym_id: Column.ofStringArray(columns.chromosome),

        Cartn_x: Column.ofFloatArray(columns.x),
        Cartn_y: Column.ofFloatArray(columns.y),
        Cartn_z: Column.ofFloatArray(columns.z),

        object_radius: Column.ofConst(objectRadius, rowCount, Column.Schema.float),
        rmsf: Column.ofConst(0, rowCount, Column.Schema.float),
        model_id: Column.ofConst(1, rowCount, Column.Schema.int),
    }, rowCount);

    return createBasic({
        entity: entityBuilder.getEntityTable(),
        ihm_model_list: Table.ofPartialColumns(BasicSchema.ihm_model_list, {
            model_id: Column.ofIntArray([1]),
            model_name: Column.ofStringArray(['3DG Model']),
        }, 1),
        ihm_sphere_obj_site
    });
}

export function trajectoryFromG3D(data: G3dDataBlock): Task<Trajectory> {
    return Task.create('Parse G3D', async ctx => {
        const basic = getBasic(data);
        return createModels(basic, { kind: 'g3d', name: 'G3D', data }, ctx);
    });
}