/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model } from '../../mol-model/structure/model';
import { Task } from '../../mol-task';
import { ModelFormat } from '../format';
import { Column, Table } from '../../mol-data/db';
import { EntityBuilder } from './common/entity';
import { File3DG } from '../../mol-io/reader/3dg/parser';
import { fillSerial } from '../../mol-util/array';
import { MoleculeType } from '../../mol-model/structure/model/types';
import { BasicSchema, createBasic } from './basic/schema';
import { createModels } from './basic/parser';

function getBasic(table: File3DG['table']) {
    const entityIds = new Array<string>(table._rowCount);
    const entityBuilder = new EntityBuilder();

    const seqIdStarts = table.position.toArray({ array: Uint32Array });
    const seqIdEnds = new Uint32Array(table._rowCount);
    const stride = seqIdStarts[1] - seqIdStarts[0];

    const objectRadius = stride / 3500;

    for (let i = 0, il = table._rowCount; i < il; ++i) {
        const chr = table.chromosome.value(i);
        const entityId = entityBuilder.getEntityId(chr, MoleculeType.DNA, chr);
        entityIds[i] = entityId;
        seqIdEnds[i] = seqIdStarts[i] + stride - 1;
    }

    const ihm_sphere_obj_site = Table.ofPartialColumns(BasicSchema.ihm_sphere_obj_site, {
        id: Column.ofIntArray(fillSerial(new Uint32Array(table._rowCount))),
        entity_id: Column.ofStringArray(entityIds),
        seq_id_begin: Column.ofIntArray(seqIdStarts),
        seq_id_end: Column.ofIntArray(seqIdEnds),
        asym_id: table.chromosome,

        Cartn_x: Column.ofFloatArray(Column.mapToArray(table.x, x => x * 10, Float32Array)),
        Cartn_y: Column.ofFloatArray(Column.mapToArray(table.y, y => y * 10, Float32Array)),
        Cartn_z: Column.ofFloatArray(Column.mapToArray(table.z, z => z * 10, Float32Array)),

        object_radius: Column.ofConst(objectRadius, table._rowCount, Column.Schema.float),
        rmsf: Column.ofConst(0, table._rowCount, Column.Schema.float),
        model_id: Column.ofConst(1, table._rowCount, Column.Schema.int),
    }, table._rowCount);

    return createBasic({
        entity: entityBuilder.getEntityTable(),
        ihm_model_list: Table.ofPartialColumns(BasicSchema.ihm_model_list, {
            model_id: Column.ofIntArray([1]),
            model_name: Column.ofStringArray(['3DG Model']),
        }, 1),
        ihm_sphere_obj_site
    });
}

//

export { Format3dg };

type Format3dg = ModelFormat<File3DG>

namespace Format3dg {
    export function is(x: ModelFormat): x is Format3dg {
        return x.kind === '3dg';
    }

    export function from3dg(file3dg: File3DG): Format3dg {
        return { kind: '3dg', name: '3DG', data: file3dg };
    }
}

export function trajectoryFrom3DG(file3dg: File3DG): Task<Model.Trajectory> {
    return Task.create('Parse 3DG', async ctx => {
        const format = Format3dg.from3dg(file3dg);
        const basic = getBasic(file3dg.table);
        return createModels(basic, format, ctx);
    });
}
