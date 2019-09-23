/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model } from '../../mol-model/structure/model';
import { Task } from '../../mol-task';
import { ModelFormat } from './format';
import { _parse_mmCif } from './mmcif/parser';
import { CifCategory, CifField } from '../../mol-io/reader/cif';
import { Column } from '../../mol-data/db';
import { mmCIF_Schema } from '../../mol-io/reader/cif/schema/mmcif';
import { EntityBuilder } from './common/entity';
import { File3DG } from '../../mol-io/reader/3dg/parser';
import { fillSerial } from '../../mol-util/array';
import { MoleculeType } from '../../mol-model/structure/model/types';

function getCategories(table: File3DG['table']) {
    const entityIds = new Array<string>(table._rowCount)
    const entityBuilder = new EntityBuilder()

    const seqIdStarts = table.position.toArray({ array: Uint32Array })
    const seqIdEnds = new Uint32Array(table._rowCount)
    const stride = seqIdStarts[1] - seqIdStarts[0]

    const objectRadius = stride / 3500

    for (let i = 0, il = table._rowCount; i < il; ++i) {
        const chr = table.chromosome.value(i)
        const entityId = entityBuilder.getEntityId(chr, MoleculeType.DNA, chr)
        entityIds[i] = entityId
        seqIdEnds[i] = seqIdStarts[i] + stride - 1
    }

    const ihm_sphere_obj_site: CifCategory.SomeFields<mmCIF_Schema['ihm_sphere_obj_site']> = {
        id: CifField.ofNumbers(fillSerial(new Uint32Array(table._rowCount))),
        entity_id: CifField.ofStrings(entityIds),
        seq_id_begin: CifField.ofNumbers(seqIdStarts),
        seq_id_end: CifField.ofNumbers(seqIdEnds),
        asym_id: CifField.ofColumn(table.chromosome),

        Cartn_x: CifField.ofNumbers(Column.mapToArray(table.x, x => x * 10, Float32Array)),
        Cartn_y: CifField.ofNumbers(Column.mapToArray(table.y, y => y * 10, Float32Array)),
        Cartn_z: CifField.ofNumbers(Column.mapToArray(table.z, z => z * 10, Float32Array)),

        object_radius: CifField.ofColumn(Column.ofConst(objectRadius, table._rowCount, Column.Schema.float)),
        rmsf: CifField.ofColumn(Column.ofConst(0, table._rowCount, Column.Schema.float)),
        model_id: CifField.ofColumn(Column.ofConst(1, table._rowCount, Column.Schema.int)),
    }

    return {
        entity: entityBuilder.getEntityCategory(),
        ihm_model_list: CifCategory.ofFields('ihm_model_list', {
            model_id: CifField.ofNumbers([1]),
            model_name: CifField.ofStrings(['3DG Model']),
        }),
        ihm_sphere_obj_site: CifCategory.ofFields('ihm_sphere_obj_site', ihm_sphere_obj_site)
    }
}

async function mmCifFrom3dg(file3dg: File3DG) {
    const categories = getCategories(file3dg.table)

    return {
        header: '3DG',
        categoryNames: Object.keys(categories),
        categories
    };
}

export function trajectoryFrom3DG(file3dg: File3DG): Task<Model.Trajectory> {
    return Task.create('Parse 3DG', async ctx => {
        await ctx.update('Converting to mmCIF');
        const cif = await mmCifFrom3dg(file3dg);
        const format = ModelFormat.mmCIF(cif);
        return _parse_mmCif(format, ctx);
    })
}
