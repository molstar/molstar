/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mmCIF_Database as mmCIF, mmCIF_Schema } from 'mol-io/reader/cif/schema/mmcif'
import CoarseGrained from '../../properties/coarse-grained'
import { Entities } from '../../properties/common';
import { Column } from 'mol-data/db';

function coarseGrainedFromIHM(data: mmCIF, entities: Entities): CoarseGrained {
    if (data.ihm_model_list._rowCount === 0) return CoarseGrained.Empty;

    const { ihm_model_list, ihm_sphere_obj_site, ihm_gaussian_obj_site } = data;
    const modelIndex = Column.createIndexer(ihm_model_list.model_id);

    return {
        isDefined: true,
        modelList: ihm_model_list,
        spheres: getSpheres(ihm_sphere_obj_site, entities, modelIndex),
        gaussians: getGaussians(ihm_gaussian_obj_site, entities, modelIndex)
    };
}

function getSpheres(data: mmCIF['ihm_sphere_obj_site'], entities: Entities, modelIndex: (id: number) => number): CoarseGrained.Spheres {
    const { Cartn_x, Cartn_y, Cartn_z, object_radius: radius, rmsf } = data;
    const x = Cartn_x.toArray({ array: Float32Array });
    const y = Cartn_y.toArray({ array: Float32Array });
    const z = Cartn_z.toArray({ array: Float32Array });
    return { count: x.length, ...getCommonColumns(data, entities, modelIndex), x, y, z, radius, rmsf };
}

function getGaussians(data: mmCIF['ihm_gaussian_obj_site'], entities: Entities, modelIndex: (id: number) => number): CoarseGrained.Gaussians {
    const { mean_Cartn_x, mean_Cartn_y, mean_Cartn_z, weight, covariance_matrix  } = data;
    const x = mean_Cartn_x.toArray({ array: Float32Array });
    const y = mean_Cartn_y.toArray({ array: Float32Array });
    const z = mean_Cartn_z.toArray({ array: Float32Array });
    return { count: x.length, ...getCommonColumns(data, entities, modelIndex), x, y, z, weight, covariance_matrix, matrix_space: mmCIF_Schema.ihm_gaussian_obj_site.covariance_matrix.space };
}

function getCommonColumns(data: mmCIF['ihm_sphere_obj_site'] | mmCIF['ihm_gaussian_obj_site'], entities: Entities, modelIndex: (id: number) => number) {
    const { model_id, entity_id, seq_id_begin, seq_id_end, asym_id } = data;

    return {
        entityKey: Column.mapToArray(entity_id, id => entities.getEntityIndex(id), Int32Array),
        modelKey: Column.mapToArray(model_id, modelIndex, Int32Array),
        asym_id,
        seq_id_begin,
        seq_id_end
    };
}

export { coarseGrainedFromIHM }