/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CoarseHierarchy, CoarseConformation, CoarseElementData, CoarseSphereConformation, CoarseGaussianConformation } from '../../../mol-model/structure/model/properties/coarse';
import { Entities } from '../../../mol-model/structure/model/properties/common';
import { Column } from '../../../mol-data/db';
import { getCoarseKeys } from '../../../mol-model/structure/model/properties/utils/coarse-keys';
import { UUID } from '../../../mol-util';
import { Segmentation, Interval } from '../../../mol-data/int';
import { Mat3, Tensor } from '../../../mol-math/linear-algebra';
import { ElementIndex, ChainIndex } from '../../../mol-model/structure/model/indexing';
import { getCoarseRanges } from '../../../mol-model/structure/model/properties/utils/coarse-ranges';
import { IhmSphereObjSite, IhmGaussianObjSite, AtomSite, BasicSchema } from './schema';
import { Model } from '../../../mol-model/structure';

export interface CoarseData {
    model_id: number,
    model_name: string,
    model_group_name: string,
    entities: Entities,
    atom_site: AtomSite,
    atom_site_sourceIndex: Column<number>,
    ihm_sphere_obj_site: IhmSphereObjSite,
    ihm_gaussian_obj_site: IhmGaussianObjSite
}

export const EmptyCoarse = { hierarchy: CoarseHierarchy.Empty, conformation: void 0 as any };

export function getCoarse(data: CoarseData, properties: Model['properties']): { hierarchy: CoarseHierarchy, conformation: CoarseConformation } {
    const { ihm_sphere_obj_site, ihm_gaussian_obj_site } = data;

    if (ihm_sphere_obj_site._rowCount === 0 && ihm_gaussian_obj_site._rowCount === 0) return EmptyCoarse;

    const sphereData = getData(ihm_sphere_obj_site);
    const sphereConformation = getSphereConformation(ihm_sphere_obj_site);
    const sphereKeys = getCoarseKeys(sphereData, data.entities);
    const sphereRanges = getCoarseRanges(sphereData, properties.chemicalComponentMap);

    const gaussianData = getData(ihm_gaussian_obj_site);
    const gaussianConformation = getGaussianConformation(ihm_gaussian_obj_site);
    const gaussianKeys = getCoarseKeys(gaussianData, data.entities);
    const gaussianRanges = getCoarseRanges(gaussianData, properties.chemicalComponentMap);

    return {
        hierarchy: {
            isDefined: true,
            spheres: { ...sphereData, ...sphereKeys, ...sphereRanges },
            gaussians: { ...gaussianData, ...gaussianKeys, ...gaussianRanges },
        },
        conformation: {
            id: UUID.create22(),
            spheres: sphereConformation,
            gaussians: gaussianConformation
        }
    };
}

function getSphereConformation(data: IhmSphereObjSite): CoarseSphereConformation {
    return {
        x: data.Cartn_x.toArray({ array: Float32Array }),
        y: data.Cartn_y.toArray({ array: Float32Array }),
        z: data.Cartn_z.toArray({ array: Float32Array }),
        radius: data.object_radius.toArray({ array: Float32Array }),
        rmsf: data.rmsf.toArray({ array: Float32Array })
    };
}

function getGaussianConformation(data: IhmGaussianObjSite): CoarseGaussianConformation {
    const matrix_space = BasicSchema.ihm_gaussian_obj_site.covariance_matrix.space;
    const covariance_matrix: Mat3[] = [];
    const { covariance_matrix: cm } = data;

    for (let i = 0, _i = cm.rowCount; i < _i; i++) {
        covariance_matrix[i] = Tensor.toMat3(Mat3(), matrix_space, cm.value(i));
    }

    return {
        x: data.mean_Cartn_x.toArray({ array: Float32Array }),
        y: data.mean_Cartn_y.toArray({ array: Float32Array }),
        z: data.mean_Cartn_z.toArray({ array: Float32Array }),
        weight: data.weight.toArray({ array: Float32Array }),
        covariance_matrix
    };
}

function getSegments(asym_id: Column<string>, seq_id_begin: Column<number>, seq_id_end: Column<number>) {
    const chainOffsets = [0 as ElementIndex];
    for (let i = 1, _i = asym_id.rowCount; i < _i; i++) {
        const newChain = !asym_id.areValuesEqual(i - 1, i);
        if (newChain) chainOffsets[chainOffsets.length] = i as ElementIndex;
    }

    return {
        chainElementSegments: Segmentation.ofOffsets<ElementIndex, ChainIndex>(chainOffsets, Interval.ofBounds(0, asym_id.rowCount))
    };
}

function getData(data: IhmSphereObjSite | IhmGaussianObjSite): CoarseElementData {
    const { entity_id, seq_id_begin, seq_id_end, asym_id } = data;
    return { count: entity_id.rowCount, entity_id, asym_id, seq_id_begin, seq_id_end, ...getSegments(asym_id, seq_id_begin, seq_id_end) };
}