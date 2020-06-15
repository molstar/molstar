/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '../../../mol-math/linear-algebra';
import { AtomicConformation } from './properties/atomic';
import { CoarseConformation } from './properties/coarse';
import { arrayMinMax } from '../../../mol-util/array';
import { Model } from './model';

export function calcModelCenter(atomicConformation: AtomicConformation, coarseConformation?: CoarseConformation) {
    let rangesX: number[] = [];
    let rangesY: number[] = [];
    let rangesZ: number[] = [];

    if (atomicConformation.x.length) {
        rangesX.push(...arrayMinMax(atomicConformation.x));
        rangesY.push(...arrayMinMax(atomicConformation.y));
        rangesZ.push(...arrayMinMax(atomicConformation.z));
    }

    if (coarseConformation) {
        if (coarseConformation.spheres.x.length) {
            rangesX.push(...arrayMinMax(coarseConformation.spheres.x));
            rangesY.push(...arrayMinMax(coarseConformation.spheres.y));
            rangesZ.push(...arrayMinMax(coarseConformation.spheres.z));
        }
        if (coarseConformation.gaussians.x.length) {
            rangesX.push(...arrayMinMax(coarseConformation.gaussians.x));
            rangesY.push(...arrayMinMax(coarseConformation.gaussians.y));
            rangesZ.push(...arrayMinMax(coarseConformation.gaussians.z));
        }
    }

    const [minX, maxX] = arrayMinMax(rangesX);
    const [minY, maxY] = arrayMinMax(rangesY);
    const [minZ, maxZ] = arrayMinMax(rangesZ);

    const x = minX + (maxX - minX) / 2;
    const y = minY + (maxY - minY) / 2;
    const z = minZ + (maxZ - minZ) / 2;

    return Vec3.create(x, y, z);
}

export function getAsymIdCount(model: Model) {
    const auth = new Set<string>();
    const label = new Set<string>();
    model.properties.structAsymMap.forEach(({ auth_id }, label_id) => {
        auth.add(auth_id);
        label.add(label_id);
    });
    return { auth: auth.size, label: label.size };
}