/**
 * Copyright (c) 2021-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../mol-util/value-cell';
import { Vec2, Vec3, Vec4 } from '../../mol-math/linear-algebra';
import { TextureImage, createTextureImage } from '../../mol-gl/renderable/util';
import { createNullTexture, Texture } from '../../mol-gl/webgl/texture';
import { Material } from '../../mol-util/material';

export type SubstanceType = 'instance' | 'groupInstance' | 'volumeInstance';

export type SubstanceData = {
    tSubstance: ValueCell<TextureImage<Uint8Array>>
    uSubstanceTexDim: ValueCell<Vec2>
    dSubstance: ValueCell<boolean>,

    tSubstanceGrid: ValueCell<Texture>,
    uSubstanceGridDim: ValueCell<Vec3>,
    uSubstanceGridTransform: ValueCell<Vec4>,
    dSubstanceType: ValueCell<string>,
    uSubstanceStrength: ValueCell<number>,
}

export function applySubstanceMaterial(array: Uint8Array, start: number, end: number, material: Material) {
    for (let i = start; i < end; ++i) {
        Material.toArray(material, array, i * 4);
        array[i * 4 + 3] = 255;
    }
    return true;
}

export function clearSubstance(array: Uint8Array, start: number, end: number) {
    array.fill(0, start * 4, end * 4);
    return true;
}

export function createSubstance(count: number, type: SubstanceType, substanceData?: SubstanceData): SubstanceData {
    const substance = createTextureImage(Math.max(1, count), 4, Uint8Array, substanceData && substanceData.tSubstance.ref.value.array);
    if (substanceData) {
        ValueCell.update(substanceData.tSubstance, substance);
        ValueCell.update(substanceData.uSubstanceTexDim, Vec2.create(substance.width, substance.height));
        ValueCell.updateIfChanged(substanceData.dSubstance, count > 0);
        ValueCell.updateIfChanged(substanceData.dSubstanceType, type);
        return substanceData;
    } else {
        return {
            tSubstance: ValueCell.create(substance),
            uSubstanceTexDim: ValueCell.create(Vec2.create(substance.width, substance.height)),
            dSubstance: ValueCell.create(count > 0),

            tSubstanceGrid: ValueCell.create(createNullTexture()),
            uSubstanceGridDim: ValueCell.create(Vec3.create(1, 1, 1)),
            uSubstanceGridTransform: ValueCell.create(Vec4.create(0, 0, 0, 1)),
            dSubstanceType: ValueCell.create(type),
            uSubstanceStrength: ValueCell.create(1),
        };
    }
}

const emptySubstanceTexture = { array: new Uint8Array(4), width: 1, height: 1 };
export function createEmptySubstance(substanceData?: SubstanceData): SubstanceData {
    if (substanceData) {
        ValueCell.update(substanceData.tSubstance, emptySubstanceTexture);
        ValueCell.update(substanceData.uSubstanceTexDim, Vec2.create(1, 1));
        return substanceData;
    } else {
        return {
            tSubstance: ValueCell.create(emptySubstanceTexture),
            uSubstanceTexDim: ValueCell.create(Vec2.create(1, 1)),
            dSubstance: ValueCell.create(false),

            tSubstanceGrid: ValueCell.create(createNullTexture()),
            uSubstanceGridDim: ValueCell.create(Vec3.create(1, 1, 1)),
            uSubstanceGridTransform: ValueCell.create(Vec4.create(0, 0, 0, 1)),
            dSubstanceType: ValueCell.create('groupInstance'),
            uSubstanceStrength: ValueCell.create(1),
        };
    }
}