/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ValueCell } from '../../mol-util';
import { TextureImage } from '../../mol-gl/renderable/util';
import { Vec2, Vec3, Vec4 } from '../../mol-math/linear-algebra';
import { LocationIterator } from '../util/location-iterator';
import { ColorTheme } from '../../mol-theme/color';
import { createNullTexture } from '../../mol-gl/webgl/texture';
import { ColorData, createColors, isColorTypeLocation } from './color-data';

export type InteriorColorData = {
    uInteriorThemeColor: ValueCell<Vec3>,
    uInteriorColorTexDim: ValueCell<Vec2>,
    tInteriorColor: ValueCell<TextureImage<Uint8Array>>,
    dInteriorColorType: ValueCell<string>,
}

const emptyInteriorColorTexture = { array: new Uint8Array(3), width: 1, height: 1 };

function createEmptyInteriorColors(interiorColorData?: InteriorColorData): InteriorColorData {
    if (interiorColorData) {
        ValueCell.update(interiorColorData.tInteriorColor, emptyInteriorColorTexture);
        ValueCell.update(interiorColorData.uInteriorColorTexDim, Vec2.set(interiorColorData.uInteriorColorTexDim.ref.value, 1, 1));
        ValueCell.updateIfChanged(interiorColorData.dInteriorColorType, 'none');
        return interiorColorData;
    }
    return {
        uInteriorThemeColor: ValueCell.create(Vec3()),
        uInteriorColorTexDim: ValueCell.create(Vec2.create(1, 1)),
        tInteriorColor: ValueCell.create(emptyInteriorColorTexture),
        dInteriorColorType: ValueCell.create('none'),
    };
}

/** Cells `createColors` writes but the interior color does not use, the palette is shared with the main color */
const unusedColorData = {
    tPalette: ValueCell.create(emptyInteriorColorTexture),
    dUsePalette: ValueCell.create(false),
    tColorGrid: ValueCell.create(createNullTexture()),
    uColorGridDim: ValueCell.create(Vec3()),
    uColorGridTransform: ValueCell.create(Vec4()),
    uPaletteDomain: ValueCell.create(Vec2()),
    uPaletteDefault: ValueCell.create(Vec3()),
};
const colorDataView = { ...unusedColorData } as ColorData;

function asColorData(data: InteriorColorData): ColorData {
    colorDataView.uColor = data.uInteriorThemeColor;
    colorDataView.uColorTexDim = data.uInteriorColorTexDim;
    colorDataView.tColor = data.tInteriorColor;
    colorDataView.dColorType = data.dInteriorColorType;
    return colorDataView;
}

export function createInteriorColors(locationIt: LocationIterator, positionIt: LocationIterator, colorTheme: ColorTheme<any, any>, themeColor: boolean, interiorColorData?: InteriorColorData): InteriorColorData {
    const { interiorColor } = colorTheme;
    const granularity = colorTheme.interiorGranularity ?? colorTheme.granularity;
    if (!themeColor || !interiorColor || !isColorTypeLocation(granularity)) {
        return createEmptyInteriorColors(interiorColorData);
    }
    const data = interiorColorData || createEmptyInteriorColors();
    createColors(locationIt, positionIt, { ...colorTheme, granularity, color: interiorColor, palette: undefined } as ColorTheme<any, any>, asColorData(data));
    return data;
}

/** Updates the interior colors of geometries that have them, a no-op for those that do not */
export function updateInteriorColors(values: object, locationIt: LocationIterator, positionIt: LocationIterator, colorTheme: ColorTheme<any, any>, themeColor: boolean) {
    if ('tInteriorColor' in values) {
        createInteriorColors(locationIt, positionIt, colorTheme, themeColor, values as InteriorColorData);
    }
}
