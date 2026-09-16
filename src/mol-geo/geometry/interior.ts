/**
 * Copyright (c) 2025-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { Vec4 } from '../../mol-math/linear-algebra/3d/vec4';
import { Color } from '../../mol-util/color/color';
import { Material } from '../../mol-util/material';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ValueCell } from '../../mol-util/value-cell';

export type InteriorData = {
    uInteriorColor: ValueCell<Vec4>,
    uInteriorSubstance: ValueCell<Vec4>,
}

export function getInteriorParam() {
    return PD.Group({
        themeColor: PD.Boolean(false, { description: 'Use the interior color of the color theme, when it provides one' }),
        color: PD.Color(Color.fromRgb(76, 76, 76)),
        colorStrength: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }),
        substance: Material.getParam(),
        substanceStrength: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
    });
}
export type InteriorParam = ReturnType<typeof getInteriorParam>
export type InteriorProps = InteriorParam['defaultValue'];

export function areInteriorPropsEquals(a: InteriorProps, b: InteriorProps): boolean {
    return a.themeColor === b.themeColor
        && a.color === b.color
        && a.colorStrength === b.colorStrength
        && Material.areEqual(a.substance, b.substance)
        && a.substanceStrength === b.substanceStrength;
}

export function hasInteriorThemeColor(props: { interior?: InteriorProps, [k: string]: any }): boolean {
    return !!props.interior?.themeColor;
}

export function getInteriorColor(props: InteriorProps, out: Vec4): Vec4 {
    Color.toArrayNormalized(props.color, out, 0);
    out[3] = props.colorStrength;
    return out;
}

export function getInteriorSubstance(props: InteriorProps, out: Vec4): Vec4 {
    Material.toArrayNormalized(props.substance, out, 0);
    out[3] = props.substanceStrength;
    return out;
}

export function createInteriorValues(props: InteriorProps) {
    return {
        uInteriorColor: ValueCell.create(getInteriorColor(props, Vec4())),
        uInteriorSubstance: ValueCell.create(getInteriorSubstance(props, Vec4())),
    };
}

export function updateInteriorValues(values: InteriorData, props: InteriorProps) {
    ValueCell.update(values.uInteriorColor, getInteriorColor(props, values.uInteriorColor.ref.value));
    ValueCell.update(values.uInteriorSubstance, getInteriorSubstance(props, values.uInteriorSubstance.ref.value));
}
