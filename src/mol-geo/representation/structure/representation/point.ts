/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UnitsRepresentation } from '..';
import { ElementPointVisual, DefaultElementPointProps } from '../visual/element-point';
import { StructureRepresentation } from '../units-representation';
import { SizeThemeProps } from 'mol-view/theme/size';

export const DefaultPointProps = {
    ...DefaultElementPointProps,

    sizeTheme: { name: 'uniform', value: 0.2 } as SizeThemeProps,
}
export type PointProps = typeof DefaultPointProps

export type PointRepresentation = StructureRepresentation<PointProps>

export function PointRepresentation(): PointRepresentation {
    return UnitsRepresentation('Point', ElementPointVisual)
}