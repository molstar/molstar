/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from 'mol-model/structure';
import { ColorThemeProps } from 'mol-view/theme/color';
import { SizeThemeProps } from 'mol-view/theme/size';
import { Representation, RepresentationProps } from '..';
import { Geometry } from '../../geometry/geometry';
import { Mesh } from '../../geometry/mesh/mesh';
import { Point } from '../../geometry/point/point';

export interface StructureRepresentation<P extends RepresentationProps = {}> extends Representation<Structure, P> { }

export const DefaultStructureProps = {
    ...Geometry.DefaultProps,
    colorTheme: { name: 'unit-index' } as ColorThemeProps,
    sizeTheme: { name: 'physical' } as SizeThemeProps,
}
export type StructureProps = typeof DefaultStructureProps

export const DefaultStructureMeshProps = {
    ...Mesh.DefaultProps,
    ...DefaultStructureProps,
}
export type StructureMeshProps = typeof DefaultStructureMeshProps

export const DefaultStructurePointProps = {
    ...Point.DefaultProps,
    ...DefaultStructureProps,
}
export type StructurePointProps = typeof DefaultStructurePointProps

export interface VisualUpdateState {
    updateTransform: boolean
    updateColor: boolean
    updateSize: boolean
    createGeometry: boolean
}
export namespace VisualUpdateState {
    export function create(): VisualUpdateState {
        return {
            updateTransform: false,
            updateColor: false,
            updateSize: false,
            createGeometry: false
        }
    }
    export function reset(state: VisualUpdateState) {
        state.updateTransform = false
        state.updateColor = false
        state.updateSize = false
        state.createGeometry = false
    }
}

export { ComplexRepresentation } from './complex-representation'
export { UnitsRepresentation } from './units-representation'
export { ComplexVisual } from './complex-visual'
export { UnitsVisual } from './units-visual'