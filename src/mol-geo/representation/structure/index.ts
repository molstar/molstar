/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from 'mol-model/structure';
import { ColorThemeName, ColorThemeOptions } from 'mol-view/theme/color';
import { SizeThemeName, SizeThemeOptions } from 'mol-view/theme/size';
import { Representation, RepresentationProps } from '..';
import { Geometry } from '../../geometry/geometry';
import { Mesh } from '../../geometry/mesh/mesh';
import { Points } from '../../geometry/points/points';
import { Lines } from '../../geometry/lines/lines';
import { SelectParam, paramDefaultValues } from 'mol-view/parameter';

export interface StructureRepresentation<P extends RepresentationProps = {}> extends Representation<Structure, P> { }

export const StructureParams = {
    ...Geometry.Params,
    colorTheme: SelectParam<ColorThemeName>('Color Theme', '', 'unit-index', ColorThemeOptions),
    sizeTheme: SelectParam<SizeThemeName>('Size Theme', '', 'physical', SizeThemeOptions),
}
export const DefaultStructureProps = paramDefaultValues(StructureParams)
export type StructureProps = typeof DefaultStructureProps

export const StructureMeshParams = {
    ...Mesh.Params,
    ...StructureParams,
}
export const DefaultStructureMeshProps = paramDefaultValues(StructureMeshParams)
export type StructureMeshProps = typeof DefaultStructureMeshProps

export const StructurePointsParams = {
    ...Points.Params,
    ...StructureParams,
}
export const DefaultStructurePointsProps = paramDefaultValues(StructurePointsParams)
export type StructurePointsProps = typeof DefaultStructurePointsProps

export const StructureLinesParams = {
    ...Lines.Params,
    ...StructureParams,
}
export const DefaultStructureLinesProps = paramDefaultValues(StructureLinesParams)
export type StructureLinesProps = typeof DefaultStructureLinesProps

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