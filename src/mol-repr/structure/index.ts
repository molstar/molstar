/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from 'mol-model/structure';
import { ColorThemeName, ColorThemeOptions } from 'mol-theme/color';
import { SizeThemeName, SizeThemeOptions } from 'mol-theme/size';
import { Representation, RepresentationProps } from '..';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Geometry } from 'mol-geo/geometry/geometry';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { Points } from 'mol-geo/geometry/points/points';
import { Lines } from 'mol-geo/geometry/lines/lines';
import { DirectVolume } from 'mol-geo/geometry/direct-volume/direct-volume';

export interface StructureRepresentation<P extends RepresentationProps = {}> extends Representation<Structure, P> { }
// export interface  StructureVisual<P extends RepresentationProps = {}> extends Visual<Structure, P> { }

export const StructureParams = {
    ...Geometry.Params,
    colorTheme: PD.SelectParam<ColorThemeName>('Color Theme', '', 'polymer-index', ColorThemeOptions),
    sizeTheme: PD.SelectParam<SizeThemeName>('Size Theme', '', 'physical', SizeThemeOptions),
}
export const DefaultStructureProps = PD.paramDefaultValues(StructureParams)
export type StructureProps = typeof DefaultStructureProps

export const StructureMeshParams = {
    ...Mesh.Params,
    ...StructureParams,
}
export const DefaultStructureMeshProps = PD.paramDefaultValues(StructureMeshParams)
export type StructureMeshProps = typeof DefaultStructureMeshProps

export const StructurePointsParams = {
    ...Points.Params,
    ...StructureParams,
}
export const DefaultStructurePointsProps = PD.paramDefaultValues(StructurePointsParams)
export type StructurePointsProps = typeof DefaultStructurePointsProps

export const StructureLinesParams = {
    ...Lines.Params,
    ...StructureParams,
}
export const DefaultStructureLinesProps = PD.paramDefaultValues(StructureLinesParams)
export type StructureLinesProps = typeof DefaultStructureLinesProps

export const StructureDirectVolumeParams = {
    ...DirectVolume.Params,
    ...StructureParams,
}
export const DefaultStructureDirectVolumeProps = PD.paramDefaultValues(StructureDirectVolumeParams)
export type StructureDirectVolumeProps = typeof DefaultStructureDirectVolumeProps

export { ComplexRepresentation } from './complex-representation'
export { UnitsRepresentation } from './units-representation'
export { ComplexVisual } from './complex-visual'
export { UnitsVisual } from './units-visual'