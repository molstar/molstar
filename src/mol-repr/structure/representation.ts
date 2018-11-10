/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from 'mol-model/structure';
import { Representation, RepresentationProps, RepresentationProvider } from '../representation';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Geometry } from 'mol-geo/geometry/geometry';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { Points } from 'mol-geo/geometry/points/points';
import { Lines } from 'mol-geo/geometry/lines/lines';
import { DirectVolume } from 'mol-geo/geometry/direct-volume/direct-volume';

export interface StructureRepresentation<P extends RepresentationProps = {}> extends Representation<Structure, P> { }

export type StructureRepresentationProvider<P extends PD.Params> = RepresentationProvider<Structure, P>

//

export const StructureParams = {
    ...Geometry.Params,
}
export const DefaultStructureProps = PD.getDefaultValues(StructureParams)
export type StructureProps = typeof DefaultStructureProps

export const StructureMeshParams = {
    ...Mesh.Params,
    ...StructureParams,
}
export const DefaultStructureMeshProps = PD.getDefaultValues(StructureMeshParams)
export type StructureMeshProps = typeof DefaultStructureMeshProps

export const StructurePointsParams = {
    ...Points.Params,
    ...StructureParams,
}
export const DefaultStructurePointsProps = PD.getDefaultValues(StructurePointsParams)
export type StructurePointsProps = typeof DefaultStructurePointsProps

export const StructureLinesParams = {
    ...Lines.Params,
    ...StructureParams,
}
export const DefaultStructureLinesProps = PD.getDefaultValues(StructureLinesParams)
export type StructureLinesProps = typeof DefaultStructureLinesProps

export const StructureDirectVolumeParams = {
    ...DirectVolume.Params,
    ...StructureParams,
}
export const DefaultStructureDirectVolumeProps = PD.getDefaultValues(StructureDirectVolumeParams)
export type StructureDirectVolumeProps = typeof DefaultStructureDirectVolumeProps

export { ComplexRepresentation } from './complex-representation'
export { UnitsRepresentation } from './units-representation'
export { ComplexVisual } from './complex-visual'
export { UnitsVisual } from './units-visual'