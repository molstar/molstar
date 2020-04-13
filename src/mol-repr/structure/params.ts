/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DirectVolume } from '../../mol-geo/geometry/direct-volume/direct-volume';
import { Lines } from '../../mol-geo/geometry/lines/lines';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { Points } from '../../mol-geo/geometry/points/points';
import { Spheres } from '../../mol-geo/geometry/spheres/spheres';
import { Text } from '../../mol-geo/geometry/text/text';
import { TextureMesh } from '../../mol-geo/geometry/texture-mesh/texture-mesh';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { UnitKind, UnitKindOptions } from './visual/util/common';

export function getUnitKindsParam(defaultValue: UnitKind[]) {
    return PD.MultiSelect<UnitKind>(defaultValue, UnitKindOptions, { description: 'For which kinds of units/chains to show the representation visuals.' });
}

export const StructureParams = {
    unitKinds: getUnitKindsParam(['atomic', 'spheres']),
};
export type StructureParams = typeof StructureParams

export const StructureMeshParams = { ...Mesh.Params };
export type StructureMeshParams = typeof StructureMeshParams

export const StructureSpheresParams = { ...Spheres.Params };
export type StructureSpheresParams = typeof StructureSpheresParams

export const StructurePointsParams = { ...Points.Params };
export type StructurePointsParams = typeof StructurePointsParams


export const StructureLinesParams = { ...Lines.Params };
export type StructureLinesParams = typeof StructureLinesParams

export const StructureTextParams = { ...Text.Params };
export type StructureTextParams = typeof StructureTextParams

export const StructureDirectVolumeParams = { ...DirectVolume.Params };
export type StructureDirectVolumeParams = typeof StructureDirectVolumeParams

export const StructureTextureMeshParams = { ...TextureMesh.Params };
export type StructureTextureMeshParams = typeof StructureTextureMeshParams