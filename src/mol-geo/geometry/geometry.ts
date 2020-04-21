/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mesh } from './mesh/mesh';
import { Points } from './points/points';
import { Text } from './text/text';
import { RenderableState } from '../../mol-gl/renderable';
import { LocationIterator } from '../util/location-iterator';
import { ColorType } from './color-data';
import { SizeType } from './size-data';
import { Lines } from './lines/lines';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { DirectVolume } from './direct-volume/direct-volume';
import { Color } from '../../mol-util/color';
import { Spheres } from './spheres/spheres';
import { arrayMax } from '../../mol-util/array';
import { TransformData } from './transform-data';
import { Theme } from '../../mol-theme/theme';
import { RenderObjectValues } from '../../mol-gl/render-object';
import { TextureMesh } from './texture-mesh/texture-mesh';
import { Image } from './image/image';

export type GeometryKind = 'mesh' | 'points' | 'spheres' | 'text' | 'lines' | 'direct-volume' | 'image' | 'texture-mesh'

export type Geometry<T extends GeometryKind = GeometryKind> =
    T extends 'mesh' ? Mesh :
        T extends 'points' ? Points :
            T extends 'spheres' ? Spheres :
                T extends 'text' ? Text :
                    T extends 'lines' ? Lines :
                        T extends 'direct-volume' ? DirectVolume :
                            T extends 'image' ? Image :
                                T extends 'texture-mesh' ? TextureMesh : never

type GeometryParams<T extends GeometryKind> =
    T extends 'mesh' ? Mesh.Params :
        T extends 'points' ? Points.Params :
            T extends 'spheres' ? Spheres.Params :
                T extends 'text' ? Text.Params :
                    T extends 'lines' ? Lines.Params :
                        T extends 'direct-volume' ? DirectVolume.Params :
                            T extends 'image' ? Image.Params :
                                T extends 'texture-mesh' ? TextureMesh.Params : never

export interface GeometryUtils<G extends Geometry, P extends PD.Params = GeometryParams<G['kind']>, V = RenderObjectValues<G['kind']>> {
    Params: P
    createEmpty(geometry?: G): G
    createValues(geometry: G, transform: TransformData, locationIt: LocationIterator, theme: Theme, props: PD.Values<P>): V
    createValuesSimple(geometry: G, props: Partial<PD.Values<P>>, colorValue: Color, sizeValue: number, transform?: TransformData): V
    updateValues(values: V, props: PD.Values<P>): void
    updateBoundingSphere(values: V, geometry: G): void
    createRenderableState(props: Partial<PD.Values<P>>): RenderableState
    updateRenderableState(state: RenderableState, props: PD.Values<P>): void
}

export namespace Geometry {
    export type Params<G extends Geometry> = GeometryParams<G['kind']>

    export function getDrawCount(geometry: Geometry): number {
        switch (geometry.kind) {
            case 'mesh': return geometry.triangleCount * 3;
            case 'points': return geometry.pointCount;
            case 'spheres': return geometry.sphereCount * 2 * 3;
            case 'text': return geometry.charCount * 2 * 3;
            case 'lines': return geometry.lineCount * 2 * 3;
            case 'direct-volume': return 12 * 3;
            case 'image': return 2 * 3;
            case 'texture-mesh': return geometry.vertexCount;
        }
    }

    export function getGroupCount(geometry: Geometry): number {
        switch (geometry.kind) {
            case 'mesh':
            case 'points':
            case 'spheres':
            case 'text':
            case 'lines':
                return getDrawCount(geometry) === 0 ? 0 : (arrayMax(geometry.groupBuffer.ref.value) + 1);
            case 'direct-volume':
                return 1;
            case 'image':
                return arrayMax(geometry.groupTexture.ref.value.array) + 1;
            case 'texture-mesh':
                return geometry.groupCount;
        }
    }

    export function getUtils<G extends Geometry>(geometry: G): GeometryUtils<G> {
        // TODO avoid casting
        switch (geometry.kind) {
            case 'mesh': return Mesh.Utils as any;
            case 'points': return Points.Utils as any;
            case 'spheres': return Spheres.Utils as any;
            case 'text': return Text.Utils as any;
            case 'lines': return Lines.Utils as any;
            case 'direct-volume': return DirectVolume.Utils as any;
            case 'image': return Image.Utils as any;
            case 'texture-mesh': return TextureMesh.Utils as any;
        }
    }

    export function getGranularity(locationIt: LocationIterator, granularity: ColorType | SizeType) {
        // Always use 'group' granularity for 'complex' location iterators,
        // i.e. for which an instance may include multiple units
        return granularity === 'instance' && locationIt.isComplex ? 'group' : granularity;
    }
}