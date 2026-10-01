/**
 * Copyright (c) 2025-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { OrderedSet } from '@molstar/core/data/int/ordered-set';
import type { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { DataLocation } from '@molstar/model/model/location';
import { DataLoci } from '@molstar/model/model/loci';
import { Grid } from '@molstar/model/model/volume/grid';
import type { Volume } from '@molstar/model/model/volume/volume';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { BoundaryHelper } from '@molstar/core/math/geometry/boundary-helper';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import type { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';

export type StreamlinePoint = Vec3
export type Streamline = StreamlinePoint[]
export type Streamlines = Streamline[]
export type StreamlinesIndex = { readonly '@type': 'streamlines-index' } & number

//

type VolumeStreamlines = { readonly volume: Volume, readonly streamlines: Streamlines }

export interface StreamlinesLocation extends DataLocation<VolumeStreamlines, { index: StreamlinesIndex, instance: Volume.InstanceIndex }> {}

export function StreamlinesLocation(streamlines: Streamlines, volume: Volume, index?: StreamlinesIndex, instance?: Volume.InstanceIndex): StreamlinesLocation {
    return DataLocation('streamlines', { volume, streamlines }, { index: index as any, instance: instance as any });
}

export function isStreamlinesLocation(x: any): x is StreamlinesLocation {
    return !!x && x.kind === 'data-location' && x.tag === 'streamlines';
}

export function areStreamlinesLocationsEqual(locA: StreamlinesLocation, locB: StreamlinesLocation) {
    return (
        locA.data.volume === locB.data.volume &&
        locA.data.streamlines === locB.data.streamlines &&
        locA.element.index === locB.element.index &&
        locA.element.instance === locB.element.instance
    );
}

export function streamlinesLocationLabel(streamlines: Streamlines, volume: Volume, index: StreamlinesIndex, instance: Volume.InstanceIndex): string {
    const label = [
        `${volume.label || 'Volume'}`,
        `Streamline #${index}`
    ];
    if (volume.instances.length > 1) {
        label.push(`Instance #${instance}`);
    }
    return label.join(' | ');
}

export function streamlinesLociLabel(streamlines: Streamlines, volume: Volume, elements: ReadonlyArray<StreamlinesElement>): string {
    const size = getStreamlinesLociSize(elements);
    const label = [
        `${volume.label || 'Volume'}`
    ];
    if (size === 0) {
        label.push('No Streamlines');
    } else if (size === 1) {
        const index = OrderedSet.start(elements[0].indices);
        label.push(`Streamline #${index}`);
        if (volume.instances.length > 1) {
            const instance = OrderedSet.start(elements[0].instances);
            label.push(`Instance #${instance}`);
        }
    } else {
        label.push(`${size} Streamlines`);
    }
    return label.join(' | ');
}

type StreamlinesElement = {
    readonly indices: OrderedSet<StreamlinesIndex>,
    readonly instances: OrderedSet<Volume.InstanceIndex>
}

export interface StreamlinesLoci extends DataLoci<VolumeStreamlines, StreamlinesElement> { }

export function StreamlinesLoci(streamlines: Streamlines, volume: Volume, elements: ReadonlyArray<StreamlinesElement>): StreamlinesLoci {
    return DataLoci('streamlines', { streamlines, volume }, elements,
        (boundingSphere) => getStreamlinesLociBoundingSphere(streamlines, volume, elements, boundingSphere),
        () => streamlinesLociLabel(streamlines, volume, elements));
}

export function isStreamlinesLoci(x: any): x is StreamlinesLoci {
    return !!x && x.kind === 'data-loci' && x.tag === 'streamlines';
}

export function getStreamlinesLociSize(elements: StreamlinesLoci['elements']): number {
    let size = 0;
    for (const e of elements) {
        size += OrderedSet.size(e.indices) * OrderedSet.size(e.instances);
    }
    return size;
}

const boundaryHelper = new BoundaryHelper('98');
const tmpBoundaryPos = Vec3();
const tmpBoundaryPos2 = Vec3();
export function getStreamlinesLociBoundingSphere(streamlines: Streamlines, volume: Volume, elements: StreamlinesLoci['elements'], boundingSphere?: Sphere3D) {
    boundaryHelper.reset();
    const transform = Grid.getGridToCartesianTransform(volume.grid);

    for (const { indices, instances } of elements) {
        for (let i = 0, _i = OrderedSet.size(indices); i < _i; i++) {
            const o = OrderedSet.getAt(indices, i);
            for (const p of streamlines[o]) {
                Vec3.transformMat4(tmpBoundaryPos, p, transform);
                for (let j = 0, _j = OrderedSet.size(instances); j < _j; j++) {
                    const instance = volume.instances[OrderedSet.getAt(instances, j)];
                    Vec3.transformMat4(tmpBoundaryPos2, tmpBoundaryPos, instance.transform);
                    boundaryHelper.includePosition(tmpBoundaryPos2);
                }
            }
        }
    }
    boundaryHelper.finishedIncludeStep();
    for (const { indices, instances } of elements) {
        for (let i = 0, _i = OrderedSet.size(indices); i < _i; i++) {
            const o = OrderedSet.getAt(indices, i);
            for (const p of streamlines[o]) {
                Vec3.transformMat4(tmpBoundaryPos, p, transform);
                for (let j = 0, _j = OrderedSet.size(instances); j < _j; j++) {
                    const instance = volume.instances[OrderedSet.getAt(instances, j)];
                    Vec3.transformMat4(tmpBoundaryPos2, tmpBoundaryPos, instance.transform);
                    boundaryHelper.radiusPosition(tmpBoundaryPos2);
                }
            }
        }
    }

    return boundaryHelper.getSphere(boundingSphere);
}

export function areStreamlinesLociEqual(a: StreamlinesLoci, b: StreamlinesLoci) {
    if (a.data.volume !== b.data.volume || a.elements.length !== b.elements.length) return false;
    for (let i = 0, il = a.elements.length; i < il; ++i) {
        const ae = a.elements[i], be = b.elements[i];
        if (!OrderedSet.areEqual(ae.instances, be.instances) ||
            !OrderedSet.areEqual(ae.indices, be.indices)) return false;
    }
    return true;
}

export function isStreamlinesLociEmpty(loci: StreamlinesLoci) {
    for (const { indices, instances } of loci.elements) {
        if (!OrderedSet.isEmpty(indices) || !OrderedSet.isEmpty(instances)) return false;
    }
    return true;
}

//

export const CommonStreamlinesParams = {
    anchorEnabled: PD.Boolean(false, { label: 'Anchor' }),
    anchorCenter: PD.Vec3(Vec3(), undefined, { hideIf: p => !p.anchorEnabled }),
    anchorRadius: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }, { hideIf: p => !p.anchorEnabled }),
    dashEnabled: PD.Boolean(false, { label: 'Dash' }),
    dashPoints: PD.Numeric(5, { min: 1, max: 25, step: 1 }, { description: 'Number of streamline points per dash/gap', hideIf: p => !p.dashEnabled }),
    dashShift: PD.Boolean(false, { description: 'Shift dashes so dashes become gaps and vice versa', hideIf: p => !p.dashEnabled }),
};
export type CommonStreamlinesParams = typeof CommonStreamlinesParams
export type CommonStreamlinesProps = PD.Values<CommonStreamlinesParams>

const tmpFilterVec = Vec3();
/**
 * Check if a streamline passes the sphere filter.
 * A streamline passes if its start or end point is within the filter sphere.
 */
export function streamlinePassesFilter(streamline: Streamline, gridToCartn: Mat4, props: CommonStreamlinesProps): boolean {
    if (streamline.length < 2) return false;
    if (!props.anchorEnabled) return true;

    // Check start point
    const start = streamline[0];
    Vec3.transformMat4(tmpFilterVec, start, gridToCartn);
    if (Vec3.distance(tmpFilterVec, props.anchorCenter) <= props.anchorRadius) return true;

    // Check end point
    const end = streamline[streamline.length - 1];
    Vec3.transformMat4(tmpFilterVec, end, gridToCartn);
    if (Vec3.distance(tmpFilterVec, props.anchorCenter) <= props.anchorRadius) return true;

    return false;
}
