/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { VolumeData, VolumeIsoValue } from './data';
import { OrderedSet } from '../../mol-data/int';
import { Sphere3D } from '../../mol-math/geometry';
import { Vec3 } from '../../mol-math/linear-algebra';

export namespace Volume {
    export type CellIndex = { readonly '@type': 'cell-index' } & number

    export interface Loci { readonly kind: 'volume-loci', readonly volume: VolumeData }
    export function Loci(volume: VolumeData): Loci { return { kind: 'volume-loci', volume }; }
    export function isLoci(x: any): x is Loci { return !!x && x.kind === 'volume-loci'; }
    export function areLociEqual(a: Loci, b: Loci) { return a.volume === b.volume; }
    export function isLociEmpty(loci: Loci) { return loci.volume.data.data.length === 0; }

    export function getBoundingSphere(volume: VolumeData, boundingSphere?: Sphere3D) {
        if (!boundingSphere) boundingSphere = Sphere3D();

        const transform = VolumeData.getGridToCartesianTransform(volume);
        const [x, y, z] = volume.data.space.dimensions;

        const cpA = Vec3.create(0, 0, 0); Vec3.transformMat4(cpA, cpA, transform);
        const cpB = Vec3.create(x, y, z); Vec3.transformMat4(cpB, cpB, transform);
        const cpC = Vec3.create(x, 0, 0); Vec3.transformMat4(cpC, cpC, transform);
        const cpD = Vec3.create(0, y, z); Vec3.transformMat4(cpD, cpC, transform);

        const cpE = Vec3.create(0, 0, z); Vec3.transformMat4(cpE, cpE, transform);
        const cpF = Vec3.create(x, 0, z); Vec3.transformMat4(cpF, cpF, transform);
        const cpG = Vec3.create(x, y, 0); Vec3.transformMat4(cpG, cpG, transform);
        const cpH = Vec3.create(0, y, 0); Vec3.transformMat4(cpH, cpH, transform);

        const center = Vec3();
        Vec3.add(center, cpA, cpB);
        Vec3.scale(center, center, 0.5);
        const d = Math.max(Vec3.distance(cpA, cpB), Vec3.distance(cpC, cpD));
        Sphere3D.set(boundingSphere, center, d / 2);
        Sphere3D.setExtrema(boundingSphere, [cpA, cpB, cpC, cpD, cpE, cpF, cpG, cpH]);

        return boundingSphere;
    }

    export namespace Isosurface {
        export interface Loci { readonly kind: 'isosurface-loci', readonly volume: VolumeData, readonly isoValue: VolumeIsoValue }
        export function Loci(volume: VolumeData, isoValue: VolumeIsoValue): Loci { return { kind: 'isosurface-loci', volume, isoValue }; }
        export function isLoci(x: any): x is Loci { return !!x && x.kind === 'isosurface-loci'; }
        export function areLociEqual(a: Loci, b: Loci) { return a.volume === b.volume && VolumeIsoValue.areSame(a.isoValue, b.isoValue, a.volume.dataStats); }
        export function isLociEmpty(loci: Loci) { return loci.volume.data.data.length === 0; }

        export function getBoundingSphere(volume: VolumeData, boundingSphere?: Sphere3D) {
            return Volume.getBoundingSphere(volume, boundingSphere);
        }
    }

    export namespace Cell {
        export interface Loci { readonly kind: 'cell-loci', readonly volume: VolumeData, readonly indices: OrderedSet<CellIndex> }
        export function Loci(volume: VolumeData, indices: OrderedSet<CellIndex>): Loci { return { kind: 'cell-loci', volume, indices }; }
        export function isLoci(x: any): x is Loci { return !!x && x.kind === 'cell-loci'; }
        export function areLociEqual(a: Loci, b: Loci) { return a.volume === b.volume && OrderedSet.areEqual(a.indices, b.indices); }
        export function isLociEmpty(loci: Loci) { return OrderedSet.size(loci.indices) === 0; }

        export function getBoundingSphere(volume: VolumeData, boundingSphere?: Sphere3D) {
            if (!boundingSphere) boundingSphere = Sphere3D();
            // TODO get sphere reflecting cell coords and size
            return Volume.getBoundingSphere(volume, boundingSphere);
        }
    }
}