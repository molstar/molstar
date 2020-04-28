/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Grid } from './grid';
import { OrderedSet } from '../../mol-data/int';
import { Sphere3D } from '../../mol-math/geometry';
import { Vec3, Mat4 } from '../../mol-math/linear-algebra';
import { BoundaryHelper } from '../../mol-math/geometry/boundary-helper';
import { CubeFormat } from '../../mol-model-formats/volume/cube';
import { equalEps } from '../../mol-math/linear-algebra/3d/common';
import { ModelFormat } from '../../mol-model-formats/format';
import { CustomProperties } from '../custom-property';

export interface Volume {
    readonly label?: string
    readonly grid: Grid
    readonly sourceData: ModelFormat

    // TODO use...
    customProperties: CustomProperties

    /**
     * Not to be accessed directly, each custom property descriptor
     * defines property accessors that use this field to store the data.
     */
    _propertyData: { [name: string]: any }

    // TODO add as customProperty?
    readonly colorVolume?: Volume
}

export namespace Volume {
    export type CellIndex = { readonly '@type': 'cell-index' } & number

    export type IsoValue = IsoValue.Absolute | IsoValue.Relative

    export namespace IsoValue {
        export type Relative = Readonly<{ kind: 'relative', relativeValue: number }>
        export type Absolute = Readonly<{ kind: 'absolute', absoluteValue: number }>

        export function areSame(a: IsoValue, b: IsoValue, stats: Grid['stats']) {
            return equalEps(toAbsolute(a, stats).absoluteValue, toAbsolute(b, stats).absoluteValue, stats.sigma / 100);
        }

        export function absolute(value: number): Absolute { return { kind: 'absolute', absoluteValue: value }; }
        export function relative(value: number): Relative { return { kind: 'relative', relativeValue: value }; }

        export function calcAbsolute(stats: Grid['stats'], relativeValue: number): number {
            return relativeValue * stats.sigma + stats.mean;
        }

        export function calcRelative(stats: Grid['stats'], absoluteValue: number): number {
            return stats.sigma === 0 ? 0 : ((absoluteValue - stats.mean) / stats.sigma);
        }

        export function toAbsolute(value: IsoValue, stats: Grid['stats']): Absolute {
            return value.kind === 'absolute' ? value : { kind: 'absolute', absoluteValue: IsoValue.calcAbsolute(stats, value.relativeValue) };
        }

        export function toRelative(value: IsoValue, stats: Grid['stats']): Relative {
            return value.kind === 'relative' ? value : { kind: 'relative', relativeValue: IsoValue.calcRelative(stats, value.absoluteValue) };
        }

        export function toString(value: IsoValue) {
            return value.kind === 'relative'
                ? `${value.relativeValue.toFixed(2)} σ`
                : `${value.absoluteValue.toPrecision(4)}`;
        }
    }

    export const One: Volume = {
        label: '',
        grid: Grid.One,
        sourceData: { kind: '', name: '', data: {} },
        customProperties: new CustomProperties(),
        _propertyData: Object.create(null),
    };

    export function areEquivalent(volA: Volume, volB: Volume) {
        return Grid.areEquivalent(volA.grid, volB.grid);
    }

    export function isEmpty(vol: Volume) {
        return Grid.isEmpty(vol.grid);
    }

    export function isOrbitals(volume: Volume) {
        if (!CubeFormat.is(volume.sourceData)) return false;
        return volume.sourceData.data.header.orbitals;
    }

    export interface Loci { readonly kind: 'volume-loci', readonly volume: Volume }
    export function Loci(volume: Volume): Loci { return { kind: 'volume-loci', volume }; }
    export function isLoci(x: any): x is Loci { return !!x && x.kind === 'volume-loci'; }
    export function areLociEqual(a: Loci, b: Loci) { return a.volume === b.volume; }
    export function isLociEmpty(loci: Loci) { return Grid.isEmpty(loci.volume.grid); }

    export function getBoundingSphere(volume: Volume, boundingSphere?: Sphere3D) {
        if (!boundingSphere) boundingSphere = Sphere3D();

        const transform = Grid.getGridToCartesianTransform(volume.grid);
        const [x, y, z] = volume.grid.cells.space.dimensions;

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
        export interface Loci { readonly kind: 'isosurface-loci', readonly volume: Volume, readonly isoValue: Volume.IsoValue }
        export function Loci(volume: Volume, isoValue: Volume.IsoValue): Loci { return { kind: 'isosurface-loci', volume, isoValue }; }
        export function isLoci(x: any): x is Loci { return !!x && x.kind === 'isosurface-loci'; }
        export function areLociEqual(a: Loci, b: Loci) { return a.volume === b.volume && Volume.IsoValue.areSame(a.isoValue, b.isoValue, a.volume.grid.stats); }
        export function isLociEmpty(loci: Loci) { return loci.volume.grid.cells.data.length === 0; }

        export function getBoundingSphere(volume: Volume, isoValue: Volume.IsoValue, boundingSphere?: Sphere3D) {
            // TODO get bounding sphere for subgrid with values >= isoValue
            return Volume.getBoundingSphere(volume, boundingSphere);
        }
    }

    export namespace Cell {
        export interface Loci { readonly kind: 'cell-loci', readonly volume: Volume, readonly indices: OrderedSet<CellIndex> }
        export function Loci(volume: Volume, indices: OrderedSet<CellIndex>): Loci { return { kind: 'cell-loci', volume, indices }; }
        export function isLoci(x: any): x is Loci { return !!x && x.kind === 'cell-loci'; }
        export function areLociEqual(a: Loci, b: Loci) { return a.volume === b.volume && OrderedSet.areEqual(a.indices, b.indices); }
        export function isLociEmpty(loci: Loci) { return OrderedSet.size(loci.indices) === 0; }

        const boundaryHelper = new BoundaryHelper('98');
        const tmpBoundaryPos = Vec3();
        export function getBoundingSphere(volume: Volume, indices: OrderedSet<CellIndex>, boundingSphere?: Sphere3D) {
            boundaryHelper.reset();
            const transform = Grid.getGridToCartesianTransform(volume.grid);
            const { getCoords } = volume.grid.cells.space;

            for (let i = 0, _i = OrderedSet.size(indices); i < _i; i++) {
                const o = OrderedSet.getAt(indices, i);
                getCoords(o, tmpBoundaryPos);
                Vec3.transformMat4(tmpBoundaryPos, tmpBoundaryPos, transform);
                boundaryHelper.includePosition(tmpBoundaryPos);
            }
            boundaryHelper.finishedIncludeStep();
            for (let i = 0, _i = OrderedSet.size(indices); i < _i; i++) {
                const o = OrderedSet.getAt(indices, i);
                getCoords(o, tmpBoundaryPos);
                Vec3.transformMat4(tmpBoundaryPos, tmpBoundaryPos, transform);
                boundaryHelper.radiusPosition(tmpBoundaryPos);
            }

            const bs = boundaryHelper.getSphere(boundingSphere);
            return Sphere3D.expand(bs, bs, Mat4.getMaxScaleOnAxis(transform) * 10);
        }
    }
}