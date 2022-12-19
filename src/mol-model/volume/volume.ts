/**
 * Copyright (c) 2020-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Grid } from './grid';
import { OrderedSet, SortedArray } from '../../mol-data/int';
import { Box3D, Sphere3D } from '../../mol-math/geometry';
import { Vec3, Mat4 } from '../../mol-math/linear-algebra';
import { BoundaryHelper } from '../../mol-math/geometry/boundary-helper';
import { CubeFormat } from '../../mol-model-formats/volume/cube';
import { equalEps } from '../../mol-math/linear-algebra/3d/common';
import { ModelFormat } from '../../mol-model-formats/format';
import { CustomProperties } from '../custom-property';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { toPrecision } from '../../mol-util/number';
import { DscifFormat } from '../../mol-model-formats/volume/density-server';

export interface Volume {
    readonly label?: string
    readonly entryId?: string,
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
    export function is(x: any): x is Volume {
        // TODO: improve
        return (
            x?.grid?.cells?.space?.dimensions?.length &&
            x?.sourceData &&
            x?.customProperties &&
            x?._propertyData
        );
    }

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

    // Converts iso value to relative if using downsample VolumeServer data
    export function adjustedIsoValue(volume: Volume, value: number, kind: 'absolute' | 'relative') {
        if (kind === 'relative') return IsoValue.relative(value);

        const absolute = IsoValue.absolute(value);
        if (DscifFormat.is(volume.sourceData)) {
            const stats = {
                min: volume.sourceData.data.volume_data_3d_info.min_source.value(0),
                max: volume.sourceData.data.volume_data_3d_info.max_source.value(0),
                mean: volume.sourceData.data.volume_data_3d_info.mean_source.value(0),
                sigma: volume.sourceData.data.volume_data_3d_info.sigma_source.value(0),
            };
            return Volume.IsoValue.toRelative(absolute, stats);
        }
        return absolute;
    }

    const defaultStats: Grid['stats'] = { min: -1, max: 1, mean: 0, sigma: 0.1 };
    export function createIsoValueParam(defaultValue: Volume.IsoValue, stats?: Grid['stats']) {
        const sts = stats || defaultStats;
        const { min, max, mean, sigma } = sts;

        // using ceil/floor could lead to "ouf of bounds" when converting
        const relMin = (min - mean) / sigma;
        const relMax = (max - mean) / sigma;

        let def = defaultValue;
        if (defaultValue.kind === 'absolute') {
            if (defaultValue.absoluteValue < min) def = Volume.IsoValue.absolute(min);
            else if (defaultValue.absoluteValue > max) def = Volume.IsoValue.absolute(max);
        } else {
            if (defaultValue.relativeValue < relMin) def = Volume.IsoValue.relative(relMin);
            else if (defaultValue.relativeValue > relMax) def = Volume.IsoValue.relative(relMax);
        }

        return PD.Conditioned(
            def,
            {
                'absolute': PD.Converted(
                    (v: Volume.IsoValue) => Volume.IsoValue.toAbsolute(v, Grid.One.stats).absoluteValue,
                    (v: number) => Volume.IsoValue.absolute(v),
                    PD.Numeric(mean, { min, max, step: toPrecision(sigma / 100, 2) }, { immediateUpdate: true })
                ),
                'relative': PD.Converted(
                    (v: Volume.IsoValue) => Volume.IsoValue.toRelative(v, Grid.One.stats).relativeValue,
                    (v: number) => Volume.IsoValue.relative(v),
                    PD.Numeric(Math.min(1, relMax), { min: relMin, max: relMax, step: toPrecision(Math.round(((max - min) / sigma)) / 100, 2) }, { immediateUpdate: true })
                )
            },
            (v: Volume.IsoValue) => v.kind === 'absolute' ? 'absolute' : 'relative',
            (v: Volume.IsoValue, c: 'absolute' | 'relative') => c === 'absolute' ? Volume.IsoValue.toAbsolute(v, sts) : Volume.IsoValue.toRelative(v, sts),
            { isEssential: true }
        );
    }

    export const IsoValueParam = createIsoValueParam(Volume.IsoValue.relative(2));
    export type IsoValueParam = typeof IsoValueParam

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
        return Grid.getBoundingSphere(volume.grid, boundingSphere);
    }

    export namespace Isosurface {
        export interface Loci { readonly kind: 'isosurface-loci', readonly volume: Volume, readonly isoValue: Volume.IsoValue }
        export function Loci(volume: Volume, isoValue: Volume.IsoValue): Loci { return { kind: 'isosurface-loci', volume, isoValue }; }
        export function isLoci(x: any): x is Loci { return !!x && x.kind === 'isosurface-loci'; }
        export function areLociEqual(a: Loci, b: Loci) { return a.volume === b.volume && Volume.IsoValue.areSame(a.isoValue, b.isoValue, a.volume.grid.stats); }
        export function isLociEmpty(loci: Loci) { return loci.volume.grid.cells.data.length === 0; }

        const bbox = Box3D();
        export function getBoundingSphere(volume: Volume, isoValue: Volume.IsoValue, boundingSphere?: Sphere3D) {
            const value = Volume.IsoValue.toAbsolute(isoValue, volume.grid.stats).absoluteValue;
            const neg = value < 0;

            const c = [0, 0, 0];
            const getCoords = volume.grid.cells.space.getCoords;
            const d = volume.grid.cells.data;
            const [xn, yn, zn] = volume.grid.cells.space.dimensions;

            let minx = xn - 1, miny = yn - 1, minz = zn - 1;
            let maxx = 0, maxy = 0, maxz = 0;
            for (let i = 0, il = d.length; i < il; ++i) {
                if ((neg && d[i] <= value) || (!neg && d[i] >= value)) {
                    getCoords(i, c);
                    if (c[0] < minx) minx = c[0];
                    if (c[1] < miny) miny = c[1];
                    if (c[2] < minz) minz = c[2];
                    if (c[0] > maxx) maxx = c[0];
                    if (c[1] > maxy) maxy = c[1];
                    if (c[2] > maxz) maxz = c[2];
                }
            }

            Vec3.set(bbox.min, minx - 1, miny - 1, minz - 1);
            Vec3.set(bbox.max, maxx + 1, maxy + 1, maxz + 1);
            const transform = Grid.getGridToCartesianTransform(volume.grid);
            Box3D.transform(bbox, bbox, transform);
            return Sphere3D.fromBox3D(boundingSphere || Sphere3D(), bbox);
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

    export namespace Segment {
        export interface Loci { readonly kind: 'segment-loci', readonly volume: Volume, readonly segments: SortedArray }
        export function Loci(volume: Volume, segments: ArrayLike<number>): Loci { return { kind: 'segment-loci', volume, segments: SortedArray.ofUnsortedArray(segments) }; }
        export function isLoci(x: any): x is Loci { return !!x && x.kind === 'segment-loci'; }
        export function areLociEqual(a: Loci, b: Loci) { return a.volume === b.volume && SortedArray.areEqual(a.segments, b.segments); }
        export function isLociEmpty(loci: Loci) { return loci.volume.grid.cells.data.length === 0 || loci.segments.length === 0; }

        const bbox = Box3D();
        export function getBoundingSphere(volume: Volume, segments: ArrayLike<number>, boundingSphere?: Sphere3D) {
            const segmentation = Volume.Segmentation.get(volume);
            if (segmentation) {
                Box3D.setEmpty(bbox);
                for (let i = 0, il = segments.length; i < il; ++i) {
                    const b = segmentation.bounds[segments[i]];
                    Box3D.add(bbox, b.min);
                    Box3D.add(bbox, b.max);
                }
                const transform = Grid.getGridToCartesianTransform(volume.grid);
                Box3D.transform(bbox, bbox, transform);
                return Sphere3D.fromBox3D(boundingSphere || Sphere3D(), bbox);
            } else {
                return Volume.getBoundingSphere(volume, boundingSphere);
            }
        }

        export interface Location {
            readonly kind: 'segment-location',
            volume: Volume
            segment: number
        }
        export function Location(volume?: Volume, segment?: number): Location {
            return { kind: 'segment-location', volume: volume as any, segment: segment as any };
        }
        export function isLocation(x: any): x is Location {
            return !!x && x.kind === 'segment-location';
        }
    }

    export type PickingGranularity = 'volume' | 'object' | 'voxel';
    export const PickingGranularity = {
        set(volume: Volume, granularity: PickingGranularity) {
            volume._propertyData['__picking_granularity__'] = granularity;
        },
        get(volume: Volume): PickingGranularity {
            return volume._propertyData['__picking_granularity__'] ?? 'voxel';
        }
    };

    export type Segmentation = {
        segments: Map<number, Set<number>>
        sets: Map<number, Set<number>>
        bounds: { [k: number]: Box3D }
        labels: { [k: number]: string }
    };
    export const Segmentation = {
        set(volume: Volume, segmentation: Segmentation) {
            volume._propertyData['__segmentation__'] = segmentation;
        },
        get(volume: Volume): Segmentation | undefined {
            return volume._propertyData['__segmentation__'];
        }
    };
}