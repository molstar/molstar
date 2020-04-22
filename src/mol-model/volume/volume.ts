/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { VolumeData, VolumeIsoValue } from './data';
import { OrderedSet } from '../../mol-data/int';

export namespace Volume {
    export type CellIndex = { readonly '@type': 'cell-index' } & number

    export interface Loci { readonly kind: 'volume-loci', readonly volume: VolumeData }
    export function Loci(volume: VolumeData): Loci { return { kind: 'volume-loci', volume }; }
    export function isLoci(x: any): x is Loci { return !!x && x.kind === 'volume-loci'; }
    export function areLociEqual(a: Loci, b: Loci) { return a.volume === b.volume; }
    export function isLociEmpty(loci: Loci) { return loci.volume.data.data.length === 0; }

    export namespace Isosurface {
        export interface Loci { readonly kind: 'isosurface-loci', readonly volume: VolumeData, readonly isoValue: VolumeIsoValue }
        export function Loci(volume: VolumeData, isoValue: VolumeIsoValue): Loci { return { kind: 'isosurface-loci', volume, isoValue }; }
        export function isLoci(x: any): x is Loci { return !!x && x.kind === 'isosurface-loci'; }
        export function areLociEqual(a: Loci, b: Loci) { return a.volume === b.volume && VolumeIsoValue.areSame(a.isoValue, b.isoValue, a.volume.dataStats); }
        export function isLociEmpty(loci: Loci) { return loci.volume.data.data.length === 0; }
    }

    export namespace Cell {
        export interface Loci { readonly kind: 'cell-loci', readonly volume: VolumeData, readonly indices: OrderedSet<CellIndex> }
        export function Loci(volume: VolumeData, indices: OrderedSet<CellIndex>): Loci { return { kind: 'cell-loci', volume, indices }; }
        export function isLoci(x: any): x is Loci { return !!x && x.kind === 'cell-loci'; }
        export function areLociEqual(a: Loci, b: Loci) { return a.volume === b.volume && OrderedSet.areEqual(a.indices, b.indices); }
        export function isLociEmpty(loci: Loci) { return OrderedSet.size(loci.indices) === 0; }
    }
}