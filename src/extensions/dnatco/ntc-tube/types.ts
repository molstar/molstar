/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { NtCTubeSegmentLabel } from './behavior';
import { DnatcoTypes } from '../types';
import { Sphere3D } from '../../../mol-math/geometry/primitives/sphere3d';
import { DataLocation } from '../../../mol-model/location';
import { DataLoci } from '../../../mol-model/loci';

export namespace NtCTubeTypes {
    const DataTag = 'dnatco-tube-segment-data';
    const DummyTag = 'dnatco-tube-dummy';

    export type Data = {
        data: DnatcoTypes.Steps,
    }

    export type TubeBlock = {
        step: DnatcoTypes.Step,
        kind: 'upper' | 'lower' | 'residue-boundary' | 'segment-boundary';
    }

    export interface Location extends DataLocation<TubeBlock> {}

    export function Location(payload: TubeBlock) {
        return DataLocation(DataTag, payload, {});
    }

    export function isLocation(x: any): x is Location {
        return !!x && x.kind === 'data-location' && x.tag === DataTag;
    }

    export interface Loci extends DataLoci<DnatcoTypes.Step[], number> {}
    export interface DummyLoci extends DataLoci<{}, number> {}

    export function Loci(data: DnatcoTypes.Step[], stepIndices: number[], elements: number[], boundingSphere?: Sphere3D): Loci {
        return DataLoci(DataTag, data, elements, boundingSphere ? () => boundingSphere : undefined, () => stepIndices[0] !== undefined ? NtCTubeSegmentLabel(data[stepIndices[0]]) : '');
    }

    export function DummyLoci(): DummyLoci {
        return DataLoci(DummyTag, {}, [], undefined, () => '');
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'data-loci' && x.tag === DataTag;
    }
}
