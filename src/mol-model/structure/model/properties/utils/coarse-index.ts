/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ElementIndex } from '../../indexing';
import { CoarseElementData, CoarseElementReference, CoarseIndex, CoarseElementKey } from '../coarse';

export function getCoarseIndex(data: { spheres: CoarseElementData, gaussians: CoarseElementData }): CoarseIndex {
    return new Index(data);
}

class EmptyIndex implements CoarseIndex {
    findElement(key: CoarseElementKey, out: CoarseElementReference): boolean {
        out.kind = undefined;
        out.index = -1 as ElementIndex;
        return false;
    }
    findSphereElement(key: CoarseElementKey): ElementIndex {
        return -1 as ElementIndex;
    }
    findGaussianElement(key: CoarseElementKey): ElementIndex {
        return -1 as ElementIndex;
    }
}

export const EmptyCoarseIndex: CoarseIndex = new EmptyIndex();

class Index implements CoarseIndex {
    private _sphereMapping: CoarseElementMapping | undefined = void 0;
    private _gaussianMapping: CoarseElementMapping | undefined = void 0;

    get sphereMapping() {
        if (!this._sphereMapping) this._sphereMapping = buildMapping(this.data.spheres);
        return this._sphereMapping;
    }

    get gaussianMapping() {
        if (!this._gaussianMapping) this._gaussianMapping = buildMapping(this.data.gaussians);
        return this._gaussianMapping;
    }

    findSphereElement(key: CoarseElementKey): ElementIndex {
        const mapping = this.sphereMapping;
        let xs: any = mapping[key.label_entity_id];
        if (!xs) return -1 as ElementIndex;
        xs = xs[key.label_asym_id];
        if (!xs) return -1 as ElementIndex;
        return xs[key.label_seq_id] ?? -1;
    }

    findGaussianElement(key: CoarseElementKey): ElementIndex {
        const mapping = this.gaussianMapping;
        let xs: any = mapping[key.label_entity_id];
        if (!xs) return -1 as ElementIndex;
        xs = xs[key.label_asym_id];
        if (!xs) return -1 as ElementIndex;
        return xs[key.label_seq_id] ?? -1;
    }

    findElement(key: CoarseElementKey, out: CoarseElementReference): boolean {
        const sphere = this.findSphereElement(key);
        if (sphere >= 0) {
            out.kind = 'spheres';
            out.index = sphere;
            return true;
        }
        const gaussian = this.findGaussianElement(key);
        if (gaussian >= 0) {
            out.kind = 'gaussians';
            out.index = gaussian;
            return true;
        }
        return false;
    }

    constructor(private data: { spheres: CoarseElementData, gaussians: CoarseElementData }) {
    }
}

type CoarseElementMapping = { [entityId: string]: { [chainId: string]: { [seqId: number]: ElementIndex } } };

function buildMapping({ count, entity_id, asym_id, seq_id_begin, seq_id_end }: CoarseElementData): CoarseElementMapping {
    const ret: CoarseElementMapping = {};
    for (let i = 0; i < count; i++) {
        const entityId = entity_id.value(i);
        const asymId = asym_id.value(i);

        if (!ret[entityId]) ret[entityId] = {};
        if (!ret[entityId][asymId]) ret[entityId][asymId] = {};

        const elements = ret[entityId][asymId];
        const seqIdBegin = seq_id_begin.value(i);
        const seqIdEnd = seq_id_end.value(i);
        for (let seqId = seqIdBegin; seqId <= seqIdEnd; seqId++) {
            elements[seqId] = i as ElementIndex;
        }
    }

    return ret;
}