/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement } from '../../element';
import { IntraUnitBonds } from './data';
import { SortedArray } from '../../../../../mol-data/int';
import { Model } from '../../../model';

export class ElementSetIntraBondCache {
    private data = new Map<number, [StructureElement.Set, IntraUnitBonds][]>();

    get(xs: StructureElement.Set): IntraUnitBonds | undefined {
        const hash = SortedArray.hashCode(xs);
        if (!this.data.has(hash)) return void 0;
        for (const [s, b] of this.data.get(hash)!) {
            if (SortedArray.areEqual(xs, s)) return b;
        }
    }

    set(xs: StructureElement.Set, bonds: IntraUnitBonds) {
        const hash = SortedArray.hashCode(xs);
        if (this.data.has(hash)) {
            const es = this.data.get(hash)!;
            for (const e of es) {
                if (SortedArray.areEqual(xs, e[0])) {
                    e[1] = bonds;
                    return;
                }
            }
            es.push([xs, bonds]);
        } else {
            this.data.set(hash, [[xs, bonds]]);
        }
    }

    static get(model: Model): ElementSetIntraBondCache {
        if (!model._dynamicPropertyData.ElementSetIntraBondCache) {
            model._dynamicPropertyData.ElementSetIntraBondCache = new ElementSetIntraBondCache();
        }
        return model._dynamicPropertyData.ElementSetIntraBondCache;
    }
}