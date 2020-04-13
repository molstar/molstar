/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createReferenceCache } from '../reference-cache';

describe('reference-cache', () => {
    it('basic', () => {
        const refCache = createReferenceCache(
            (x: number) => x.toString(),
            (x) => x,
            () => {}
        );
        expect(refCache.count).toBe(0);

        const ref2a = refCache.get(2);
        expect(refCache.count).toBe(1);

        const ref2b = refCache.get(2);
        expect(refCache.count).toBe(1);
        expect(ref2b.value).toBe(2);

        const ref3 = refCache.get(3);
        expect(refCache.count).toBe(2);
        expect(ref3.value).toBe(3);

        ref2a.free();
        refCache.clear();
        expect(refCache.count).toBe(2);
        ref2b.free();
        refCache.clear();
        expect(refCache.count).toBe(1);

        refCache.dispose();
        expect(refCache.count).toBe(0);
    });
});