/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { deflate, inflate, unzip, zip } from '../zip/zip';
import { SyncRuntimeContext } from '../../mol-task/execution/synchronous';
import { describe, it, expect } from 'vitest';

describe('zip', () => {
    it('roundtrip deflate/inflate', async () => {
        const data = new Uint8Array([1, 2, 3, 4, 5, 6, 7]);
        const deflated = await deflate(SyncRuntimeContext, data);
        const inflated = await inflate(SyncRuntimeContext, deflated);
        expect(inflated).toEqual(data);
    });

    it('roundtrip zip/unzip', async () => {
        const data = {
            'test.foo': new Uint8Array([1, 2, 3, 4, 5, 6, 7])
        };
        const zipped = await zip(SyncRuntimeContext, data);
        const unzipped = await unzip(SyncRuntimeContext, zipped);
        expect(unzipped).toEqual(data);
    });
});