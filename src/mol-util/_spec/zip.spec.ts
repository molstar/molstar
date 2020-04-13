/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { deflate, inflate, unzip, zip } from '../zip/zip';
import { SyncRuntimeContext } from '../../mol-task/execution/synchronous';

describe('zip', () => {
    it('roundtrip deflate/inflate', async () => {
        const data = new Uint8Array([1, 2, 3, 4, 5, 6, 7]);
        const deflated = deflate(data);
        const inflated = await inflate(SyncRuntimeContext, deflated);
        expect(inflated).toEqual(data);
    });

    it('roundtrip zip', async () => {
        const data = {
            'test.foo': new Uint8Array([1, 2, 3, 4, 5, 6, 7])
        };
        const zipped = zip(data);
        const unzipped = await unzip(SyncRuntimeContext, zipped);
        expect(unzipped).toEqual(data);
    });
});