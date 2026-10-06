/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { PNG } from 'pngjs';
import { getAreaTexture, getSearchTexture } from '../../../mol-canvas3d/passes/smaa-textures';
import { getSmaaLookupData } from '../smaa-lookups';

it('native SMAA lookup bytes exactly match both shared PNG tables without DOM decoding', async () => {
    const data = await getSmaaLookupData();
    for (const [image, url] of [[data.area, getAreaTexture()], [data.search, getSearchTexture()]] as const) {
        const decoded = PNG.sync.read(Buffer.from(url.split(',')[1], 'base64'));
        expect([image.width, image.height]).toEqual([decoded.width, decoded.height]);
        expect(Buffer.from(image.array).equals(decoded.data)).toBe(true);
    }
});
