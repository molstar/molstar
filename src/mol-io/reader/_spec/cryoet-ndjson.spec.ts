/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author OpenAI
 */

import { parseCryoEtDataPortalNdjson } from '../cryoet/ndjson';

test('parses CryoET Data Portal ndjson records', async () => {
    const data = [
        JSON.stringify({ type: 'point', location: { x: 1, y: 2, z: 3 } }),
        JSON.stringify({ type: 'orientedPoint', location: { x: 4, y: 5, z: 6 }, xyz_rotation_matrix: [[1, 0, 0], [0, 0, -1], [0, 1, 0]], instance_id: 7 }),
    ].join('\n');

    const parsed = await parseCryoEtDataPortalNdjson(data).run();
    if (parsed.isError) throw new Error(parsed.message);

    expect(parsed.result.records).toHaveLength(2);
    expect(parsed.result.records[0].type).toBe('point');
    expect(parsed.result.records[1].type).toBe('orientedPoint');
    expect(parsed.result.records[1].instance_id).toBe(7);
    expect(parsed.result.records[1].location).toEqual({ x: 4, y: 5, z: 6 });
});
