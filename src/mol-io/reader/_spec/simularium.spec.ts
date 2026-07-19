/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { parseSimularium } from '../simularium/parser';
import { createParticleListFromSimularium, getSimulariumAgentTypeNames, getSimulariumFrameCount } from '../../../mol-model-formats/particles/simularium';

const trajectoryInfo = {
    version: 3,
    totalSteps: 1,
    timeStepSize: 1,
    spatialUnits: { magnitude: 1, name: 'nm' },
    size: { x: 100, y: 100, z: 100 },
    typeMapping: {
        '0': { name: 'A', geometry: { displayType: 'SPHERE' } },
        '1': { name: 'B', geometry: { displayType: 'FIBER' } },
    },
};

// agent 0: DEFAULT sphere at (1,2,3), no rotation, radius 5, 0 subpoints (11 floats)
// agent 1: FIBER at origin, no rotation, radius 1, 6 subpoint floats = 2 points (17 floats)
const agents = [
    1000, 0, 0, 1, 2, 3, 0, 0, 0, 5, 0,
    1001, 1, 1, 0, 0, 0, 0, 0, 0, 1, 6, 10, 0, 0, 20, 0, 0,
];

function buildJson(): Uint8Array {
    const json = {
        trajectoryInfo,
        spatialData: {
            version: 1,
            msgType: 1,
            bundleStart: 0,
            bundleSize: 1,
            bundleData: [{ frameNumber: 0, time: 0, data: agents }],
        },
        plotData: { version: 1, data: [] },
    };
    return new TextEncoder().encode(JSON.stringify(json));
}

const SIGNATURE = 'SIMULARIUMBINARY';

function buildBinary(): Uint8Array {
    const trajJson = new TextEncoder().encode(JSON.stringify(trajectoryInfo));
    const trajPadded = Math.ceil(trajJson.length / 4) * 4;

    const nBlocks = 2;
    const headerLength = SIGNATURE.length + 12 + nBlocks * 12; // 16 + 12 + 24 = 52

    const block0Offset = headerLength;
    const block0Size = 8 + trajPadded;

    const frameLength = 12 + agents.length * 4; // frameNumber + time + nAgents + agents
    const frameOffsetInBlock = 8 + 8 + 8; // block header + (version + nFrames) + 1-frame toc
    const block1Content = 8 + 8 + frameLength; // version + nFrames + toc + frame
    const block1Offset = block0Offset + block0Size;
    const block1Size = 8 + block1Content;

    const buf = new Uint8Array(block1Offset + block1Size);
    const view = new DataView(buf.buffer);

    for (let i = 0; i < SIGNATURE.length; ++i) buf[i] = SIGNATURE.charCodeAt(i);

    let p = SIGNATURE.length;
    view.setUint32(p, headerLength, true); p += 4;
    view.setUint32(p, 1, true); p += 4; // header version
    view.setUint32(p, nBlocks, true); p += 4;
    view.setUint32(p, block0Offset, true); p += 4;
    view.setUint32(p, 1, true); p += 4; // block type: trajectory info json
    view.setUint32(p, block0Size, true); p += 4;
    view.setUint32(p, block1Offset, true); p += 4;
    view.setUint32(p, 3, true); p += 4; // block type: spatial data binary
    view.setUint32(p, block1Size, true); p += 4;

    view.setUint32(block0Offset, 1, true);
    view.setUint32(block0Offset + 4, block0Size, true);
    buf.set(trajJson, block0Offset + 8);

    view.setUint32(block1Offset, 3, true);
    view.setUint32(block1Offset + 4, block1Size, true);
    let q = block1Offset + 8;
    view.setUint32(q, 1, true); q += 4; // spatial data version
    view.setUint32(q, 1, true); q += 4; // nFrames
    view.setUint32(q, frameOffsetInBlock, true); q += 4; // frame offset (relative to block start)
    view.setUint32(q, frameLength, true); q += 4; // frame length
    view.setUint32(q, 0, true); q += 4; // frameNumber
    view.setFloat32(q, 0, true); q += 4; // time
    view.setUint32(q, 2, true); q += 4; // nAgents
    for (let i = 0; i < agents.length; ++i) { view.setFloat32(q, agents[i], true); q += 4; }

    return buf;
}

test('parses JSON simularium', async () => {
    const parsed = await parseSimularium(buildJson()).run();
    if (parsed.isError) throw new Error(parsed.message);

    const file = parsed.result;
    expect(getSimulariumFrameCount(file)).toBe(1);
    expect(file.trajectoryInfo.spatialUnits!.name).toBe('nm');
    expect(file.frames[0].data.length).toBe(agents.length);
    expect(Array.from(file.frames[0].data)).toEqual(agents);
});

test('parses binary simularium', async () => {
    const parsed = await parseSimularium(buildBinary()).run();
    if (parsed.isError) throw new Error(parsed.message);

    const file = parsed.result;
    expect(getSimulariumFrameCount(file)).toBe(1);
    expect(file.trajectoryInfo.typeMapping['1'].name).toBe('B');
    expect(file.frames[0].data.length).toBe(agents.length);
    expect(Array.from(file.frames[0].data)).toEqual(agents);
});

test('binary and JSON frames are identical', async () => {
    const fromJson = await parseSimularium(buildJson()).run();
    const fromBinary = await parseSimularium(buildBinary()).run();
    if (fromJson.isError) throw new Error(fromJson.message);
    if (fromBinary.isError) throw new Error(fromBinary.message);
    expect(Array.from(fromBinary.result.frames[0].data)).toEqual(Array.from(fromJson.result.frames[0].data));
});

test('builds particle list with exploded fibers', async () => {
    const parsed = await parseSimularium(buildJson()).run();
    if (parsed.isError) throw new Error(parsed.message);

    const list = createParticleListFromSimularium(parsed.result, { frameIndex: 0, scale: 1 });

    // 1 default sphere + 2 fiber subpoint spheres
    expect(list.count).toBe(3);

    // default agent position
    expect(list.coordinates[0]).toBeCloseTo(1);
    expect(list.coordinates[1]).toBeCloseTo(2);
    expect(list.coordinates[2]).toBeCloseTo(3);

    // fiber subpoint positions (identity rotation, origin agent)
    expect(list.coordinates[3]).toBeCloseTo(10);
    expect(list.coordinates[6]).toBeCloseTo(20);

    // fiber connectivity
    expect(list.fibers).toBeDefined();
    expect(list.fibers!.count).toBe(1);
    expect(Array.from(list.fibers!.offsets)).toEqual([0, 2]);
    expect(Array.from(list.fibers!.indices)).toEqual([1, 2]);

    // two distinct entities (A and B)
    expect(list.entityInfo!.size).toBe(2);
});

test('applies spatial unit auto-scaling (nm to angstrom)', async () => {
    const parsed = await parseSimularium(buildJson()).run();
    if (parsed.isError) throw new Error(parsed.message);

    const list = createParticleListFromSimularium(parsed.result, { frameIndex: 0 });
    // nm -> angstrom factor of 10
    expect(list.coordinates[0]).toBeCloseTo(10);
    expect(list.coordinates[1]).toBeCloseTo(20);
    expect(list.coordinates[2]).toBeCloseTo(30);
});

test('lists agent type names', async () => {
    const parsed = await parseSimularium(buildJson()).run();
    if (parsed.isError) throw new Error(parsed.message);

    const types = getSimulariumAgentTypeNames(parsed.result);
    expect(types).toEqual([{ id: 0, name: 'A' }, { id: 1, name: 'B' }]);
});
