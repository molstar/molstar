/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { parseVtp } from '../parser';

// Builds a minimal uncompressed raw-appended VTP in memory.
// Layout: positions at offset 0, connectivity at offset 40, offsets at offset 56.
// Each block: [UInt32 nbytes LE][raw data bytes].
function makeTriangleVtp(): Uint8Array {
    const header =
        '<?xml version="1.0"?>\n' +
        '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian" header_type="UInt32">\n' +
        '  <PolyData>\n' +
        '    <Piece NumberOfPoints="3" NumberOfPolys="1">\n' +
        '      <Points>\n' +
        '        <DataArray type="Float32" Name="Points" NumberOfComponents="3" format="appended" offset="0"/>\n' +
        '      </Points>\n' +
        '      <Polys>\n' +
        '        <DataArray type="Int32" Name="connectivity" format="appended" offset="40"/>\n' +
        '        <DataArray type="Int32" Name="offsets" format="appended" offset="56"/>\n' +
        '      </Polys>\n' +
        '    </Piece>\n' +
        '  </PolyData>\n' +
        '  <AppendedData encoding="raw">\n_';
    const footer = '\n  </AppendedData>\n</VTKFile>\n';

    // Positions: (0,0,0) (1,0,0) (0,1,0) — 9 Float32 = 36 bytes
    const pos = new Float32Array([0, 0, 0, 1, 0, 0, 0, 1, 0]);
    // Connectivity: triangle [0,1,2] — 3 Int32 = 12 bytes; offset 40 = 4+36
    const conn = new Int32Array([0, 1, 2]);
    // Offsets: [3] — 1 Int32 = 4 bytes; offset 56 = 40+4+12
    const offs = new Int32Array([3]);

    return buildRawAppended(header, footer, [pos, conn, offs]);
}

// Builds a VTP with a quad face to exercise fan-triangulation (quad → 2 triangles).
function makeQuadVtp(): Uint8Array {
    const header =
        '<?xml version="1.0"?>\n' +
        '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian" header_type="UInt32">\n' +
        '  <PolyData>\n' +
        '    <Piece NumberOfPoints="4" NumberOfPolys="1">\n' +
        '      <Points>\n' +
        '        <DataArray type="Float32" Name="Points" NumberOfComponents="3" format="appended" offset="0"/>\n' +
        '      </Points>\n' +
        '      <Polys>\n' +
        '        <DataArray type="Int32" Name="connectivity" format="appended" offset="52"/>\n' +
        '        <DataArray type="Int32" Name="offsets" format="appended" offset="72"/>\n' +
        '      </Polys>\n' +
        '    </Piece>\n' +
        '  </PolyData>\n' +
        '  <AppendedData encoding="raw">\n_';
    const footer = '\n  </AppendedData>\n</VTKFile>\n';

    // 4 vertices: (0,0,0) (1,0,0) (1,1,0) (0,1,0) — 12 Float32 = 48 bytes; block = 4+48 = 52
    const pos = new Float32Array([0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0]);
    // Connectivity: quad [0,1,2,3] — 4 Int32 = 16 bytes; offset 52, block = 4+16 = 20
    const conn = new Int32Array([0, 1, 2, 3]);
    // Offsets: [4] — 1 Int32 = 4 bytes; offset 72
    const offs = new Int32Array([4]);

    return buildRawAppended(header, footer, [pos, conn, offs]);
}

// Builds a VTP with one triangle and one CellData Float32 scalar.
function makeCellDataVtp(): Uint8Array {
    // Layout: positions offset=0 (block=40), connectivity offset=40 (block=16),
    //         offsets offset=56 (block=8), curvature offset=64 (block=8)
    const header =
        '<?xml version="1.0"?>\n' +
        '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian" header_type="UInt32">\n' +
        '  <PolyData>\n' +
        '    <Piece NumberOfPoints="3" NumberOfPolys="1">\n' +
        '      <CellData>\n' +
        '        <DataArray type="Float32" Name="curvature" NumberOfComponents="1" format="appended" offset="64" RangeMin="-1.0" RangeMax="1.0"/>\n' +
        '      </CellData>\n' +
        '      <Points>\n' +
        '        <DataArray type="Float32" Name="Points" NumberOfComponents="3" format="appended" offset="0"/>\n' +
        '      </Points>\n' +
        '      <Polys>\n' +
        '        <DataArray type="Int32" Name="connectivity" format="appended" offset="40"/>\n' +
        '        <DataArray type="Int32" Name="offsets" format="appended" offset="56"/>\n' +
        '      </Polys>\n' +
        '    </Piece>\n' +
        '  </PolyData>\n' +
        '  <AppendedData encoding="raw">\n_';
    const footer = '\n  </AppendedData>\n</VTKFile>\n';

    const pos = new Float32Array([0, 0, 0, 1, 0, 0, 0, 1, 0]);
    const conn = new Int32Array([0, 1, 2]);
    const offs = new Int32Array([3]);
    const curv = new Float32Array([0.5]);

    return buildRawAppended(header, footer, [pos, conn, offs, curv]);
}

function buildRawAppended(header: string, footer: string, arrays: (Float32Array | Int32Array)[]): Uint8Array {
    const blocks = arrays.map(arr => {
        const data = new Uint8Array(arr.buffer, arr.byteOffset, arr.byteLength);
        const block = new Uint8Array(4 + data.length);
        new DataView(block.buffer).setUint32(0, data.length, true);
        block.set(data, 4);
        return block;
    });

    const enc = new TextEncoder();
    const headerBytes = enc.encode(header);
    const footerBytes = enc.encode(footer);
    const totalBinary = blocks.reduce((s, b) => s + b.length, 0);

    const out = new Uint8Array(headerBytes.length + totalBinary + footerBytes.length);
    let pos = 0;
    out.set(headerBytes, pos); pos += headerBytes.length;
    for (const block of blocks) { out.set(block, pos); pos += block.length; }
    out.set(footerBytes, pos);
    return out;
}

describe('VTP parser', () => {
    it('parses a minimal triangle', async () => {
        const result = await parseVtp(makeTriangleVtp()).run();
        expect(result.isError).toBe(false);
        if (result.isError) return;

        const vtp = result.result;
        expect(vtp.numberOfPoints).toBe(3);
        expect(vtp.numberOfCells).toBe(1);
        expect(vtp.numberOfTriangles).toBe(1);
        expect(vtp.positions.length).toBe(9);
        expect(vtp.connectivity.length).toBe(3);

        expect(vtp.positions[0]).toBeCloseTo(0);
        expect(vtp.positions[3]).toBeCloseTo(1);
        expect(vtp.positions[7]).toBeCloseTo(1);

        expect(vtp.pointData.size).toBe(0);
        expect(vtp.cellData.size).toBe(0);

        expect.assertions(11);
    });

    it('fan-triangulates a quad into 2 triangles', async () => {
        const result = await parseVtp(makeQuadVtp()).run();
        expect(result.isError).toBe(false);
        if (result.isError) return;

        const vtp = result.result;
        expect(vtp.numberOfPoints).toBe(4);
        expect(vtp.numberOfCells).toBe(1);
        expect(vtp.numberOfTriangles).toBe(2);
        expect(vtp.connectivity.length).toBe(6);

        // Fan from v0: tri0 = (0,1,2), tri1 = (0,2,3)
        expect(vtp.connectivity[0]).toBe(0);
        expect(vtp.connectivity[1]).toBe(1);
        expect(vtp.connectivity[2]).toBe(2);
        expect(vtp.connectivity[3]).toBe(0);
        expect(vtp.connectivity[4]).toBe(2);
        expect(vtp.connectivity[5]).toBe(3);

        expect.assertions(11);
    });

    it('parses CellData scalar array', async () => {
        const result = await parseVtp(makeCellDataVtp()).run();
        expect(result.isError).toBe(false);
        if (result.isError) return;

        const vtp = result.result;
        expect(vtp.cellData.size).toBe(1);
        expect(vtp.pointData.size).toBe(0);

        const curv = vtp.cellData.get('curvature');
        expect(curv).toBeDefined();
        if (!curv) return;

        expect(curv.values.length).toBe(1);
        expect(curv.values[0]).toBeCloseTo(0.5);
        expect(curv.desc.rangeMin).toBeCloseTo(-1.0);
        expect(curv.desc.rangeMax).toBeCloseTo(1.0);

        expect.assertions(8);
    });
});
