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

// Builds a VTP with 2 quads and per-cell CellData to exercise triangleCellIndex
// for non-triangle polygons.  Each quad fans into 2 triangles, so 4 triangles total.
// Cell 0 → triangles 0,1; cell 1 → triangles 2,3.
function makeQuadCellDataVtp(): Uint8Array {
    // 6 unique points laid out as two quads sharing an edge:
    //   quad0: v0(0,0,0) v1(1,0,0) v2(1,1,0) v3(0,1,0)
    //   quad1: v4(2,0,0) v5(2,1,0); reuses v1,v2
    // connectivity: [0,1,2,3] [1,4,5,2]  offsets: [4, 8]
    // Layout: positions offset=0 (block=4+72=76), connectivity offset=76 (block=4+32=36),
    //         offsets offset=112 (block=4+8=12), cellval offset=124 (block=4+8=12)
    const header =
        '<?xml version="1.0"?>\n' +
        '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian" header_type="UInt32">\n' +
        '  <PolyData>\n' +
        '    <Piece NumberOfPoints="6" NumberOfPolys="2">\n' +
        '      <CellData>\n' +
        '        <DataArray type="Float32" Name="cellval" NumberOfComponents="1" format="appended" offset="124" RangeMin="0.0" RangeMax="1.0"/>\n' +
        '      </CellData>\n' +
        '      <Points>\n' +
        '        <DataArray type="Float32" Name="Points" NumberOfComponents="3" format="appended" offset="0"/>\n' +
        '      </Points>\n' +
        '      <Polys>\n' +
        '        <DataArray type="Int32" Name="connectivity" format="appended" offset="76"/>\n' +
        '        <DataArray type="Int32" Name="offsets" format="appended" offset="112"/>\n' +
        '      </Polys>\n' +
        '    </Piece>\n' +
        '  </PolyData>\n' +
        '  <AppendedData encoding="raw">\n_';
    const footer = '\n  </AppendedData>\n</VTKFile>\n';

    // 6 points × 3 floats = 18 Float32 = 72 bytes; block = 4+72 = 76
    const pos = new Float32Array([0,0,0, 1,0,0, 1,1,0, 0,1,0, 2,0,0, 2,1,0]);
    // connectivity: 8 Int32 = 32 bytes; block = 4+32 = 36; offset=76
    const conn = new Int32Array([0,1,2,3, 1,4,5,2]);
    // offsets: 2 Int32 = 8 bytes; block = 4+8 = 12; offset=112
    const offs = new Int32Array([4, 8]);
    // cellval: 2 Float32 = 8 bytes; block = 4+8 = 12; offset=124
    const cellval = new Float32Array([0.25, 0.75]);

    return buildRawAppended(header, footer, [pos, conn, offs, cellval]);
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

function makeAsciiTriangleVtp(): Uint8Array {
    const xml = [
        '<?xml version="1.0"?>',
        '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">',
        '  <PolyData>',
        '    <Piece NumberOfPoints="3" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="1">',
        '      <Points>',
        '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">',
        '          0.0 0.0 0.0  1.0 0.0 0.0  0.0 1.0 0.0',
        '        </DataArray>',
        '      </Points>',
        '      <CellData Scalars="val">',
        '        <DataArray type="Float32" Name="val" NumberOfComponents="1" format="ascii" RangeMin="0.5" RangeMax="0.5">',
        '          0.5',
        '        </DataArray>',
        '      </CellData>',
        '      <Polys>',
        '        <DataArray type="Int32" Name="connectivity" format="ascii">',
        '          0 1 2',
        '        </DataArray>',
        '        <DataArray type="Int32" Name="offsets" format="ascii">',
        '          3',
        '        </DataArray>',
        '      </Polys>',
        '    </Piece>',
        '  </PolyData>',
        '</VTKFile>',
    ].join('\n');
    return new TextEncoder().encode(xml);
}

function makeAsciiMultiComponentVtp(): Uint8Array {
    // Single triangle with a 3-component PointData vector (normals)
    const xml = [
        '<?xml version="1.0"?>',
        '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">',
        '  <PolyData>',
        '    <Piece NumberOfPoints="3" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="1">',
        '      <Points>',
        '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">',
        '          0.0 0.0 0.0  1.0 0.0 0.0  0.0 1.0 0.0',
        '        </DataArray>',
        '      </Points>',
        '      <PointData>',
        '        <DataArray type="Float32" Name="normals" NumberOfComponents="3" format="ascii">',
        '          0.0 0.0 1.0  0.0 0.0 1.0  0.0 0.0 1.0',
        '        </DataArray>',
        '      </PointData>',
        '      <Polys>',
        '        <DataArray type="Int32" Name="connectivity" format="ascii">',
        '          0 1 2',
        '        </DataArray>',
        '        <DataArray type="Int32" Name="offsets" format="ascii">',
        '          3',
        '        </DataArray>',
        '      </Polys>',
        '    </Piece>',
        '  </PolyData>',
        '</VTKFile>',
    ].join('\n');
    return new TextEncoder().encode(xml);
}

describe('VTP parser', () => {
    it('parses a minimal triangle', async () => {
        expect.assertions(11);
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
    });

    it('fan-triangulates a quad into 2 triangles', async () => {
        expect.assertions(11);
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
    });

    it('parses CellData scalar array', async () => {
        expect.assertions(8);
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
    });

    it('triangleCellIndex maps quad fan-triangles to their originating cells', async () => {
        expect.assertions(10);
        const result = await parseVtp(makeQuadCellDataVtp()).run();
        expect(result.isError).toBe(false);
        if (result.isError) return;

        const vtp = result.result;
        // 2 quads × 2 triangles each = 4 triangles
        expect(vtp.numberOfTriangles).toBe(4);
        expect(vtp.triangleCellIndex.length).toBe(4);

        // Triangles 0 and 1 originate from cell 0
        expect(vtp.triangleCellIndex[0]).toBe(0);
        expect(vtp.triangleCellIndex[1]).toBe(0);
        // Triangles 2 and 3 originate from cell 1
        expect(vtp.triangleCellIndex[2]).toBe(1);
        expect(vtp.triangleCellIndex[3]).toBe(1);

        // CellData must have 2 entries with correct values
        const cellval = vtp.cellData.get('cellval');
        expect(cellval).toBeDefined();
        if (!cellval) return;
        expect(cellval.values[0]).toBeCloseTo(0.25);
        expect(cellval.values[1]).toBeCloseTo(0.75);
    });

    it('parses an ASCII-format VTP file', async () => {
        expect.assertions(9);
        const result = await parseVtp(makeAsciiTriangleVtp()).run();
        expect(result.isError).toBe(false);
        if (result.isError) return;

        const vtp = result.result;
        expect(vtp.numberOfPoints).toBe(3);
        expect(vtp.numberOfCells).toBe(1);
        expect(vtp.numberOfTriangles).toBe(1);
        expect(vtp.positions[0]).toBeCloseTo(0);
        expect(vtp.positions[3]).toBeCloseTo(1);

        const val = vtp.cellData.get('val');
        expect(val).toBeDefined();
        if (!val) return;
        expect(val.values[0]).toBeCloseTo(0.5);
        expect(val.desc.format).toBe('ascii');
    });

    it('stores multi-component PointData arrays', async () => {
        expect.assertions(8);
        const result = await parseVtp(makeAsciiMultiComponentVtp()).run();
        expect(result.isError).toBe(false);
        if (result.isError) return;

        const vtp = result.result;
        const normals = vtp.pointData.get('normals');
        expect(normals).toBeDefined();
        if (!normals) return;

        expect(normals.desc.numberOfComponents).toBe(3);
        // 3 points × 3 components = 9 values
        expect(normals.values.length).toBe(9);
        // Each normal is (0,0,1) — third component of first vertex
        expect(normals.values[0]).toBeCloseTo(0);
        expect(normals.values[1]).toBeCloseTo(0);
        expect(normals.values[2]).toBeCloseTo(1);
        expect(normals.desc.format).toBe('ascii');
    });
});
