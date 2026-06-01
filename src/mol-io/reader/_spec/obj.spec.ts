/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { parseObj } from '../obj/parser';

// Simple triangle
const objTriangle = `# simple triangle
v 0.0 0.0 0.0
v 1.0 0.0 0.0
v 0.5 1.0 0.0
f 1 2 3
`;

// Quad that gets fan-triangulated into 2 triangles
const objQuad = `# quad fan-triangulated
v -1.0 -1.0 0.0
v  1.0 -1.0 0.0
v  1.0  1.0 0.0
v -1.0  1.0 0.0
f 1 2 3 4
`;

// Vertex normals
const objWithNormals = `# vertex normals
v 0.0 0.0 0.0
v 1.0 0.0 0.0
v 0.5 1.0 0.0
vn 0.0 0.0 1.0
vn 0.0 0.0 1.0
vn 0.0 0.0 1.0
f 1//1 2//2 3//3
`;

// v/vt/vn format (texture coords are ignored but should not break parsing)
const objWithTexture = `v 0.0 0.0 0.0
v 1.0 0.0 0.0
v 0.5 1.0 0.0
vn 0.0 0.0 1.0
f 1/1/1 2/2/1 3/3/1
`;

// Multiple materials / usemtl groups
const objMultiMaterial = `v 0.0 0.0 0.0
v 1.0 0.0 0.0
v 0.5 1.0 0.0
v 2.0 0.0 0.0
v 2.5 1.0 0.0
usemtl red
f 1 2 3
usemtl green
f 2 4 5
`;

// Negative indices (relative addressing)
const objNegativeIndices = `v 0.0 0.0 0.0
v 1.0 0.0 0.0
v 0.5 1.0 0.0
f -3 -2 -1
`;

// Comments and blank lines should be ignored
const objWithComments = `# header comment
# another comment

v 0.0 0.0 0.0
# inline comment after data

v 1.0 0.0 0.0
v 0.5 1.0 0.0

f 1 2 3
`;

// Unsupported directives (g, o, s, mtllib, vt, vp) should be silently skipped
const objUnsupportedDirectives = `mtllib material.mtl
o MyObject
g mygroup
s 1
v 0.0 0.0 0.0
vt 0.0 0.0
vp 0.0 1.0
v 1.0 0.0 0.0
v 0.5 1.0 0.0
f 1 2 3
`;

// Cube (6 faces × 2 triangles = 12 triangles)
const objCube = `# unit cube
v -1 -1 -1
v  1 -1 -1
v  1  1 -1
v -1  1 -1
v -1 -1  1
v  1 -1  1
v  1  1  1
v -1  1  1
# bottom (-z)
f 1 2 3
f 1 3 4
# top (+z)
f 5 6 7
f 5 7 8
# front (+x)
f 2 6 7
f 2 7 3
# back (-x)
f 5 1 4
f 5 4 8
# left (-y)
f 1 5 6
f 1 6 2
# right (+y)
f 4 3 7
f 4 7 8
`;

// CRLF line endings
const objCRLF = '# crlf triangle\r\nv 0.0 0.0 0.0\r\nv 1.0 0.0 0.0\r\nv 0.5 1.0 0.0\r\nf 1 2 3\r\n';

// Tabs and leading whitespace before keywords
const objLeadingWhitespace = '\tv 0.0 0.0 0.0\n  v 1.0 0.0 0.0\n\t v 0.5 1.0 0.0\n\tf 1 2 3\n';

// Degenerate face (fewer than 3 vertices) should be skipped with a warning
const objDegenerateFace = `v 0.0 0.0 0.0
v 1.0 0.0 0.0
v 0.5 1.0 0.0
f 1 2
f 1 2 3
`;

// Mixed face-vertices: some reference a normal, some do not, within one mesh
const objMixedNormals = `v 0.0 0.0 0.0
v 1.0 0.0 0.0
v 0.5 1.0 0.0
vn 0.0 0.0 1.0
f 1//1 2 3//1
`;

// Negative normal indices (relative addressing for normals)
const objNegativeNormalIndices = `v 0.0 0.0 0.0
v 1.0 0.0 0.0
v 0.5 1.0 0.0
vn 0.0 0.0 1.0
f 1//-1 2//-1 3//-1
`;

// usemtl reuse: a material name is referenced again after another material
const objReusedMaterial = `v 0.0 0.0 0.0
v 1.0 0.0 0.0
v 0.5 1.0 0.0
v 2.0 0.0 0.0
v 2.5 1.0 0.0
usemtl red
f 1 2 3
usemtl green
f 2 4 5
usemtl red
f 1 3 4
`;

// Empty file
const objEmpty = '';

// File with vertices but no faces
const objNoFaces = `v 0.0 0.0 0.0
v 1.0 0.0 0.0
v 0.5 1.0 0.0
`;

describe('obj reader', () => {
    it('parses a simple triangle', async () => {
        const parsed = await parseObj(objTriangle).run();
        if (parsed.isError) throw new Error(parsed.message);
        const obj = parsed.result;

        expect(obj.positionCount).toBe(3);
        expect(obj.triangleCount).toBe(1);
        // First vertex
        expect(obj.positions[0]).toBeCloseTo(0.0);
        expect(obj.positions[1]).toBeCloseTo(0.0);
        expect(obj.positions[2]).toBeCloseTo(0.0);
        // Triangle indices (0-based)
        expect(Array.from(obj.positionIndices)).toEqual([0, 1, 2]);
    });

    it('fan-triangulates a quad into two triangles', async () => {
        const parsed = await parseObj(objQuad).run();
        if (parsed.isError) throw new Error(parsed.message);
        const obj = parsed.result;

        expect(obj.positionCount).toBe(4);
        expect(obj.triangleCount).toBe(2);
        // Fan from vertex 0: (0,1,2) and (0,2,3)
        expect(Array.from(obj.positionIndices)).toEqual([0, 1, 2, 0, 2, 3]);
    });

    it('parses vertex normals with v//vn format', async () => {
        const parsed = await parseObj(objWithNormals).run();
        if (parsed.isError) throw new Error(parsed.message);
        const obj = parsed.result;

        expect(obj.normalCount).toBe(3);
        expect(obj.triangleCount).toBe(1);
        expect(Array.from(obj.normalIndices)).toEqual([0, 1, 2]);
        // Normal z component of first normal
        expect(obj.normals[2]).toBeCloseTo(1.0);
    });

    it('parses v/vt/vn format (texture coords ignored)', async () => {
        const parsed = await parseObj(objWithTexture).run();
        if (parsed.isError) throw new Error(parsed.message);
        const obj = parsed.result;

        expect(obj.positionCount).toBe(3);
        expect(obj.normalCount).toBe(1);
        expect(obj.triangleCount).toBe(1);
        expect(Array.from(obj.normalIndices)).toEqual([0, 0, 0]);
    });

    it('assigns material groups via usemtl', async () => {
        const parsed = await parseObj(objMultiMaterial).run();
        if (parsed.isError) throw new Error(parsed.message);
        const obj = parsed.result;

        expect(obj.triangleCount).toBe(2);
        // 'default' is always first; 'red' and 'green' added on use
        expect(obj.groups[1]).toBe('red');
        expect(obj.groups[2]).toBe('green');
        // First triangle belongs to 'red' (index 1), second to 'green' (index 2)
        expect(obj.groupIndices[0]).toBe(1);
        expect(obj.groupIndices[1]).toBe(2);
    });

    it('handles negative (relative) vertex indices', async () => {
        const parsed = await parseObj(objNegativeIndices).run();
        if (parsed.isError) throw new Error(parsed.message);
        const obj = parsed.result;

        expect(obj.triangleCount).toBe(1);
        // -3, -2, -1 with posCount=3 → 0, 1, 2
        expect(Array.from(obj.positionIndices)).toEqual([0, 1, 2]);
    });

    it('ignores comments and blank lines', async () => {
        const parsed = await parseObj(objWithComments).run();
        if (parsed.isError) throw new Error(parsed.message);
        const obj = parsed.result;

        expect(obj.positionCount).toBe(3);
        expect(obj.triangleCount).toBe(1);
    });

    it('silently skips unsupported directives', async () => {
        const parsed = await parseObj(objUnsupportedDirectives).run();
        if (parsed.isError) throw new Error(parsed.message);
        const obj = parsed.result;

        expect(obj.positionCount).toBe(3);
        expect(obj.triangleCount).toBe(1);
    });

    it('parses a cube (12 triangles, 8 vertices)', async () => {
        const parsed = await parseObj(objCube).run();
        if (parsed.isError) throw new Error(parsed.message);
        const obj = parsed.result;

        expect(obj.positionCount).toBe(8);
        expect(obj.triangleCount).toBe(12);
        expect(obj.positionIndices.length).toBe(36); // 12 * 3
    });

    it('returns no normals when none are defined', async () => {
        const parsed = await parseObj(objTriangle).run();
        if (parsed.isError) throw new Error(parsed.message);
        const obj = parsed.result;

        expect(obj.normalCount).toBe(0);
        // All normal indices should be -1
        expect(Array.from(obj.normalIndices)).toEqual([-1, -1, -1]);
    });

    it('default group is always present', async () => {
        const parsed = await parseObj(objTriangle).run();
        if (parsed.isError) throw new Error(parsed.message);
        const obj = parsed.result;

        expect(obj.groups[0]).toBe('default');
        expect(obj.groupIndices[0]).toBe(0);
    });

    it('parses CRLF line endings', async () => {
        const parsed = await parseObj(objCRLF).run();
        if (parsed.isError) throw new Error(parsed.message);
        const obj = parsed.result;

        expect(obj.positionCount).toBe(3);
        expect(obj.triangleCount).toBe(1);
        expect(Array.from(obj.positionIndices)).toEqual([0, 1, 2]);
    });

    it('handles tabs and leading whitespace before keywords', async () => {
        const parsed = await parseObj(objLeadingWhitespace).run();
        if (parsed.isError) throw new Error(parsed.message);
        const obj = parsed.result;

        expect(obj.positionCount).toBe(3);
        expect(obj.triangleCount).toBe(1);
        expect(Array.from(obj.positionIndices)).toEqual([0, 1, 2]);
    });

    it('skips a degenerate face and emits a warning', async () => {
        const parsed = await parseObj(objDegenerateFace).run();
        if (parsed.isError) throw new Error(parsed.message);
        const obj = parsed.result;

        // Only the valid triangle survives
        expect(obj.triangleCount).toBe(1);
        expect(Array.from(obj.positionIndices)).toEqual([0, 1, 2]);
        // A warning was recorded for the degenerate face
        expect(parsed.warnings.length).toBeGreaterThan(0);
        expect(parsed.warnings.some(w => w.includes('degenerate'))).toBe(true);
    });

    it('parses faces with mixed normal/no-normal vertices', async () => {
        const parsed = await parseObj(objMixedNormals).run();
        if (parsed.isError) throw new Error(parsed.message);
        const obj = parsed.result;

        expect(obj.positionCount).toBe(3);
        expect(obj.normalCount).toBe(1);
        expect(obj.triangleCount).toBe(1);
        // First and third vertices reference normal 0; the middle has none (-1)
        expect(Array.from(obj.normalIndices)).toEqual([0, -1, 0]);
    });

    it('handles negative (relative) normal indices', async () => {
        const parsed = await parseObj(objNegativeNormalIndices).run();
        if (parsed.isError) throw new Error(parsed.message);
        const obj = parsed.result;

        expect(obj.normalCount).toBe(1);
        expect(obj.triangleCount).toBe(1);
        // -1 with normCount=1 → 0
        expect(Array.from(obj.normalIndices)).toEqual([0, 0, 0]);
    });

    it('reuses an already-seen usemtl material name', async () => {
        const parsed = await parseObj(objReusedMaterial).run();
        if (parsed.isError) throw new Error(parsed.message);
        const obj = parsed.result;

        expect(obj.triangleCount).toBe(3);
        // 'red' (1) and 'green' (2) are the only added groups; no duplicate 'red'
        expect(obj.groups[1]).toBe('red');
        expect(obj.groups[2]).toBe('green');
        expect(obj.groups.indexOf('red')).toBe(1);
        expect(obj.groups.lastIndexOf('red')).toBe(1);
        // Third face reuses 'red' (index 1)
        expect(Array.from(obj.groupIndices)).toEqual([1, 2, 1]);
    });

    it('parses an empty file', async () => {
        const parsed = await parseObj(objEmpty).run();
        if (parsed.isError) throw new Error(parsed.message);
        const obj = parsed.result;

        expect(obj.positionCount).toBe(0);
        expect(obj.normalCount).toBe(0);
        expect(obj.triangleCount).toBe(0);
        expect(obj.groups[0]).toBe('default');
    });

    it('parses a file with vertices but no faces', async () => {
        const parsed = await parseObj(objNoFaces).run();
        if (parsed.isError) throw new Error(parsed.message);
        const obj = parsed.result;

        expect(obj.positionCount).toBe(3);
        expect(obj.triangleCount).toBe(0);
        expect(obj.positionIndices.length).toBe(0);
    });
});
