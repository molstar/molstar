/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { parseMtl } from '../../../mol-io/reader/obj/mtl-parser';
import { ObjFile } from '../../../mol-io/reader/obj/schema';
import { createObjShapeParams } from '../obj';
import { Color } from '../../../mol-util/color';
import { ColorNames } from '../../../mol-util/color/names';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';

// ─── helpers ─────────────────────────────────────────────────────────────────

function makeObjFile(materialNames: string[], vertexColorsLength = 0): ObjFile {
    return {
        positions: new Float32Array(0),
        normals: new Float32Array(0),
        positionIndices: new Int32Array(0),
        normalIndices: new Int32Array(0),
        positionCount: 0,
        normalCount: 0,
        triangleCount: 0,
        vertexColors: new Float32Array(vertexColorsLength),
        materialNames,
        faceGroups: new Int32Array(0),
    };
}

function coloringOptionKeys(params: ReturnType<typeof createObjShapeParams>): string[] {
    return (params.coloring as PD.Mapped<any>).select.options.map((o => o[0]));
}

// ─── createObjShapeParams ────────────────────────────────────────────────────

describe('createObjShapeParams', () => {
    it('no materials, no MTL → uniform is default, no given or custom options', () => {
        const params = createObjShapeParams(makeObjFile([]));
        const keys = coloringOptionKeys(params);
        expect(params.coloring.defaultValue.name).toBe('uniform');
        expect(keys).toContain('uniform');
        expect(keys).not.toContain('given');
        // custom is always present
        expect(keys).toContain('custom');
    });

    it('materials, no MTL → custom is default, no given option', () => {
        const obj = makeObjFile(['red', 'green']);
        const params = createObjShapeParams(obj);
        const keys = coloringOptionKeys(params);
        expect(params.coloring.defaultValue.name).toBe('custom');
        expect(keys).not.toContain('given');
        expect(keys).toContain('custom');
    });

    it('materials + MTL → given is default, both given and custom present', () => {
        const obj = makeObjFile(['red', 'green']);
        const mtl = parseMtl('newmtl red\nKd 1 0 0\nnewmtl green\nKd 0 1 0\n');
        const params = createObjShapeParams(obj, mtl);
        const keys = coloringOptionKeys(params);
        expect(params.coloring.defaultValue.name).toBe('given');
        expect(keys).toContain('given');
        expect(keys).toContain('custom');
    });

    it('MTL provided but no materials in OBJ → no given option, uniform is default', () => {
        const obj = makeObjFile([]);
        const mtl = parseMtl('newmtl red\nKd 1 0 0\n');
        const params = createObjShapeParams(obj, mtl);
        const keys = coloringOptionKeys(params);
        expect(params.coloring.defaultValue.name).toBe('uniform');
        expect(keys).not.toContain('given');
    });

    it('vertex colors present → vertex is default, given still added when MTL provided', () => {
        const obj = makeObjFile(['red'], 9 /* 3 vertices × 3 components */);
        const mtl = parseMtl('newmtl red\nKd 1 0 0\n');
        const params = createObjShapeParams(obj, mtl);
        const keys = coloringOptionKeys(params);
        expect(params.coloring.defaultValue.name).toBe('vertex');
        expect(keys).toContain('given');
        expect(keys).toContain('vertex');
    });

    it('custom option uses distinctColors defaults, not MTL colors', () => {
        const obj = makeObjFile(['mat1']);
        const mtl = parseMtl('newmtl mat1\nKd 1 0 0\n');
        const params = createObjShapeParams(obj, mtl);
        // custom group params for 'mat1' should be grey (distinctColors fallback for 1 material)
        // not the MTL red (1,0,0)
        const customGroup = (params.coloring as PD.Mapped<any>).map('custom') as PD.Group<any>;
        const customDefault = customGroup.defaultValue as Record<string, Color>;
        expect(customDefault['mat1']).toBe(ColorNames.grey);
        expect(customDefault['mat1']).not.toBe(Color.fromNormalizedRgb(1, 0, 0));
    });

    it('given option has no user params (empty group)', () => {
        const obj = makeObjFile(['red']);
        const mtl = parseMtl('newmtl red\nKd 1 0 0\n');
        const params = createObjShapeParams(obj, mtl);
        const givenGroup = (params.coloring as PD.Mapped<any>).map('given') as PD.Group<any>;
        expect(Object.keys(givenGroup.params)).toHaveLength(0);
    });
});

// ─── given coloring mode: color resolution ───────────────────────────────────

describe('given coloring mode via PD.getDefaultValues', () => {
    it('default values with given mode use MTL Kd colors', () => {
        const obj = makeObjFile(['red', 'green']);
        const mtl = parseMtl('newmtl red\nKd 1 0 0\nnewmtl green\nKd 0 1 0\n');
        const params = createObjShapeParams(obj, mtl);
        const defaults = PD.getDefaultValues(params);

        // Default coloring is 'given'
        expect(defaults.coloring.name).toBe('given');

        // The 'given' default value reflects the MTL colors indirectly through the
        // coloring name — actual color lookup happens at render time via getMaterialColors.
        // Here we just verify the option structure is correct.
        expect(defaults.coloring.params).toEqual({});
    });

    it('default values with custom mode expose per-material color params', () => {
        const obj = makeObjFile(['alpha', 'beta']);
        const params = createObjShapeParams(obj); // no MTL
        const defaults = PD.getDefaultValues(params);

        expect(defaults.coloring.name).toBe('custom');
        expect(defaults.coloring.params).toHaveProperty('alpha');
        expect(defaults.coloring.params).toHaveProperty('beta');
    });
});
