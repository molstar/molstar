/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { ComputeShaders, ShaderModules } from './compute-shaders.js';
import { GraphicsShaders } from './graphics-shaders.js';

const SourceRoot = 'packages/graphics/src';

/** schema and shader code are generated per call, see the 'grid3d compute' tests */
const Dynamic = new Set(['gl/compute/grid3d.ts']);

function getSourceFiles(dir: string, files: string[] = []) {
    for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
        const p = path.join(dir, entry.name);
        if (entry.isDirectory()) {
            if (entry.name !== '_test' && entry.name !== 'shader') getSourceFiles(p, files);
        } else if (entry.name.endsWith('.ts')) {
            files.push(p);
        }
    }
    return files;
}

describe('shader list', () => {
    const definitions: [file: string, name: string][] = [];
    const unlisted: string[] = [];
    for (const p of getSourceFiles(SourceRoot)) {
        const file = path.relative(SourceRoot, p).split(path.sep).join('/');
        const lines = fs.readFileSync(p, 'utf8').split('\n');
        lines.forEach((line, i) => {
            if (!/\bShaderCode\(/.test(line) || /function ShaderCode\(/.test(line) || Dynamic.has(file)) return;
            const m = line.match(/^export const (\w+) = ShaderCode\(/);
            if (m) definitions.push([file, m[1]]);
            else unlisted.push(`${file}:${i + 1}`);
        });
    }

    it('shader codes are exported module-level consts', () => {
        expect(unlisted).toEqual([]);
    });

    it('all shader codes are listed', () => {
        const listed = new Set([...GraphicsShaders, ...ComputeShaders].map(([, shaderCode]) => shaderCode));
        const defined = new Set(
            definitions.map(([file, name]) => (ShaderModules[file] as Record<string, any> | undefined)?.[name]),
        );
        const missing = definitions
            .filter(([file, name]) => !listed.has((ShaderModules[file] as Record<string, any> | undefined)?.[name]))
            .map(([file, name]) => `${file}#${name}`);
        expect(missing).toEqual([]);
        expect([...listed].filter((shaderCode) => !defined.has(shaderCode)).map((s) => s.name)).toEqual([]);
    });

    it('names are unique', () => {
        const names = [...GraphicsShaders, ...ComputeShaders].map(([name]) => name);
        expect(names.length).toBe(new Set(names).size);
    });

    it('program names are unique', () => {
        // the name keys the per-backend source registries
        const shaderCodes = new Set([...GraphicsShaders, ...ComputeShaders].map(([, shaderCode]) => shaderCode));
        const names = [...shaderCodes].map((shaderCode) => shaderCode.name);
        expect(names.filter((name, i) => names.indexOf(name) !== i)).toEqual([]);
    });
});
