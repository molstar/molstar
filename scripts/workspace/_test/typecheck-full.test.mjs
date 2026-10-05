import assert from 'node:assert/strict';
import { mkdtemp, mkdir, writeFile, readFile, readdir, rm } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import path from 'node:path';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { test } from 'node:test';

test('publish type check rechecks imported declarations despite normal build caches', async () => {
    const dir = await mkdtemp(path.join(tmpdir(), 'molstar-full-check-'));
    const root = fileURLToPath(new URL('../../../', import.meta.url));
    try {
        const config = path.join(dir, 'tsconfig.json');
        await writeFile(config, JSON.stringify({
            compilerOptions: { composite: true, skipLibCheck: true, types: [], outDir: './lib', tsBuildInfoFile: './lib/tsconfig.tsbuildinfo' },
            files: ['index.ts', 'dependency.d.ts'],
        }));
        await writeFile(path.join(dir, 'index.ts'), "import type { Value } from './dependency';\nexport const value: Value = 1;\n");
        await writeFile(path.join(dir, 'dependency.d.ts'), 'export type Value = number;\nexport declare const invalid: MissingType;\n');
        const normal = () => spawnSync(process.execPath, [path.join(root, 'node_modules/typescript/bin/tsc'), '-b', config], { encoding: 'utf8' });
        assert.equal(normal().status, 0);
        assert.equal(normal().status, 0);
        const full = spawnSync(process.execPath, [path.join(root, 'scripts/workspace/typecheck-full.mjs'), config], { encoding: 'utf8' });
        assert.notEqual(full.status, 0);
        assert.match(full.stdout + full.stderr, /MissingType/);
    } finally {
        await rm(dir, { recursive: true, force: true });
    }
});

test('full check overrides inherited settings throughout references and preserves normal caches', async () => {
    const dir = await mkdtemp(path.join(tmpdir(), 'molstar-full-graph-'));
    const root = fileURLToPath(new URL('../../../', import.meta.url));
    const compiler = path.join(root, 'node_modules/typescript/bin/tsc');
    const fullCheck = path.join(root, 'scripts/workspace/typecheck-full.mjs');
    try {
        await mkdir(path.join(dir, 'dependency'));
        await mkdir(path.join(dir, 'consumer'));
        await writeFile(path.join(dir, 'base.json'), JSON.stringify({
            compilerOptions: { composite: true, skipLibCheck: true, types: [] },
        }));
        await writeFile(path.join(dir, 'tsconfig.json'), JSON.stringify({
            files: [], references: [{ path: './consumer' }],
        }));
        for (const project of ['dependency', 'consumer']) {
            await writeFile(path.join(dir, project, 'tsconfig.json'), JSON.stringify({
                extends: '../base.json',
                compilerOptions: { outDir: './lib', tsBuildInfoFile: './lib/tsconfig.tsbuildinfo' },
                files: project === 'dependency' ? ['index.ts', 'external.d.ts'] : ['index.ts'],
                references: project === 'consumer' ? [{ path: '../dependency/tsconfig.json' }] : [],
            }));
        }
        await writeFile(path.join(dir, 'dependency/index.ts'), 'export const value = 1;\n');
        await writeFile(path.join(dir, 'dependency/external.d.ts'), 'export declare const invalid: MissingType;\n');
        await writeFile(path.join(dir, 'consumer/index.ts'), "export { value } from '../dependency';\n");
        const config = path.join(dir, 'tsconfig.json');
        const normal = spawnSync(process.execPath, [compiler, '-b', config], { encoding: 'utf8' });
        assert.equal(normal.status, 0, normal.stdout + normal.stderr);
        const caches = await Promise.all(['dependency', 'consumer'].map(project => readFile(path.join(dir, project, 'lib/tsconfig.tsbuildinfo'), 'utf8')));
        const runFull = () => spawnSync(process.execPath, [fullCheck, config], { encoding: 'utf8' });
        const failed = runFull();
        assert.notEqual(failed.status, 0);
        assert.match(failed.stdout + failed.stderr, /MissingType/);
        await writeFile(path.join(dir, 'dependency/external.d.ts'), 'export declare const valid: number;\n');
        const passed = runFull();
        assert.equal(passed.status, 0, passed.stdout + passed.stderr);
        for (const [index, project] of ['dependency', 'consumer'].entries()) {
            assert.equal(await readFile(path.join(dir, project, 'lib/tsconfig.tsbuildinfo'), 'utf8'), caches[index]);
            assert((await readdir(path.join(dir, project, 'lib'))).includes('tsconfig.full.tsbuildinfo'));
        }
        for (const project of ['', 'dependency', 'consumer']) {
            assert(!(await readdir(path.join(dir, project))).some(name => name.startsWith('.tsconfig.full-')));
        }
    } finally {
        await rm(dir, { recursive: true, force: true });
    }
});
