import assert from 'node:assert/strict';
import { mkdtemp, writeFile, rm } from 'node:fs/promises';
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
