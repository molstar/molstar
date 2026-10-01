import assert from 'node:assert/strict';
import { mkdtemp, mkdir, writeFile, rm, realpath } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import path from 'node:path';
import { spawnSync } from 'node:child_process';
import { test } from 'node:test';
import { pathToFileURL } from 'node:url';
import { compactExports, expandExports, resolveExport } from '../exports.mjs';

const moduleEntry = (name, ext = 'ts') => ({ 'molstar-src': `./src/${name}.${ext}`, types: `./lib/${name}.d.ts`, import: `./lib/${name}.js` });

test('compaction preserves nested exports, aliases, mixed TS/TSX, Sass and private exclusions', () => {
    const original = {
        '.': moduleEntry('index'),
        './index': moduleEntry('index'),
        './deep/one': moduleEntry('deep/one'),
        './view': moduleEntry('view', 'tsx'),
        './task': moduleEntry('task/index'),
        './task/index': moduleEntry('task/index'),
        './ambient': { types: './src/ambient.d.ts' },
    };
    for (const skin of ['light', 'dark']) {
        original[`./skin/${skin}.scss`] = { 'molstar-src': `./src/skin/${skin}.scss`, sass: `./src/skin/${skin}.scss`, default: `./lib/skin/${skin}.scss` };
        original[`./skin/${skin}.css`] = `./lib/skin/${skin}.css`;
    }
    const sources = ['index.ts', 'deep/one.ts', 'view.tsx', 'task/index.ts', 'ambient.d.ts', 'deep/_test/private.test.ts', 'skin/light.scss', 'skin/dark.scss'];
    const files = sources.map(file => `src/${file}`).concat(['index', 'deep/one', 'view', 'task/index'].flatMap(name => [`lib/${name}.js`, `lib/${name}.d.ts`]), ['light', 'dark'].map(name => `lib/skin/${name}.css`));
    const compact = compactExports(original, sources);
    assert(compact['./*']);
    assert(compact['./*.scss']);
    assert(compact['./*.css']);
    assert.deepEqual(expandExports(compact, files), original);
    assert.equal(resolveExport(compact, './deep/_test/private.test'), null);
    assert.equal(resolveExport(compact, './ambient.d'), null);
    assert.deepEqual(compactExports(compact, sources), compact);
});

test('wildcard resolution matches Node for exact overrides, nested patterns and null exclusions', async () => {
    const dir = await realpath(await mkdtemp(path.join(tmpdir(), 'molstar-export-patterns-')));
    try {
        const exports = {
            './*': moduleEntry('*'),
            './ui/*': moduleEntry('ui/*', 'tsx'),
            './ui/special': moduleEntry('special'),
            './*.scss': { sass: './src/*.scss', default: './lib/*.scss' },
            './private/*': null,
            './ui/_test/*': null,
        };
        await mkdir(path.join(dir, 'node_modules/fixture'), { recursive: true });
        await writeFile(path.join(dir, 'node_modules/fixture/package.json'), JSON.stringify({ name: 'fixture', type: 'module', exports }));
        const keys = ['./deep/nested', './ui/view', './ui/special', './skin/light.scss', './private/hidden', './ui/_test/view'];
        for (const condition of ['import', 'molstar-src', 'types', 'sass']) {
            const script = `const keys = ${JSON.stringify(keys)}; console.log(JSON.stringify(keys.map(key => { try { return import.meta.resolve('fixture/' + key.slice(2)); } catch (error) { return error.code; } })));`;
            const result = spawnSync(process.execPath, [`--conditions=${condition}`, '--input-type=module', '-e', script], { cwd: dir, encoding: 'utf8' });
            assert.equal(result.status, 0, result.stderr);
            const expected = keys.map(key => {
                const entry = resolveExport(exports, key);
                if (entry === null) return 'ERR_PACKAGE_PATH_NOT_EXPORTED';
                const target = entry[condition] ?? entry.import ?? entry.default;
                return new URL(target, pathToFileURL(path.join(dir, 'node_modules/fixture/'))).href;
            });
            assert.deepEqual(JSON.parse(result.stdout), expected);
        }
    } finally { await rm(dir, { recursive: true, force: true }); }
});

test('expansion rejects empty patterns and targets outside the package', () => {
    assert.throws(() => expandExports({ './*': { import: './lib/*.js' } }, ['src/a.ts']), /no matching files/);
    assert.throws(() => expandExports({ './*': '../*.js' }, ['src/a.ts']), /package-local/);
});
