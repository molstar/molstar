import assert from 'node:assert/strict';
import { mkdir, mkdtemp, rm, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import path from 'node:path';
import { test } from 'node:test';
import { buildEntry, checkMetafile, importChain, outputBytes } from '../slim-bundle.mjs';

const input = (...imports) => ({ bytes: 1, imports: imports.map((p) => ({ path: p, kind: 'import-statement' })) });

const metafile = {
  inputs: {
    'app/index.ts': input('lib/a.ts', 'lib/b.ts'),
    'lib/a.ts': input('lib/c.ts'),
    'lib/b.ts': input('lib/c.ts', 'lib/heavy/x.ts'),
    'lib/c.ts': input(),
    'lib/heavy/x.ts': input('lib/heavy/y.ts'),
    'lib/heavy/y.ts': input(),
  },
  outputs: {
    'out/index.js': { bytes: 100, entryPoint: 'app/index.ts', inputs: {} },
    'out/chunk.js': { bytes: 20, inputs: {} },
    'out/index.js.map': { bytes: 999, inputs: {} },
  },
};

test('importChain finds the shortest chain from the entry point', () => {
  assert.deepEqual(importChain(metafile, 'lib/heavy/y.ts'), [
    'app/index.ts',
    'lib/b.ts',
    'lib/heavy/x.ts',
    'lib/heavy/y.ts',
  ]);
  assert.deepEqual(importChain(metafile, 'lib/c.ts'), ['app/index.ts', 'lib/a.ts', 'lib/c.ts']);
});

test('outputBytes sums the JavaScript outputs', () => {
  assert.equal(outputBytes(metafile), 120);
});

test('checkMetafile reports excluded inputs with their chain and skips known leaks', () => {
  const exclusions = { excluded: [{ name: 'Heavy', paths: ['lib/heavy/*.ts'] }], knownLeaks: [] };
  const { errors } = checkMetafile(metafile, exclusions, 'test');
  assert.equal(errors.length, 2);
  assert.match(
    errors[1],
    /^test: excluded module lib\/heavy\/y\.ts \(Heavy\) is in the bundle, imported through: app\/index\.ts -> lib\/b\.ts -> lib\/heavy\/x\.ts -> lib\/heavy\/y\.ts$/,
  );

  exclusions.knownLeaks = [{ path: 'lib/heavy/*.ts', reason: 'r', decision: 'd' }];
  const known = checkMetafile(metafile, exclusions, 'test');
  assert.deepEqual(known.errors, []);
  assert.deepEqual([...known.known.keys()], ['lib/heavy/x.ts', 'lib/heavy/y.ts']);
});

test('buildEntry bundles with esbuild in both modes and the metafile shows what it contains', async () => {
  const root = await mkdtemp(path.join(tmpdir(), 'molstar slim bundle '));
  try {
    const files = {
      'version.json': JSON.stringify({ version: '0.0.0' }),
      'tsconfig.json': '{}',
      'app/index.ts': "import { a } from './a.js';\nconsole.log(a);\nexport const lazy = () => import('./lazy.js');\n",
      'app/a.ts': "import { heavy } from './heavy.js';\nexport const a = heavy;\n",
      'app/heavy.ts': 'export const heavy = 1;\n',
      'app/lazy.ts': 'export const lazy = 2;\n',
      'app/style.scss': 'a { color: red; }\n',
    };
    for (const [file, text] of Object.entries(files)) {
      await mkdir(path.dirname(path.join(root, file)), { recursive: true });
      await writeFile(path.join(root, file), text);
    }
    const exclusions = {
      excluded: [{ name: 'Heavy', paths: ['app/heavy.ts', 'app/lazy.ts'] }],
      knownLeaks: [],
    };
    for (const mode of ['split ESM', 'single file']) {
      const result = await buildEntry({ root, entry: 'app/index.ts', mode });
      assert.ok(outputBytes(result.metafile) > 0);
      const { errors } = checkMetafile(result.metafile, exclusions, mode);
      assert.equal(errors.length, 2, errors.join('\n'));
      assert.match(
        errors.find((e) => /heavy\.ts/.test(e)),
        /app\/index\.ts -> app\/a\.ts -> app\/heavy\.ts$/,
      );
      // a dynamic import is part of both builds: a chunk when splitting, inlined otherwise
      assert.match(
        errors.find((e) => /lazy\.ts/.test(e)),
        /app\/index\.ts -> app\/lazy\.ts$/,
      );
    }
    const split = await buildEntry({ root, entry: 'app/index.ts', mode: 'split ESM' });
    assert.ok(Object.keys(split.metafile.outputs).length > 1, 'splitting creates a chunk for the dynamic import');
  } finally {
    await rm(root, { recursive: true, force: true });
  }
});
