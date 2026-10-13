import assert from 'node:assert/strict';
import { mkdir, mkdtemp, rm, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import path from 'node:path';
import { test } from 'node:test';
import { checkMigrationRecords, exportNames } from '../migration-records.mjs';

const sources = {
  'packages/a/src/one.ts':
    'export const Alpha = 1;\nexport function beta() {}\nconst hidden = 2;\nexport { hidden as Gamma };\n',
  'packages/a/src/two.ts': "export * from './one.js';\nexport type Delta = number;\nexport namespace Ns {}\n",
};

async function withFixture(map, symbols, run) {
  const root = await mkdtemp(path.join(tmpdir(), 'molstar migration records '));
  try {
    const files = {
      ...sources,
      '.v6/plans/migration-map.json': JSON.stringify(map),
      '.v6/plans/migration-symbols.json': JSON.stringify(symbols),
    };
    for (const [file, text] of Object.entries(files)) {
      await mkdir(path.dirname(path.join(root, file)), { recursive: true });
      await writeFile(path.join(root, file), text);
    }
    return await run(root);
  } finally {
    await rm(root, { recursive: true, force: true });
  }
}

test('exportNames scans declarations, aliases and relative star re-exports', async () => {
  await withFixture({}, {}, (root) => {
    const names = exportNames(path.join(root, 'packages/a/src/two.ts'));
    assert.deepEqual([...names].sort(), ['Alpha', 'Delta', 'Gamma', 'Ns', 'beta']);
  });
});

test('accepts strings, arrays, nulls and exported symbols', async () => {
  const map = {
    'src/a.ts': 'packages/a/src/one.ts',
    'src/b.ts': ['packages/a/src/one.ts', 'packages/a/src/two.ts'],
    'src/c.ts': null,
  };
  const symbols = {
    'src/b.ts': {
      Alpha: 'packages/a/src/one.ts',
      Delta: 'packages/a/src/two.ts',
      Removed: null,
      'Legacy.Ns': 'packages/a/src/two.ts',
      'Legacy.BuiltIn': 'packages/a/src/one.ts',
    },
  };
  await withFixture(map, symbols, (root) => assert.deepEqual(checkMigrationRecords({ root }).errors, []));
});

test('reports missing targets, missing exports and malformed values', async () => {
  const map = {
    'src/a.ts': 'packages/a/src/gone.ts',
    'src/b.ts': ['packages/a/src/one.ts', 'packages/a/src/missing.ts'],
    'src/c.ts': ['packages/a/src/one.ts'],
  };
  const symbols = {
    'src/b.ts': {
      Alpha: 'packages/a/src/gone.ts',
      Nope: 'packages/a/src/one.ts',
      'Legacy.Nope': 'packages/a/src/one.ts',
      Bad: 5,
    },
  };
  await withFixture(map, symbols, (root) => {
    const { errors } = checkMigrationRecords({ root });
    const has = (text) =>
      assert.ok(
        errors.some((e) => e.includes(text)),
        `expected error containing "${text}"\n${errors.join('\n')}`,
      );
    has('src/a.ts: target does not exist: packages/a/src/gone.ts');
    has('src/b.ts: target does not exist: packages/a/src/missing.ts');
    has('src/c.ts: an array needs at least two distinct targets');
    has('src/b.ts: Alpha: target does not exist: packages/a/src/gone.ts');
    has('src/b.ts: Nope: not exported by packages/a/src/one.ts');
    has('src/b.ts: Legacy.Nope: not exported by packages/a/src/one.ts');
    has('src/b.ts: Bad: invalid target 5');
    assert.equal(errors.length, 7);
  });
});
