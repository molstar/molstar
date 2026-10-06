import assert from 'node:assert/strict';
import { mkdir, mkdtemp, rm, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import path from 'node:path';
import { test } from 'node:test';
import { checkImportGraph, classify } from '../import-graph.mjs';

const LIB = 'packages/lib';
const HIGH = 'packages/high';
const APP = 'apps/app';
const packages = [
  { name: '@molstar/lib', path: LIB, kind: 'library' },
  { name: '@molstar/high', path: HIGH, kind: 'library' },
  { name: '@molstar/app', path: APP, kind: 'app' },
];

const baseManifest = () => ({
  layers: ['@molstar/lib', '@molstar/high'],
  catalogs: [`${LIB}/src/catalog.ts`],
  defaultComposition: [`${LIB}/src/default-spec.ts`],
  baseEntryPoints: [`${LIB}/src/context.ts`],
  allowlist: [],
});

/** Writes a synthetic workspace: package manifests with a wildcard `molstar-src` export plus the given sources. */
async function withFixture(files, manifest, run) {
  const root = await mkdtemp(path.join(tmpdir(), 'molstar import graph '));
  try {
    for (const pkg of packages)
      files[`${pkg.path}/package.json`] = JSON.stringify({
        name: pkg.name,
        exports: { './*': { 'molstar-src': './src/*.ts', import: './lib/*.js' } },
      });
    for (const [file, text] of Object.entries(files)) {
      await mkdir(path.dirname(path.join(root, file)), { recursive: true });
      await writeFile(path.join(root, file), text);
    }
    return await run(checkImportGraph({ root, packages, manifest }));
  } finally {
    await rm(root, { recursive: true, force: true });
  }
}

const clean = () => ({
  [`${LIB}/src/catalog.ts`]: 'export const Catalog = [1];\n',
  [`${LIB}/src/default-spec.ts`]: "import { Catalog } from './catalog.js';\nexport const Spec = Catalog;\n",
  [`${LIB}/src/context.ts`]: 'export const Context = 1;\n',
  [`${APP}/src/index.ts`]: "import { Spec } from '@molstar/lib/default-spec';\nconsole.log(Spec);\n",
});

test('classification: catalog, default-composition, app and base modules', () => {
  const manifest = baseManifest();
  assert.equal(classify(`${LIB}/src/catalog.ts`, manifest), 'catalog');
  assert.equal(classify(`${LIB}/src/default-spec.ts`, manifest), 'default-composition');
  for (const dir of ['apps/a', 'examples/b', 'cli/c', 'servers/d', 'smoke/e'])
    assert.equal(classify(`${dir}/src/x.ts`, manifest), 'app');
  assert.equal(classify(`${LIB}/src/context.ts`, manifest), 'base');
  assert.equal(classify('extensions/x/src/index.ts', manifest), 'base');
});

test('catalogs, default-composition modules and apps may value-import catalogs', async () => {
  const files = clean();
  files[`${APP}/src/direct.ts`] = "import { Catalog } from '@molstar/lib/catalog';\nconsole.log(Catalog);\n";
  files[`${LIB}/src/other-catalog.ts`] = "import { Catalog } from './catalog.js';\nexport const Other = Catalog;\n";
  const manifest = baseManifest();
  manifest.catalogs.push(`${LIB}/src/other-catalog.ts`);
  await withFixture(files, manifest, ({ errors }) => assert.deepEqual(errors, []));
});

test('a value import of a catalog from a base module is a violation, a type import is not', async () => {
  const files = clean();
  files[`${LIB}/src/base.ts`] = "import { Catalog } from './catalog.js';\nexport const B = Catalog;\n";
  files[`${LIB}/src/typed.ts`] = "import type { Catalog } from './catalog.js';\nexport type T = typeof Catalog;\n";
  files[`${LIB}/src/inline.ts`] =
    "import { type Catalog } from './catalog.js';\nexport function f(c: typeof Catalog) { return c; }\n";
  await withFixture(files, baseManifest(), ({ errors }) => {
    assert.equal(errors.length, 1, errors.join('\n'));
    assert.match(errors[0], /rule a: packages\/lib\/src\/base\.ts value-imports catalog module/);
  });
});

test('package subpath imports resolve to sources and are checked', async () => {
  const files = clean();
  files[`${HIGH}/src/uses.ts`] = "import { Catalog } from '@molstar/lib/catalog';\nexport const U = Catalog;\n";
  await withFixture(files, baseManifest(), ({ errors }) => {
    assert.equal(errors.length, 1, errors.join('\n'));
    assert.match(
      errors[0],
      /packages\/high\/src\/uses\.ts value-imports catalog module packages\/lib\/src\/catalog\.ts/,
    );
  });
});

test('base entry points must not reach default-composition modules transitively', async () => {
  const files = clean();
  files[`${LIB}/src/context.ts`] = "import { helper } from './helper.js';\nexport const Context = helper;\n";
  files[`${LIB}/src/helper.ts`] = "import { Spec } from './default-spec.js';\nexport const helper = Spec;\n";
  await withFixture(files, baseManifest(), ({ errors }) => {
    assert.equal(errors.length, 1, errors.join('\n'));
    assert.match(errors[0], /rule b: base entry point packages\/lib\/src\/context\.ts reaches default-composition/);
    assert.match(errors[0], /context\.ts -> packages\/lib\/src\/helper\.ts -> packages\/lib\/src\/default-spec\.ts/);
  });
});

test('allowlisted violations pass', async () => {
  const files = clean();
  files[`${LIB}/src/base.ts`] = "import { Catalog } from './catalog.js';\nexport const B = Catalog;\n";
  files[`${LIB}/src/context.ts`] = "import { Spec } from './default-spec.js';\nexport const Context = Spec;\n";
  const manifest = baseManifest();
  manifest.allowlist = [
    {
      rule: 'a',
      from: `${LIB}/src/base.ts`,
      to: `${LIB}/src/catalog.ts`,
      reason: 'preload',
      step: 'step 3',
    },
    {
      rule: 'b',
      from: `${LIB}/src/context.ts`,
      to: `${LIB}/src/default-spec.ts`,
      reason: 'fallback',
      step: 'step 3',
    },
  ];
  await withFixture(files, manifest, ({ errors, violations }) => {
    assert.deepEqual(errors, []);
    assert.equal(violations.length, 2);
  });
});

test('stale allowlist entries, duplicates and entries without a reason fail', async () => {
  const manifest = baseManifest();
  const entry = {
    rule: 'a',
    from: `${LIB}/src/base.ts`,
    to: `${LIB}/src/catalog.ts`,
    reason: 'preload',
    step: 'step 3',
  };
  manifest.allowlist = [entry, entry, { ...entry, from: `${LIB}/src/gone.ts`, reason: '' }];
  await withFixture(clean(), manifest, ({ errors }) => {
    assert.ok(
      errors.some((e) => /stale allowlist entry.*base\.ts/.test(e)),
      errors.join('\n'),
    );
    assert.ok(
      errors.some((e) => /stale allowlist entry.*gone\.ts/.test(e)),
      errors.join('\n'),
    );
    assert.ok(
      errors.some((e) => /duplicate allowlist entry/.test(e)),
      errors.join('\n'),
    );
    assert.ok(
      errors.some((e) => /needs a reason and a step.*gone\.ts/.test(e)),
      errors.join('\n'),
    );
  });
});

test('type-only imports to a higher package are reported (rule d)', async () => {
  const files = clean();
  files[`${LIB}/src/lower.ts`] = "import type { High } from '@molstar/high/high';\nexport type L = typeof High;\n";
  files[`${HIGH}/src/high.ts`] = 'export const High = 1;\n';
  files[`${HIGH}/src/upper.ts`] =
    "import type { Context } from '@molstar/lib/context';\nexport type U = typeof Context;\n";
  await withFixture(files, baseManifest(), ({ errors }) => {
    assert.equal(errors.length, 1, errors.join('\n'));
    assert.match(
      errors[0],
      /rule d: packages\/lib\/src\/lower\.ts \(@molstar\/lib\) type-imports packages\/high\/src\/high\.ts/,
    );
  });
  const manifest = baseManifest();
  manifest.allowlist = [
    { rule: 'd', from: `${LIB}/src/lower.ts`, to: `${HIGH}/src/high.ts`, reason: 'planned', step: 'step 2' },
  ];
  await withFixture(files, manifest, ({ errors }) => assert.deepEqual(errors, []));
});

test('manifest paths must exist and every catalog.ts must be listed', async () => {
  const files = clean();
  files[`${LIB}/src/nested/catalog.ts`] = 'export const Nested = 1;\n';
  const manifest = baseManifest();
  manifest.catalogs.push(`${LIB}/src/missing.ts`);
  await withFixture(files, manifest, ({ errors }) => {
    assert.ok(
      errors.some((e) => /does not exist.*missing\.ts/.test(e)),
      errors.join('\n'),
    );
    assert.ok(
      errors.some((e) => /missing from the manifest.*nested\/catalog\.ts/.test(e)),
      errors.join('\n'),
    );
  });
});

test('ignored directories do not contribute edges', async () => {
  const files = clean();
  files[`${LIB}/src/_test/base.ts`] = "import { Catalog } from '../catalog.js';\nexport const B = Catalog;\n";
  files[`${LIB}/lib/base.ts`] = "import { Catalog } from '../src/catalog.js';\nexport const B = Catalog;\n";
  await withFixture(files, baseManifest(), ({ errors }) => assert.deepEqual(errors, []));
});
