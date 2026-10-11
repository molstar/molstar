/**
 * Copyright (c) 2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** Ensure MolQL authoring/validation never loads the molecular query runtime or unrelated text languages. */
import assert from 'node:assert/strict';
import { build } from 'esbuild';

async function graph(entry) {
  const result = await build({
    entryPoints: [entry],
    bundle: true,
    platform: 'node',
    format: 'esm',
    conditions: ['molstar-src'],
    write: false,
    metafile: true,
    logLevel: 'silent',
  });
  const inputs = Object.keys(result.metafile.inputs).map((path) => path.replaceAll('\\', '/'));
  const forbidden = inputs.filter((path) => /^packages\/(model|graphics|plugin|mvs\/runtime)\//.test(path));
  assert.deepEqual(forbidden, [], `${entry} loads molecular/rendering runtime modules`);
  return inputs;
}

for (const entry of [
  'packages/query-language/src/language/builder.ts',
  'packages/query-language/src/language/validation.ts',
  'packages/mvs/builder/src/mvs-data.ts',
]) {
  const inputs = await graph(entry);
  assert.deepEqual(
    inputs.filter((path) =>
      /^packages\/query-language\/src\/(transpilers|script\/|script\.ts|compile\.ts|language\/parser\.ts)/.test(path),
    ),
    [],
    `${entry} loads text parsers/transpilers`,
  );
}
for (const language of ['pymol', 'vmd', 'jmol']) {
  const inputs = await graph(`packages/query-language/src/transpilers/${language}/parser.ts`);
  const unrelated = inputs.filter((path) =>
    ['pymol', 'vmd', 'jmol'].some((other) => other !== language && path.includes(`/transpilers/${other}/`)),
  );
  assert.deepEqual(unrelated, [], `${language} loads another text language`);
}
const authoring = await graph('packages/mvs/builder/src/molql.ts');
for (const language of ['pymol', 'vmd', 'jmol']) {
  assert(
    authoring.some((path) => path.endsWith(`/transpilers/${language}/parser.ts`)),
    `MVS authoring is missing ${language}`,
  );
}
assert(
  authoring.some((path) => path.endsWith('/language/parser.ts')),
  'MVS authoring is missing MolScript',
);
console.log('MolQL builder, validation, and text-language import boundaries passed.');
