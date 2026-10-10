/**
 * Copyright (c) 2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import assert from 'node:assert/strict';
import { createRequire } from 'node:module';
import { spawnSync } from 'node:child_process';
import { writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { MVSData } from '@molstar/mvs-builder';
import { MolScriptBuilder as B, compileScript } from '@molstar/mvs-builder/molql';
import { expressionValidationIssues } from '@molstar/query-language/language/validation';

const require = createRequire(import.meta.url);
for (const name of [
  '@molstar/model/model/structure',
  '@molstar/graphics/canvas3d/canvas3d',
  '@molstar/plugin/context',
  'gl',
  'canvas',
]) {
  assert.throws(
    () => require.resolve(name),
    { code: 'MODULE_NOT_FOUND' },
    `Standalone MVS authoring installed ${name}`,
  );
}
function state(expression) {
  const builder = MVSData.createBuilder();
  builder
    .download({ url: 'example.bcif' })
    .parse({ format: 'bcif' })
    .modelStructure()
    .component({ selector: { molql: expression } });
  return builder.getState();
}
const files = [];
for (const language of ['mol-script', 'pymol', 'vmd', 'jmol']) {
  const expression = compileScript(language, language === 'mol-script' ? '(sel.atom.all)' : 'all');
  assert.equal(expressionValidationIssues(expression), undefined);
  const data = state(expression);
  assert(MVSData.isValid(MVSData.fromMVSJ(MVSData.toMVSJ(data))));
  const file = `${language}.mvsj`;
  writeFileSync(file, MVSData.toMVSJ(data));
  files.push(file);
}
assert.equal(expressionValidationIssues(B.re('^C')), undefined);
assert.throws(() => compileScript('pymol', '('));
const cli = fileURLToPath(new URL('../bin/mvs-validate.mjs', import.meta.resolve('@molstar/mvs-builder')));
const valid = spawnSync(process.execPath, [cli, ...files], { encoding: 'utf8' });
assert.equal(valid.status, 0, valid.stderr);
assert.equal(valid.stdout.split('\n').filter((line) => line.startsWith('OK')).length, 4);
const invalid = [
  [
    'unknown.mvsj',
    MVSData.toMVSJ(state(B.struct.generator.atomGroups({ 'atom-test': { head: { name: 'unknown' } } }))),
  ],
  ['argument.mvsj', MVSData.toMVSJ(state(B.struct.generator.atomGroups({ 'atom-tset': true })))],
  ['syntax.mvsj', MVSData.toMVSJ(state({ head: { name: 'core.math.add' }, arguments: [] }))],
  ['snapshot.mvsj', MVSData.toMVSJ(MVSData.stateToStates(state({ head: { name: 'unknown' } })))],
  ['broken.mvsj', '{'],
];
for (const [file, data] of invalid) writeFileSync(file, data);
const failed = spawnSync(process.execPath, [cli, ...invalid.map(([file]) => file), 'missing.mvsj', files[0]], {
  encoding: 'utf8',
});
assert.equal(failed.status, 1, failed.stderr);
for (const [file] of invalid) assert(failed.stdout.includes(`FAILED ${file}`), failed.stdout);
assert.match(failed.stderr, /Unknown callable symbol 'unknown'/);
assert.match(failed.stderr, /Unknown argument 'atom-tset'/);
assert.match(failed.stderr, /Unknown application field 'arguments'/);
assert.match(failed.stderr, /snapshots\[0\]\.root/);
assert(failed.stdout.includes('FAILED missing.mvsj'), failed.stdout);
assert(failed.stdout.includes(`OK     ${files[0]}`), failed.stdout);
console.log('Packed standalone MVS MolQL authoring and CLI validation passed for all four languages.');
