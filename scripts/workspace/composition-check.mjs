#!/usr/bin/env node
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Acceptance checks of the plugin-composition refactor (.v6/plans/plugin-composition.md, section 6) that need the
 * Viewer's own spec, which the jest setup cannot load (its UI extensions are `.tsx` modules), and the CLI tools:
 *   - a snapshot made by the default plugin restores in a plugin with the Viewer spec, and one that needs an
 *     unregistered transformer fails before changing the Viewer plugin's state,
 *   - `customFormats` of the Viewer spec can override a built-in format, and the custom provider loads the data,
 *   - `loadTrajectory({ preset: 'all-models' })` works through the Viewer spec,
 *   - the Viewer spec evaluates a PyMOL script,
 *   - `cli/state-docs` runs and `cli/mvs-render --help` prints its usage.
 *
 * Usage: node scripts/workspace/composition-check.mjs
 *
 * Requires compiled library output (`pnpm build:lib`). Nothing is written to the repository.
 */

import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { mkdirSync, mkdtempSync, readFileSync, rmSync } from 'node:fs';
import { registerHooks } from 'node:module';
import { tmpdir } from 'node:os';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const root = resolve(dirname(fileURLToPath(import.meta.url)), '../..');
const lib = (path) => pathToFileURL(join(root, path)).href;
const fixture = (path) => readFileSync(join(root, path), 'utf8');

const failures = [];
async function check(name, fn) {
  try {
    await fn();
    console.log(`ok   ${name}`);
  } catch (e) {
    failures.push(name);
    console.error(`FAIL ${name}\n     ${e?.stack ?? e}`);
  }
}

// Node needs two shims for browser-oriented modules: `window`, and image/style imports as URL strings.
globalThis.window ??= globalThis;
registerHooks({
  load(url, context, nextLoad) {
    if (/\.(jpe?g|png|gif|svg|webp|css|scss)$/.test(url)) {
      return { format: 'module', source: `export default ${JSON.stringify(url)};`, shortCircuit: true };
    }
    return nextLoad(url, context);
  },
});

const [{ DefaultPluginSpec }, { PluginContext }, { createViewerSpec }, { ViewerAutoPreset }, { PluginUIContext }] =
  await Promise.all([
    import(lib('packages/plugin/core/lib/default-spec.js')),
    import(lib('packages/plugin/core/lib/context.js')),
    import(lib('apps/viewer/lib/plugin-spec.js')),
    import(lib('apps/viewer/lib/presets.js')),
    import(lib('packages/plugin/ui/lib/context.js')),
  ]);
const loaders = await import(lib('extensions/plugin/lib/loaders.js'));
const { Script } = await import(lib('packages/model/lib/script/script.js'));
const { StructureElement } = await import(lib('packages/model/lib/model/structure.js'));
const { StateTransform } = await import(lib('packages/core/lib/state/transform.js'));
const { StructureRepresentation3D } = await import(
  lib('packages/plugin/core/lib/state/transforms/structure/representation.js')
);
const { DefaultFormats } = await import(lib('packages/plugin/core/lib/default-registry.js'));

const crambin = fixture('data/examples/1crn.cif');
const tinyPdb = fixture('smoke/fixtures/tiny.pdb');

async function createDefaultPlugin() {
  const plugin = new PluginContext(DefaultPluginSpec());
  await plugin.init();
  return plugin;
}

/** A plugin with the Viewer's spec, set up as `Viewer.create` does before the UI renders. */
async function createViewerPlugin(options = {}) {
  const plugin = new PluginUIContext(createViewerSpec(options));
  await plugin.init();
  plugin.builders.structure.representation.registerPreset(ViewerAutoPreset);
  return plugin;
}

async function loadCrambin(plugin) {
  const data = await plugin.builders.data.rawData({ data: crambin });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
  await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
}

const takeSnapshot = (plugin) =>
  JSON.parse(
    JSON.stringify(
      plugin.state.getSnapshot({
        camera: false,
        canvas3d: false,
        canvas3dContext: false,
        animation: false,
        startAnimation: false,
      }),
    ),
  );
const cells = (plugin) =>
  [...plugin.state.data.cells.values()].filter((c) => c.transform.ref !== StateTransform.RootRef);
const refs = (plugin) =>
  cells(plugin)
    .map((c) => c.transform.ref)
    .sort();
const reprTypes = (plugin) =>
  cells(plugin)
    .filter((c) => c.transform.transformer === StructureRepresentation3D)
    .map((c) => c.transform.params.type.name)
    .sort();

await check('a default-plugin snapshot restores in a plugin with the Viewer spec', async () => {
  const source = await createDefaultPlugin();
  await loadCrambin(source);
  const snapshot = takeSnapshot(source);
  const expectedRefs = refs(source);
  const expectedTypes = reprTypes(source);
  source.dispose();
  assert.ok(expectedTypes.includes('cartoon'));

  const viewer = await createViewerPlugin();
  await viewer.state.setSnapshot(snapshot);
  assert.deepEqual(refs(viewer), expectedRefs);
  assert.deepEqual(reprTypes(viewer), expectedTypes);
  assert.ok(cells(viewer).every((c) => c.status === 'ok'));
  assert.equal(viewer.managers.structure.hierarchy.current.structures.length, 1);

  // a Viewer snapshot restores in a Viewer as well
  const again = takeSnapshot(viewer);
  viewer.dispose();
  const other = await createViewerPlugin();
  await other.state.setSnapshot(again);
  assert.deepEqual(refs(other), expectedRefs);
  other.dispose();
});

await check('a snapshot needing an unregistered transformer fails before changing the Viewer plugin', async () => {
  const source = await createDefaultPlugin();
  await loadCrambin(source);
  const broken = takeSnapshot(source);
  source.dispose();
  const transforms = broken.data.tree.transforms;
  const last = transforms[transforms.length - 1];
  (Array.isArray(last) ? last[0] : last).transformer = 'ms-plugin.not-imported-in-this-plugin';

  const viewer = await createViewerPlugin();
  await loadCrambin(viewer);
  const before = refs(viewer);
  assert.ok(before.length > 0);
  await assert.rejects(
    viewer.state.setSnapshot(broken),
    /Snapshot uses transformers that are not available in this plugin: ms-plugin\.not-imported-in-this-plugin/,
  );
  assert.deepEqual(refs(viewer), before);
  assert.equal(viewer.managers.structure.hierarchy.current.structures.length, 1);
  viewer.dispose();
});

await check(
  'Viewer customFormats overrides a built-in format and loads the data with the custom provider',
  async () => {
    const pdb = DefaultFormats.formats.find((f) => f.name === 'pdb');
    let parsed = 0;
    const custom = {
      ...pdb,
      name: undefined,
      label: 'Custom pdb',
      parse: (...args) => {
        parsed++;
        return pdb.parse(...args);
      },
    };
    const viewer = await createViewerPlugin({ customFormats: [['pdb', custom]] });
    assert.equal(viewer.dataFormats.get('pdb').label, 'Custom pdb');
    assert.equal(viewer.dataFormats.list.filter((f) => f.name === 'pdb').length, 1);
    await loaders.loadStructureFromData(viewer, tinyPdb, 'pdb');
    assert.equal(parsed, 1);
    assert.equal(viewer.managers.structure.hierarchy.current.structures.length, 1);
    viewer.dispose();
  },
);

await check("Viewer loadTrajectory({ preset: 'all-models' }) works", async () => {
  const lammpstrj = [];
  for (let f = 0; f < 2; f++) {
    lammpstrj.push('ITEM: TIMESTEP', `${f}`, 'ITEM: NUMBER OF ATOMS', '3', 'ITEM: BOX BOUNDS pp pp pp');
    lammpstrj.push('0 10', '0 10', '0 10', 'ITEM: ATOMS id type x y z');
    for (let i = 0; i < 3; i++) lammpstrj.push(`${i + 1} 1 ${i * 1.45 + f} 0 0`);
  }
  const viewer = await createViewerPlugin();
  const { preset } = await loaders.loadTrajectory(viewer, {
    model: { kind: 'model-data', data: tinyPdb, format: 'pdb' },
    coordinates: { kind: 'coordinates-data', data: lammpstrj.join('\n') + '\n', format: 'lammpstrj' },
    preset: 'all-models',
  });
  assert.equal(preset.structures.length, 2);
  viewer.dispose();
});

await check('the Viewer spec evaluates a PyMOL script', async () => {
  const viewer = await createViewerPlugin();
  const data = await viewer.builders.data.rawData({ data: crambin });
  const trajectory = await viewer.builders.structure.parseTrajectory(data, 'mmcif');
  const model = await viewer.builders.structure.createModel(trajectory);
  const structure = (await viewer.builders.structure.createStructure(model)).obj.data;
  const loci = Script.toLoci({ language: 'pymol', expression: 'resn ALA' }, structure);
  assert.ok(StructureElement.Loci.size(loci) > 0);
  viewer.dispose();
});

function runCli(script, args, cwd = root) {
  return spawnSync(process.execPath, [join(root, script), ...args], { cwd, encoding: 'utf8' });
}

await check('cli/mvs-render prints its usage', () => {
  const r = runCli('cli/mvs-render/lib/mvs-render.js', ['--help']);
  assert.equal(r.status, 0, r.stderr);
  assert.match(r.stdout, /usage: mvs-render/);
  assert.match(r.stdout, /--input/);
});

await check('cli/state-docs writes the reference listing every registered representation and theme', async () => {
  const dir = mkdtempSync(join(tmpdir(), 'molstar-state-docs-'));
  try {
    mkdirSync(join(dir, 'docs/state'), { recursive: true });
    const r = runCli('cli/state-docs/lib/index.js', [], dir);
    assert.equal(r.status, 0, r.stderr);
    const md = readFileSync(join(dir, 'docs/state/transforms.md'), 'utf8');
    assert.ok((md.match(/^## <a name=/gm) ?? []).length > 100);
    const plugin = await createDefaultPlugin();
    for (const scope of ['structure', 'volume', 'particles']) {
      const { registry, themes } = plugin.representation[scope];
      const names = [
        ...registry.list.map((p) => p.name),
        ...themes.colorThemeRegistry.list.map((p) => p.name),
        ...themes.sizeThemeRegistry.list.map((p) => p.name),
      ];
      for (const name of names) assert.ok(md.includes(`**${name}**`), `${scope} '${name}' missing from the docs`);
    }
    plugin.dispose();
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
});

if (failures.length > 0) {
  console.error(`\n${failures.length} composition check(s) failed.`);
  process.exit(1);
}
console.log('\nAll composition checks passed.');
// Some optional dependencies keep handles open.
process.exit(0);
