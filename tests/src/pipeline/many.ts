/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Headless bench of one scene of N separate structures (many-fixtures.ts: 1.2k–11.5k atoms each, nothing shared),
// through a real PluginContext with the default spec: the state tree, the default hierarchy and representation
// presets, the structure hierarchy manager and the interactivity behaviors. Here the costs are fixed per object
// (state cell, structure, component, representation, visual), not per atom. Results: .v6/perf-audit/14-many-structures.md.
//
// There is no WebGL in Node: `plugin.canvas3d` is a stub that does what Canvas3D does with representations (add,
// remove, mark every representation, update) but has no scene; it counts the calls that reach the scene, whose cost
// (frozen gl internals) is measured in the browser. The graphics context is the CPU stub of bench.ts (impostors on).
//
// usage (through scripts/pipeline-bench.mjs, which bundles this file):
//   node scripts/pipeline-bench.mjs many-fixtures --count 500                     # once
//   node scripts/pipeline-bench.mjs many --count 300 [--runs 2]                   # stage table
//   node scripts/pipeline-bench.mjs many --count 300 --profile <dir> [--alloc]    # one profile per stage
//   node scripts/pipeline-bench.mjs analyze <dir>
//
// options: --count N (structures)  --transaction none|each|stage (default none: every commit notifies the managers)
//          --only a,b / --skip a,b (as bench.ts)  --preset <representation preset> (auto)
//          --represet <representation preset> (polymer-cartoon: the `represet` stage applies it to every structure, as
//          the UI's preset menu does; each structure's commit removes its old components and representations)
//          --mesh (no impostors)  --out <file.json>  --quality <visual quality>
//
// Stages over all N structures; the per-item column divides by N (or by the number of interactive events). Single
// stages (`*:one`) act on one structure in the full scene: what should not depend on N.

import { mkdirSync, readFileSync, existsSync } from 'node:fs';
import path from 'node:path';
import { OrderedSet } from '@molstar/core/data/int';
import { type ButtonsType, ModifiersKeys } from '@molstar/core/util/input/input-observer';
import { Color } from '@molstar/core/util/color';
import { EmptyLoci } from '@molstar/model/model/loci';
import { type Structure, StructureElement, Unit } from '@molstar/model/model/structure';
import { SecondaryStructureProvider } from '@molstar/model/props/computed/secondary-structure';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { setStructureOverpaint, clearStructureOverpaint } from '@molstar/plugin/state/helpers/structure-overpaint';
import type { StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import { graphicsStub, stubCanvas } from './canvas-stub.js';
import { ROOT } from './common.js';
import { createStageRunner, exitOnFailures, type StageResult, summarize } from './stage.js';

const args = process.argv.slice(2);
const arg = (n: string, d?: string) => {
  const i = args.indexOf(`--${n}`);
  return i >= 0 ? args[i + 1] : d;
};
const list = (n: string, d = '') => (arg(n, d) || '').split(',').filter(Boolean);
const count = +arg('count', '100')!;
const runs = +arg('runs', '1')!;
const profileDir = arg('profile');
const allocMode = args.includes('--alloc');
const meshMode = args.includes('--mesh');
const preset = arg('preset', 'auto')!;
const represet = arg('represet', 'polymer-cartoon')!;
/** none: every commit notifies (plugin API scripts); each: one data transaction per structure; stage: one per stage (Prism's load.files) */
const transaction = arg('transaction', 'none') as 'none' | 'each' | 'stage';
const quality = arg('quality');
const fixtureDir = arg('fixtures', path.join(ROOT, 'tests/build/pipeline/fixtures/many'))!;
const runner = createStageRunner({ only: list('only'), skip: list('skip'), profileDir, alloc: allocMode });
const { stage } = runner;

const files = Array.from({ length: count }, (_, i) => path.join(fixtureDir, `${String(i).padStart(3, '0')}.bcif`));
if (!existsSync(files[count - 1]))
  throw new Error(`missing ${files[count - 1]}: run 'pipeline-bench.mjs many-fixtures'`);
const data = files.map((f) => new Uint8Array(readFileSync(f)));

const graphics = graphicsStub(meshMode);

type Loaded = { data: any; trajectory: any; structureProperties?: any };

const hover = (plugin: PluginContext, loci: any) =>
  plugin.behaviors.interaction.hover.next({
    current: { loci },
    buttons: 0 as ButtonsType,
    button: 0 as ButtonsType.Flag,
    modifiers: ModifiersKeys.None,
  });
/** what SelectLoci does on a click (`clickToggle`); the click event itself also focuses the camera, which needs a scene */
const click = (plugin: PluginContext, loci: any) => plugin.managers.interactivity.lociSelects.toggle({ loci });

/** a residue of the first polymer unit of a structure */
function residueLoci(s: Structure, k = 0) {
  const unit = s.units.find((u) => Unit.isAtomic(u) && u.polymerElements.length > 0) ?? s.units[0];
  const { index, offsets } = unit.model.atomicHierarchy.residueAtomSegments;
  const rI = index[unit.elements[0]] + k;
  const start = OrderedSet.findPredecessorIndex(unit.elements as any, offsets[rI]);
  const end = OrderedSet.findPredecessorIndex(unit.elements as any, offsets[rI + 1]);
  return StructureElement.Loci(s, [{ unit, indices: OrderedSet.ofBounds(start, Math.max(start + 1, end)) as any }]);
}

function sceneInfo(c: ReturnType<typeof stubCanvas>, plugin: PluginContext) {
  let visuals = 0;
  for (const r of c.reprs.keys()) visuals += (r as any).visualCount ?? 0;
  return `cells=${plugin.state.data.cells.size} reprs=${c.reprs.size} ROs=${c.renderObjects}`;
}

async function pipeline(results: Record<string, StageResult>) {
  const plugin = new PluginContext(DefaultPluginSpec());
  await plugin.init();
  const c = stubCanvas(graphics);
  (plugin as any).canvas3d = c.canvas;
  const builder = plugin.builders.structure;
  const loaded: Loaded[] = [];
  const reprParams = quality ? { quality } : undefined;
  const inStage = (f: () => Promise<void>) => (transaction === 'stage' ? plugin.dataTransaction(f) : f());
  const inEach = (f: () => Promise<void>) => (transaction === 'each' ? plugin.dataTransaction(f) : f());

  await stage(
    'load',
    results,
    async () => {
      await inStage(async () => {
        for (let i = 0; i < count; i++) {
          await inEach(async () => {
            const d = await plugin.builders.data.rawData({ data: data[i], label: path.basename(files[i]) });
            const trajectory = await builder.parseTrajectory(d, 'mmcif');
            loaded.push({ data: d, trajectory });
          });
        }
      });
      return `cells=${plugin.state.data.cells.size}`;
    },
    { required: true, count },
  );
  // the default hierarchy preset (hierarchy-presets/default.ts) without its representation preset
  await stage(
    'model+structure',
    results,
    async () => {
      await inStage(async () => {
        for (const l of loaded) {
          await inEach(async () => {
            const model = await builder.createModel(l.trajectory);
            const modelProperties = await builder.insertModelProperties(model);
            const structure = await builder.createStructure(modelProperties || model);
            l.structureProperties = await builder.insertStructureProperties(structure);
            await builder.tryCreateUnitcell(modelProperties, undefined, { isHidden: true });
          });
        }
      });
      const atoms = loaded.reduce((a, l) => a + l.structureProperties.data.elementCount, 0);
      return `atoms=${atoms} cells=${plugin.state.data.cells.size}`;
    },
    { required: true, count },
  );
  const structures = () => loaded.map((l) => l.structureProperties.data as Structure);
  await stage(
    'bonds',
    results,
    () => {
      let intra = 0;
      let inter = 0;
      for (const s of structures()) {
        for (const u of s.units) if (Unit.isAtomic(u)) intra += u.bonds.edgeCount;
        inter += s.interUnitBonds.edgeCount;
      }
      return `intra=${intra} inter=${inter}`;
    },
    { required: true, count },
  );
  await stage(
    'secondary-structure',
    results,
    async () => {
      const ctx = {
        runtime: (await import('@molstar/core/task')).RuntimeContext.Synchronous,
        assetManager: plugin.managers.asset,
      };
      for (const s of structures()) await SecondaryStructureProvider.attach(ctx, s, void 0, true);
    },
    { required: true, count },
  );
  await stage(
    `repr:${preset}`,
    results,
    async () => {
      await inStage(async () => {
        for (const l of loaded)
          await inEach(async () => {
            await plugin.builders.structure.representation.applyPreset(
              l.structureProperties,
              preset as any,
              reprParams as any,
            );
          });
      });
      return sceneInfo(c, plugin);
    },
    { required: true, count },
  );

  const hierarchy = () => plugin.managers.structure.hierarchy.current.structures;
  const allComponents = () => hierarchy().flatMap((s) => s.components);
  const mid = (xs: StructureRef[]) => xs[Math.floor(xs.length / 2)];

  await stage(
    'theme:all',
    results,
    () => plugin.managers.structure.component.updateRepresentationsTheme(allComponents(), { color: 'chain-id' }),
    { count },
  );
  await stage('theme:one', results, () =>
    plugin.managers.structure.component.updateRepresentationsTheme(mid(hierarchy()).components, {
      color: 'element-symbol',
    }),
  );

  // interactive: events on the middle structure, through the behaviors (labels, marks on every representation)
  const target = mid(hierarchy()).cell.obj!.data as Structure;
  const hovers = 20;
  const before = { ...c.counts };
  await stage(
    'hover',
    results,
    () => {
      for (let k = 0; k < hovers; k++) hover(plugin, residueLoci(target, k));
      hover(plugin, EmptyLoci);
      return `marks=${c.counts.mark - before.mark} sceneUpdate(ROs)=${c.counts.sceneUpdateAll - before.sceneUpdateAll}`;
    },
    { count: hovers + 1 },
  );
  const clicks = 10;
  await stage(
    'select',
    results,
    () => {
      for (let k = 0; k < clicks; k++) click(plugin, residueLoci(target, k));
    },
    { count: clicks },
  );
  await stage('deselect-all', results, () => plugin.managers.interactivity.lociSelects.deselectAll());

  const overpaintOne = mid(hierarchy()).components;
  await stage('overpaint:one', results, () =>
    setStructureOverpaint(plugin, overpaintOne, Color(0xff0000), async (s) => residueLoci(s, 0)),
  );
  await stage('overpaint-clear:one', results, () => clearStructureOverpaint(plugin, overpaintOne));
  await stage(
    'overpaint:all',
    results,
    () => setStructureOverpaint(plugin, allComponents(), Color(0xff0000), async (s) => residueLoci(s, 0)),
    { count },
  );
  await stage('overpaint-clear:all', results, () => clearStructureOverpaint(plugin, allComponents()), { count });

  await stage('hide:one', results, () =>
    plugin.managers.structure.hierarchy.toggleVisibility([mid(hierarchy())], 'hide'),
  );
  await stage('show:one', results, () =>
    plugin.managers.structure.hierarchy.toggleVisibility([mid(hierarchy())], 'show'),
  );
  await stage('hide:all', results, () => plugin.managers.structure.hierarchy.toggleVisibility(hierarchy(), 'hide'), {
    count,
  });
  await stage('show:all', results, () => plugin.managers.structure.hierarchy.toggleVisibility(hierarchy(), 'show'), {
    count,
  });

  await stage(
    'represet',
    results,
    async () => {
      const provider = plugin.builders.structure.representation.resolveProvider(represet as any);
      if (!provider) throw new Error(`unknown representation preset '${represet}'`);
      await plugin.managers.structure.component.applyPreset(hierarchy(), provider);
      return sceneInfo(c, plugin);
    },
    { count },
  );

  await stage('remove:one', results, async () => {
    const one = loaded.splice(Math.floor(loaded.length / 2), 1)[0];
    await plugin.managers.structure.hierarchy.remove([one.data.ref]);
    return sceneInfo(c, plugin);
  });
  await stage(
    'remove:all',
    results,
    async () => {
      await plugin.managers.structure.hierarchy.remove(loaded.map((l) => l.data.ref));
      return sceneInfo(c, plugin);
    },
    { count: count - 1 },
  );
  plugin.dispose();
}

if (profileDir) {
  mkdirSync(profileDir, { recursive: true });
  if (!args.includes('--no-warmup')) await pipeline({});
  await runner.startProfiling();
  const res: Record<string, StageResult> = {};
  await pipeline(res);
  runner.writeProfiles(res);
} else {
  const all: Record<string, StageResult[]> = {};
  for (let i = 0; i < runs; ++i) {
    if (runs > 1) console.error(`# run ${i + 1}`);
    const res: Record<string, StageResult> = {};
    await pipeline(res);
    for (const [k, v] of Object.entries(res)) (all[k] ??= []).push(v);
  }
  const out = {
    bench: 'many',
    count,
    preset,
    transaction,
    mesh: meshMode,
    quality: quality ?? 'auto',
    runs,
    rows: summarize(all),
  };
  const outFile = arg('out');
  if (outFile) (await import('node:fs')).writeFileSync(outFile, JSON.stringify(out, null, 1));
  else console.log(JSON.stringify(out));
}
exitOnFailures(runner.failureCount());
