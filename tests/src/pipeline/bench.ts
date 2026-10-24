/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Headless CPU and memory profile of the data pipeline: parse → model → structure → bonds → props → representation
// geometry → themes, and the interactive paths on one chain (hover and select marks, overpaint, transparency,
// components). Method and results: .v6/perf-audit/00-profile.md, 13-large-scale.md.
//
// usage (through scripts/pipeline-bench.mjs, which bundles this file):
//   node scripts/pipeline-bench.mjs run --fixture asm-192 --runs 2             # stage table, median of runs
//   node scripts/pipeline-bench.mjs run --fixture asm-192 --profile <dir>      # one .cpuprofile per stage
//   node scripts/pipeline-bench.mjs analyze <dir>                              # hotspots per stage
//
// options: --only a,b / --skip a,b (stage names or prefixes ending in ':' or '*', e.g. 'repr:molecular-surface*')
//          --reprs cartoon,ball-and-stick,...  --mesh (no impostors)  --out <file.json>  --no-warmup (with --profile)
//          --alloc (with --profile: sampled JS-heap allocations per stage, .heapprofile, instead of CPU)
//          --quality <visual quality> (default: the representations' 'auto', which drops to 'lowest' from 1M atoms)
//          --heap <MB> (runner: V8 old space, default 14000)
//
// Time: wall ms per stage; GC: pause time inside it. Alloc: V8 total_allocated_bytes delta (JS heap, Node 25+;
// typed-array stores not counted). Retained: heapUsed + arrayBuffers after a full GC minus before, stage outputs
// alive; negative = given back. RSS: peak resident set size so far. Each representation is destroyed in its
// `release:<repr>` stage after its theme and interactive stages, so the peak holds one representation at a time.

import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { Session } from 'node:inspector/promises';
import path from 'node:path';
import { PerformanceObserver } from 'node:perf_hooks';
import v8 from 'node:v8';
import { Interval } from '@molstar/core/data/int';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { RuntimeContext } from '@molstar/core/task';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { Color } from '@molstar/core/util/color';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { InteractionsRepresentationProvider } from '@molstar/graphics/props/computed/representations/interactions';
import { AccessibleSurfaceAreaColorThemeProvider } from '@molstar/graphics/props/computed/themes/accessible-surface-area';
import { InteractionTypeColorThemeProvider } from '@molstar/graphics/props/computed/themes/interaction-type';
import type { Representation } from '@molstar/graphics/repr/representation';
import { BallAndStickRepresentationProvider } from '@molstar/graphics/repr/structure/representation/ball-and-stick';
import { CartoonRepresentationProvider } from '@molstar/graphics/repr/structure/representation/cartoon';
import { GaussianSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-surface';
import { MolecularSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/molecular-surface';
import { SpacefillRepresentationProvider } from '@molstar/graphics/repr/structure/representation/spacefill';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { BuiltInColorThemes } from '@molstar/graphics/theme/color/catalog';
import { Overpaint } from '@molstar/graphics/theme/overpaint';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { BuiltInSizeThemes } from '@molstar/graphics/theme/size/catalog';
import { Theme } from '@molstar/graphics/theme/theme';
import { Transparency } from '@molstar/graphics/theme/transparency';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { EveryLoci } from '@molstar/model/model/loci';
import { Structure, StructureElement, StructureProperties, Unit } from '@molstar/model/model/structure';
import { Queries, QueryContext, StructureSelection } from '@molstar/model/model/structure/query';
import { StructureSymmetry } from '@molstar/model/model/structure/structure/symmetry';
import type { CustomProperty } from '@molstar/model/props/common/custom-property';
import { AccessibleSurfaceAreaProvider } from '@molstar/model/props/computed/accessible-surface-area';
import { InteractionsProvider } from '@molstar/model/props/computed/interactions';
import { SecondaryStructureProvider } from '@molstar/model/props/computed/secondary-structure';
import { ValenceModelProvider } from '@molstar/model/props/computed/valence-model';
import { ROOT, SOURCE_1TUP } from './common.js';

type Fixture = { file: string; symmetry?: [number, number, number, number, number, number]; assembly?: string };
const generated = (file: string) => path.join(ROOT, 'tests/build/pipeline/fixtures', file);
const FIXTURES: Record<string, Fixture> = {
  '7qpd': { file: path.join(ROOT, 'data/examples/7qpd.fw2.cif') },
  '1tup-2x2x2': { file: path.join(ROOT, SOURCE_1TUP), symmetry: [0, 0, 0, 1, 1, 1] },
  // made by fixtures.ts: N copies of 1tup (5,828 atoms, 9 units) as an assembly of N operators ('asm', 9 shared
  // symmetry groups) or written as one model ('flat', one group per unit). flat-32 is 00-profile's big-bcif.
  'asm-32': { file: generated('1tup-asm-32.cif'), assembly: '1' },
  'asm-192': { file: generated('1tup-asm-192.cif'), assembly: '1' },
  'asm-864': { file: generated('1tup-asm-864.cif'), assembly: '1' },
  'flat-32': { file: generated('1tup-flat-32.bcif') },
  'flat-192': { file: generated('1tup-flat-192.bcif') },
  'flat-192-cif': { file: generated('1tup-flat-192.cif') },
  'flat-864': { file: generated('1tup-flat-864.bcif') },
};

const args = process.argv.slice(2);
const arg = (n: string, d?: string) => {
  const i = args.indexOf(`--${n}`);
  return i >= 0 ? args[i + 1] : d;
};
const list = (n: string, d = '') => (arg(n, d) || '').split(',').filter(Boolean);
const fixtureName = arg('fixture', '7qpd')!;
const fixture = FIXTURES[fixtureName] ?? { file: path.resolve(fixtureName), assembly: arg('assembly') };
const runs = +arg('runs', '1')!;
const profileDir = arg('profile');
const allocMode = args.includes('--alloc');
const only = list('only');
const skip = list('skip');
const meshMode = args.includes('--mesh');
const quality = arg('quality');
const reprNames = list('reprs', 'cartoon,ball-and-stick,spacefill,gaussian-surface,molecular-surface,interactions');

/** the parts of a GraphicsContext the CPU geometry paths read: impostors on, no GPU compute (tryUseGpu: false) */
const graphicsStub = meshMode
  ? undefined
  : ({
      capabilities: { fragDepth: true, flatInterpolation: true, texture: true },
      resources: {
        limits: { maxTextureSize: 16384, max3dTextureSize: 2048 },
        // a texture that ignores uploads (the CPU work before the upload is what is measured)
        texture: () => new Proxy({} as any, { get: (t, k) => (k === 'then' ? undefined : k in t ? t[k] : () => {}) }),
        ownsTexture: () => true,
      },
      namedTextures: {},
    } as any);

const v8gc = (globalThis as any).gc as (options?: object) => void;
if (!v8gc) throw new Error('run with --expose-gc');
/**
 * A full collection whose freed array-buffer backing stores are also accounted: they are swept concurrently, and
 * `memoryUsage().arrayBuffers` drops only after the sweep, so without the wait the next stage gets the credit
 */
async function gc() {
  v8gc({ type: 'major', execution: 'sync', flavor: 'last-resort' });
  await new Promise((r) => setTimeout(r, 20));
  v8gc({ type: 'major', execution: 'sync', flavor: 'last-resort' });
}

const registries = (() => {
  const colorThemeRegistry = ColorTheme.createRegistry();
  for (const p of Object.values(BuiltInColorThemes)) colorThemeRegistry.add(p as ColorTheme.Provider<any, any, any>);
  colorThemeRegistry.add(AccessibleSurfaceAreaColorThemeProvider as any);
  colorThemeRegistry.add(InteractionTypeColorThemeProvider as any);
  const sizeThemeRegistry = SizeTheme.createRegistry();
  for (const p of Object.values(BuiltInSizeThemes)) sizeThemeRegistry.add(p as SizeTheme.Provider<any, any>);
  return { colorThemeRegistry, sizeThemeRegistry };
})();
const propCtx: CustomProperty.Context = {
  runtime: RuntimeContext.Synchronous,
  assetManager: new DefaultAssetManager(),
};

// cumulative JS-heap allocation: `total_allocated_bytes` is in Node 25+; without it the Alloc column is NaN
const hasAllocCounter = typeof (v8.getHeapStatistics() as any).total_allocated_bytes === 'number';
if (!hasAllocCounter) console.error('warning: Node 25+ is needed for the Alloc column (v8 total_allocated_bytes)');
function mem() {
  const m = process.memoryUsage();
  return {
    heap: m.heapUsed,
    ab: m.arrayBuffers,
    alloc: hasAllocCounter ? ((v8.getHeapStatistics() as any).total_allocated_bytes as number) : NaN,
  };
}

type StageResult = {
  ms: number;
  allocMB: number;
  retainedHeapMB: number;
  retainedAbMB: number;
  rssMB: number;
  gcMs: number;
  info?: string;
  /** the stage threw: its columns measure the failure, not the work (perf-audit-3 12/B4) */
  error?: string;
};
type Ctx = { parsed?: any; model?: any; structure?: Structure };

/** GC pauses (start, duration); a stage sums those that start inside it */
const gcPauses: [number, number][] = [];
new PerformanceObserver((list) => {
  for (const e of list.getEntries()) gcPauses.push([e.startTime, e.duration]);
}).observe({ entryTypes: ['gc'] });
const gcTime = (t0: number, t1: number) => gcPauses.reduce((a, [t, d]) => (t >= t0 && t <= t1 ? a + d : a), 0);

let session: Session | undefined;
const stageProfiles: Record<string, any> = {};
/** stages that threw: the run exits non-zero */
let failures = 0;

const matches = (name: string, pattern: string) =>
  pattern.endsWith('*') || pattern.endsWith(':') ? name.startsWith(pattern.replace(/\*$/, '')) : name === pattern;
const enabled = (name: string) =>
  !skip.some((p) => matches(name, p)) && (!only.length || only.some((p) => matches(name, p)));

/**
 * A `required` stage runs unmeasured when --only/--skip leave it out, so that later stages do not pay for it: parse,
 * model, structure, lookup, bonds, secondary structure and the component producers. Valence model, interactions and
 * ASA are not required: a representation that needs them (interactions) computes them inside its own stage.
 */
async function stage(
  name: string,
  results: Record<string, StageResult>,
  fn: () => Promise<string | void> | string | void,
  required = false,
) {
  if (!enabled(name)) {
    if (required) await fn();
    return;
  }
  await gc();
  const before = mem();
  if (session && allocMode) {
    // every sampled allocation, also those a later GC collects: what the stage allocates, not what it keeps
    await session.post('HeapProfiler.startSampling', {
      samplingInterval: 64 * 1024,
      includeObjectsCollectedByMajorGC: true,
      includeObjectsCollectedByMinorGC: true,
    });
  } else if (session) await session.post('Profiler.start');
  const t0 = performance.now();
  let info: string | void = undefined;
  let error: string | undefined;
  try {
    info = await fn();
  } catch (e) {
    error = (e as Error)?.message ?? String(e);
    failures++;
    console.error(e);
  }
  const t1 = performance.now();
  await new Promise((r) => setTimeout(r, 0)); // gc entries are delivered asynchronously
  const stageGcMs = gcTime(t0, t1);
  gcPauses.length = 0;
  if (session) {
    const { profile } = allocMode
      ? await session.post('HeapProfiler.stopSampling')
      : await session.post('Profiler.stop');
    stageProfiles[name] = profile;
  }
  const mid = mem();
  await gc();
  const after = mem();
  results[name] = {
    ms: t1 - t0,
    allocMB: (mid.alloc - before.alloc) / 2 ** 20,
    retainedHeapMB: (after.heap - before.heap) / 2 ** 20,
    retainedAbMB: (after.ab - before.ab) / 2 ** 20,
    rssMB: process.resourceUsage().maxRSS / 1024,
    gcMs: stageGcMs,
    info: info || undefined,
    error,
  };
  if (!profileDir) {
    const r = results[name];
    if (r.error) {
      console.error(`${name.padEnd(44)} FAILED: ${r.error}`);
      return;
    }
    console.error(
      `${name.padEnd(44)} ${r.ms.toFixed(1).padStart(9)} ms ${r.allocMB.toFixed(0).padStart(6)} MB alloc ${(r.retainedHeapMB + r.retainedAbMB).toFixed(0).padStart(6)} MB ret ${r.rssMB.toFixed(0).padStart(6)} MB rss ${r.gcMs.toFixed(0).padStart(7)} ms gc  ${r.info ?? ''}`,
    );
  }
}

/** sum of typed-array bytes of a representation's render-object values, and the top buffers by size */
function reprBuffers(repr: Representation.Any) {
  const byName = new Map<string, number>();
  let total = 0;
  let instances = 0;
  const seen = new Set<ArrayBufferLike>();
  for (const ro of repr.renderObjects) {
    const values = ro.values as Record<string, any>;
    instances += values.instanceCount?.ref.value ?? 0;
    for (const [k, cell] of Object.entries(values)) {
      const v = cell?.ref?.value;
      let arr: ArrayBufferView | undefined;
      if (ArrayBuffer.isView(v)) arr = v;
      else if (v && ArrayBuffer.isView(v.array)) arr = v.array; // TextureImage
      if (!arr || seen.has(arr.buffer)) continue;
      seen.add(arr.buffer);
      total += arr.byteLength;
      byName.set(k, (byName.get(k) ?? 0) + arr.byteLength);
    }
  }
  const top = [...byName.entries()]
    .sort((a, b) => b[1] - a[1])
    .slice(0, 4)
    .map(([k, b]) => `${k}=${(b / 2 ** 20).toFixed(1)}`)
    .join(' ');
  return `ROs=${repr.renderObjects.length} inst=${instances} buf=${(total / 2 ** 20).toFixed(1)}MB [${top}]`;
}

const reprProviders = new WeakMap<object, any>();
async function makeRepr(provider: any, structure: Structure) {
  await provider.ensureCustomProperties?.attach(propCtx, structure);
  const themeProps = {
    colorTheme: { name: provider.defaultColorTheme.name, params: provider.defaultColorTheme.props ?? {} },
    sizeTheme: { name: provider.defaultSizeTheme.name, params: provider.defaultSizeTheme.props ?? {} },
  };
  await Theme.ensureDependencies(propCtx, registries, { structure }, themeProps);
  const repr = provider.factory({ graphics: graphicsStub, ...registries }, provider.getParams) as Representation.Any;
  repr.setTheme(Theme.create(registries, { structure }, themeProps));
  const params = PD.getDefaultValues(provider.getParams(registries, structure));
  await repr.createOrUpdate({ ...params, ...(quality ? { quality } : {}), tryUseGpu: false }, structure).run();
  reprProviders.set(repr, provider);
  return repr;
}

/** a color theme applied to existing geometry (the size theme is kept: a size change rebuilds geometry) */
async function retheme(repr: Representation.Any, structure: Structure, colorTheme: string) {
  const provider = reprProviders.get(repr);
  const themeProps = {
    colorTheme: { name: colorTheme, params: {} },
    sizeTheme: { name: provider.defaultSizeTheme.name, params: provider.defaultSizeTheme.props ?? {} },
  };
  await Theme.ensureDependencies(propCtx, registries, { structure }, themeProps);
  repr.setTheme(Theme.create(registries, { structure }, themeProps));
  await repr.createOrUpdate({}, structure).run();
}

/** evaluate a color theme for every element of every unit (what createColors does at group granularity) */
function evalTheme(structure: Structure, name: string) {
  const theme = registries.colorThemeRegistry.create(name, { structure }, {});
  const color = (theme as any).color as (l: any, s: boolean) => number;
  const loc = StructureElement.Location.create(structure);
  let acc = 0;
  for (const unit of structure.units) {
    loc.unit = unit;
    const { elements } = unit;
    for (let i = 0, il = elements.length; i < il; ++i) {
      loc.element = elements[i];
      acc ^= color(loc, false);
    }
  }
  return acc;
}

const chainId = (s: Structure, u: Unit) =>
  StructureProperties.chain.label_asym_id(StructureElement.Location.create(s, u, u.elements[0]));
const wholeUnits = (s: Structure, units: readonly Unit[]) =>
  StructureElement.Loci(
    s,
    units.map((unit) => ({ unit, indices: Interval.ofBounds(0, unit.elements.length) })),
  );

/**
 * One chain (the first atomic polymer unit, one operator; on a coarse-only entry the first unit) and the same chain
 * in every copy
 */
function selections(s: Structure) {
  const first = s.units.find((u) => Unit.isAtomic(u) && u.polymerElements.length > 0) ?? s.units[0];
  const id = chainId(s, first);
  const sameChain = s.units.filter(
    (u) => u.model === first.model && chainId(s, u) === id && u.elements.length === first.elements.length,
  );
  return {
    chain: wholeUnits(s, [first]),
    chainEvery: wholeUnits(s, sameChain),
    label: `${id}: 1 and ${sameChain.length} units`,
  };
}

/** what a representation built on the component touches first: boundary, lookup, inter-unit bonds */
function touchDerived(c: Structure) {
  void c.boundary;
  void c.lookup3d;
  return `units=${c.units.length} groups=${c.unitSymmetryGroups.length} elements=${c.elementCount} interBonds=${c.interUnitBonds.edgeCount}`;
}

async function interactive(name: string, repr: Representation.Any, s: Structure, results: Record<string, StageResult>) {
  const { chain, chainEvery } = selections(s);
  await stage(`hover:${name}`, results, () => {
    repr.mark(chain, MarkerAction.Highlight);
    repr.mark(chain, MarkerAction.RemoveHighlight);
  });
  await stage(`select:${name}`, results, () => void repr.mark(chain, MarkerAction.Select));
  await stage(`select-every-copy:${name}`, results, () => void repr.mark(chainEvery, MarkerAction.Select));
  repr.mark(EveryLoci, MarkerAction.Deselect);
  const overpaint = Overpaint.ofBundle(
    [{ bundle: StructureElement.Bundle.fromLoci(chain), color: Color(0xff0000), clear: false }],
    s,
  );
  await stage(`overpaint:${name}`, results, () => repr.setState({ overpaint }));
  await stage(`overpaint-clear:${name}`, results, () => repr.setState({ overpaint: Overpaint.Empty }));
  const transparency = Transparency.ofBundle([{ bundle: StructureElement.Bundle.fromLoci(chain), value: 0.5 }], s);
  await stage(`transparency:${name}`, results, () => repr.setState({ transparency }));
  await stage(`transparency-clear:${name}`, results, () => repr.setState({ transparency: Transparency.Empty }));
}

const REPRS: Record<string, { provider: any; themes: string[] }> = {
  cartoon: { provider: CartoonRepresentationProvider, themes: ['secondary-structure'] },
  'ball-and-stick': { provider: BallAndStickRepresentationProvider, themes: ['chain-id'] },
  spacefill: { provider: SpacefillRepresentationProvider, themes: ['chain-id'] },
  'gaussian-surface': { provider: GaussianSurfaceRepresentationProvider, themes: ['element-symbol'] },
  'molecular-surface': { provider: MolecularSurfaceRepresentationProvider, themes: ['residue-name'] },
  interactions: { provider: InteractionsRepresentationProvider, themes: [] },
};

async function pipeline(results: Record<string, StageResult>): Promise<void> {
  const ctx: Ctx = {};
  const binary = fixture.file.endsWith('.bcif');
  const data = binary ? new Uint8Array(readFileSync(fixture.file)) : readFileSync(fixture.file, 'utf8');

  await stage(
    'parse',
    results,
    async () => {
      const r = binary ? await CIF.parseBinary(data as Uint8Array).run() : await CIF.parse(data as string).run();
      if (r.isError) throw new Error(r.message);
      ctx.parsed = r.result;
    },
    true,
  );
  await stage(
    'model',
    results,
    async () => {
      const traj = await trajectoryFromMmCIF(ctx.parsed.blocks[0], ctx.parsed).run();
      ctx.model = traj.representative;
    },
    true,
  );
  await stage(
    'structure',
    results,
    async () => {
      let s = Structure.ofModel(ctx.model);
      if (fixture.symmetry) {
        const y = fixture.symmetry;
        s = await StructureSymmetry.buildSymmetryRange(
          s,
          Vec3.create(y[0], y[1], y[2]),
          Vec3.create(y[3], y[4], y[5]),
        ).run();
      } else if (fixture.assembly) {
        s = await StructureSymmetry.buildAssembly(s, fixture.assembly).run();
      }
      ctx.structure = s;
      return `elements=${s.elementCount} units=${s.units.length} groups=${s.unitSymmetryGroups.length} polyRes=${s.polymerResidueCount}`;
    },
    true,
  );
  const s = ctx.structure;
  if (!s) return;
  await stage(
    'lookup+boundary',
    results,
    () => {
      void s.boundary;
      void s.lookup3d;
      for (const u of s.units) void u.lookup3d;
    },
    true,
  );
  await stage(
    'bonds-intra',
    results,
    () => {
      let n = 0;
      for (const u of s.units) if (Unit.isAtomic(u)) n += u.bonds.edgeCount;
      return `edges(all units)=${n}`;
    },
    true,
  );
  await stage('bonds-inter', results, () => `edges=${s.interUnitBonds.edgeCount}`, true);
  await stage('secondary-structure', results, () => SecondaryStructureProvider.attach(propCtx, s, void 0, true), true);
  await stage('valence-model', results, () => ValenceModelProvider.attach(propCtx, s, void 0, true));
  await stage('interactions', results, async () => {
    await InteractionsProvider.attach(propCtx, s, void 0, true);
    const v = InteractionsProvider.get(s).value;
    return v ? `contacts inter=${v.contacts.edgeCount}` : 'n/a';
  });
  await stage('asa', results, () => AccessibleSurfaceAreaProvider.attach(propCtx, s, void 0, true));

  for (const t of ['element-symbol', 'chain-id', 'secondary-structure', 'residue-name']) {
    await stage(`theme-eval:${t}`, results, () => void evalTheme(s, t));
  }

  // components from a selection (plugin `StructureComponent` with a bundle, or a static query), and what a
  // representation on them computes first
  const { chain, chainEvery, label } = selections(s);
  const bundleComponent = (loci: StructureElement.Loci) =>
    StructureElement.Bundle.toStructure(StructureElement.Bundle.fromLoci(loci), s);
  let c: Structure | undefined;
  await stage(
    'component:chain',
    results,
    () => {
      c = bundleComponent(chain);
      return label;
    },
    true,
  );
  await stage('component:chain-derived', results, () => (c ? touchDerived(c) : undefined));
  await stage('component:every-copy', results, () => void (c = bundleComponent(chainEvery)), true);
  await stage('component:every-copy-derived', results, () => (c ? touchDerived(c) : undefined));
  await stage(
    'component:polymer',
    results,
    () => {
      const polymer = Queries.generators.atoms({
        entityTest: (l) => StructureProperties.entity.type(l.element) === 'polymer',
      });
      c = StructureSelection.unionStructure(polymer(new QueryContext(s)));
    },
    true,
  );
  await stage('component:polymer-derived', results, () => (c ? touchDerived(c) : undefined));
  // a rigid move (`TransformStructureConformation`: superposition, Prism's structure move): the moved structure and its
  // polymer component, made again at every move, and the inter-unit bonds their representations need
  const move = Mat4.fromRotation(Mat4(), Math.PI / 6, Vec3.create(0, 0, 1));
  let moved = s;
  await stage('transform', results, () => void (moved = Structure.transform(s, move)), true);
  await stage('transform:bonds', results, () => `edges=${moved.interUnitBonds.edgeCount}`);
  await stage(
    'transform:polymer',
    results,
    () => {
      const polymer = Queries.generators.atoms({
        entityTest: (l) => StructureProperties.entity.type(l.element) === 'polymer',
      });
      c = StructureSelection.unionStructure(polymer(new QueryContext(moved)));
    },
    true,
  );
  await stage('transform:polymer-bonds', results, () => (c ? `edges=${c.interUnitBonds.edgeCount}` : undefined));
  c = undefined;

  for (const name of reprNames) {
    if (!REPRS[name]) throw new Error(`unknown representation '${name}' (${Object.keys(REPRS).join(', ')})`);
    const { provider, themes } = REPRS[name];
    if (!provider.isApplicable(s) || !enabled(`repr:${name}`)) continue;
    let repr: Representation.Any | undefined;
    await stage(`repr:${name}`, results, async () => {
      repr = await makeRepr(provider, s);
      return reprBuffers(repr);
    });
    if (!repr) continue;
    for (const t of themes) await stage(`theme-apply:${name}:${t}`, results, () => retheme(repr!, s, t));
    if (name !== 'interactions') await interactive(name, repr, s, results);
    // destroying gives the representation's memory back: negative retained, credited here and not to the next stage
    await stage(
      `release:${name}`,
      results,
      () => {
        repr!.destroy();
        repr = undefined;
      },
      true,
    );
  }
}

const median = (xs: number[]) => {
  const a = [...xs].sort((x, y) => x - y);
  return a.length % 2 ? a[(a.length - 1) / 2] : (a[a.length / 2 - 1] + a[a.length / 2]) / 2;
};

if (profileDir) {
  mkdirSync(profileDir, { recursive: true });
  // warm-up run so the profile shows steady-state code, not the first compile
  if (!args.includes('--no-warmup')) await pipeline({});
  session = new Session();
  session.connect();
  if (allocMode) await session.post('HeapProfiler.enable');
  else {
    await session.post('Profiler.enable');
    await session.post('Profiler.setSamplingInterval', { interval: 100 });
  }
  const res: Record<string, StageResult> = {};
  await pipeline(res);
  for (const [name, p] of Object.entries(stageProfiles))
    writeFileSync(
      `${profileDir}/${name.replace(/[:/]/g, '_')}.${allocMode ? 'heapprofile' : 'cpuprofile'}`,
      JSON.stringify(p),
    );
  writeFileSync(`${profileDir}/stages.json`, JSON.stringify(res, null, 1));
  console.log('profiles written to', profileDir);
} else {
  const all: Record<string, StageResult[]> = {};
  for (let i = 0; i < runs; ++i) {
    if (runs > 1) console.error(`# run ${i + 1}`);
    const res: Record<string, StageResult> = {};
    await pipeline(res);
    for (const [k, v] of Object.entries(res)) (all[k] ??= []).push(v);
  }
  // a stage that failed in any run has no numbers, so that a comparison does not read the failure as a speedup
  const rows = Object.entries(all).map(([k, v]) => {
    const error = v.find((x) => x.error)?.error;
    if (error) return { stage: k, error, coldMs: null, medMs: null, allocMB: null, info: v[0].info ?? '' };
    return {
      stage: k,
      coldMs: +v[0].ms.toFixed(1),
      medMs: +median(v.map((x) => x.ms)).toFixed(1),
      allocMB: +median(v.map((x) => x.allocMB)).toFixed(1),
      retHeapMB: +median(v.map((x) => x.retainedHeapMB)).toFixed(1),
      retAbMB: +median(v.map((x) => x.retainedAbMB)).toFixed(1),
      rssMB: +Math.max(...v.map((x) => x.rssMB)).toFixed(0),
      gcMs: +median(v.map((x) => x.gcMs)).toFixed(0),
      info: v[0].info ?? '',
    };
  });
  const out = { fixture: fixtureName, mesh: meshMode, quality: quality ?? 'auto', runs, rows };
  const outFile = arg('out');
  if (outFile) writeFileSync(outFile, JSON.stringify(out, null, 1));
  else console.log(JSON.stringify(out));
}
if (failures) {
  console.error(`${failures} stage${failures === 1 ? '' : 's'} failed`);
  process.exitCode = 1;
}
