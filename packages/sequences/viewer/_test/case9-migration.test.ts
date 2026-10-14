/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST2's migration test, requests: case 9 (PF00042's lazy 18,000 rows) through the tiles provider asks the server the tiles
// SW1 recorded of the old client (`parity/baselines/tiles.json`: per zoom the tile paths asked and every tile held with a
// hash of its bytes), at the same three zooms, after the passes settle. Without a margin the sets are equal; with the
// default margin (ST2 gave remote alignments one, as 1-D layers have) the baseline's are a subset: the extra tiles are
// the ring around the view. The pixels are `canvas/case9-pixels.test.ts`'s.

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { case9 } from '../../core/_test/format/cases.js';
import { fixtureLoaders, MolstarFixture } from '../../core/_test/format/load-fixtures.js';
import {
  createSeqViewer,
  DefaultRegistry,
  DefaultSeqViewerSpec,
  manualScheduler,
  RegionState,
  SeqViewerConfig,
} from '../src/index.js';
import { preparedRenderers } from './layout-fixtures.js';

// biome-ignore lint/suspicious/noExplicitAny: the manager's cache is read through its regions.
type J = any;
type Zoomed = Record<string, { asked: string[]; held: string[] }>;
const baseline: Zoomed = JSON.parse(readFileSync(new URL('./parity/baselines/tiles.json', import.meta.url), 'utf8'))[
  'case 9'
];
const flush = () => new Promise((ok) => setTimeout(ok, 0));

/** A tile's bytes as one number. */
function hash(t: { cls: Uint8Array; share: Uint8Array }) {
  let h = 0;
  for (let i = 0; i < t.cls.length; i++) h = (Math.imul(h, 31) + t.cls[i] * 256 + t.share[i]) | 0;
  return h;
}

/** Case 9 at the baselines' three zooms: per zoom the tiles asked since the last and every tile held. */
async function run(margin: number): Promise<Zoomed> {
  const registry = { ...DefaultRegistry, labels: undefined, renderers: undefined };
  const spec = {
    ...DefaultSeqViewerSpec(),
    behaviors: [],
    registry: [registry, { svs: MolstarFixture }, { renderers: preparedRenderers() }],
  };
  const loaders = fixtureLoaders();
  const ctx = createSeqViewer(spec as never, { loaders, scheduler: manualScheduler() });
  ctx.config.set(SeqViewerConfig.Stream.Margin, margin);
  const done = ctx.setDocument(case9());
  for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(1000, 500, 1);
  await done;
  await ctx.graph.update();
  // Through the graph's own driver (not `compile`): nothing reported, the alignment and its statistic both streamed.
  expect(ctx.graph.issues.value).toEqual([]);
  const set = ctx.graph.assembled.value?.set as J;
  const streamed = ['pf00042', 'cons'].map((k) =>
    set.annotations.concat(set.alignments).find((x: J) => set.strings.get(x.key) === k),
  );
  expect(streamed.map((x: J) => [x?.streamed ?? x?.lazy])).toEqual([[true], [true]]);
  const vp = ctx.viewports.get('msa.main') as J;
  const zooms: Record<string, () => void> = {
    opened: () => {},
    'every column, rows at 0.01': () => {
      vp.zoomTo(0, 0, 320);
      vp.setRowScale(0.01, 0);
    },
    'letters at HBA_HUMAN': () => {
      vp.setRowScale(1, 0);
      vp.zoomTo(0, 40, 90);
      vp.scrollTo(vp.live.current.scrollX, 4711 * 18, 'api');
    },
  };
  const m = ctx.graph.assembled.value?.view as J;
  let layer = '';
  for (let k = 0; m && k < m.resolved.layers.count; k++)
    if (m.set.strings.get(m.resolved.layers.type[k]) === 'rows') layer = m.set.strings.get(m.resolved.layers.ref[k]);
  const got: Zoomed = {};
  for (const [name, zoom] of Object.entries(zooms)) {
    const from = loaders.tiles.length;
    zoom();
    ctx.layout.get('msa.main');
    for (let i = 0, was = -1; i < 40; i++) {
      (ctx.streams as J).pass();
      await flush();
      const s = ctx.streams.stats.current;
      const now = s.requests * 1e6 + s.computed;
      if (now === was && s.inFlight === 0) break;
      was = now;
    }
    const asked = [...new Set(loaders.tiles.slice(from).filter((p: string) => p.startsWith('tiles/')))].sort();
    const held: string[] = [];
    for (const r of (ctx.streams.get(layer) as J)?.all ?? [])
      if (r.tile && r.level < 256) held.push(`${r.level >> 4}/${r.level & 15}/${r.ty}/${r.tx} #${hash(r.tile)}`);
    got[name] = { asked, held: held.sort() };
  }
  ctx.dispose();
  return got;
}

describe("case 9 through the tiles provider (SW1's baselines)", () => {
  it('with no margin asks the tiles the old client asked and holds the bytes it held, at each zoom', async () => {
    expect(await run(0)).toEqual(baseline);
  });

  it("draws conservation from its served statistic over the alignment's columns (blank before)", async () => {
    const registry = { ...DefaultRegistry, labels: undefined, renderers: undefined };
    const renderers = { renderers: preparedRenderers() };
    const spec = { ...DefaultSeqViewerSpec(), behaviors: [], registry: [registry, { svs: MolstarFixture }, renderers] };
    const loaders = fixtureLoaders();
    const ctx = createSeqViewer(spec as never, { loaders, scheduler: manualScheduler() });
    const done = ctx.setDocument(case9());
    for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(1000, 500, 1);
    await done;
    await ctx.graph.update();
    const m = ctx.graph.assembled.value?.view as J;
    let layer = '';
    for (let k = 0; k < m.resolved.layers.count; k++)
      if (m.set.strings.get(m.resolved.layers.type[k]) === 'values')
        layer = m.set.strings.get(m.resolved.layers.ref[k]);
    ctx.layout.get('msa.main');
    for (let i = 0; i < 20 && ctx.streams.state(layer, 0, 5) !== RegionState.loaded; i++) {
      (ctx.streams as J).pass();
      await flush();
    }
    // 320 columns over 1000 px: the data (one value a column), (c % 10) / 10 as the fixture serves it.
    expect(ctx.streams.state(layer, 0, 5)).toBe(RegionState.loaded);
    expect(loaders.tiles.filter((p) => p.startsWith('stats/conservation/0/'))).toEqual([
      'stats/conservation/0/0',
      'stats/conservation/0/1',
    ]);
    const out = new Float64Array(4);
    out[0] = 0.32;
    expect(ctx.streams.get(layer)?.values?.(0, 20, 30, out)).toBe(true);
    expect([out[0], out[1], out[3]]).toEqual([0, 0.9, 10]);
    expect(out[2]).toBeCloseTo(4.5);
    ctx.dispose();
  });

  it('with the margin asks the same tiles and the ring around the view, holding the same bytes', async () => {
    const got = await run(1);
    const extra = Object.fromEntries(
      Object.entries(got).map(([name, z]) => [name, z.asked.filter((p) => !baseline[name].asked.includes(p))]),
    );
    for (const [name, z] of Object.entries(baseline)) {
      expect(got[name].asked, name).toEqual(expect.arrayContaining(z.asked));
      expect(got[name].held, name).toEqual(expect.arrayContaining(z.held));
    }
    // The view's four tiles at 3/0 lie in rows 0 to 2,047; the ring adds the band of rows below.
    expect(extra).toEqual({
      opened: [],
      'every column, rows at 0.01': ['tiles/3/0/2/0', 'tiles/3/0/2/1'],
      'letters at HBA_HUMAN': [],
    });
  });
});
