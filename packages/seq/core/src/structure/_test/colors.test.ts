/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A layer's colors and an annotation's rows on units (_spec/README.md "Structures", Colors; Rows): the colors are the
// drawn track's (its renderer's compile: a colourless layer's categorical place, a field's categories, an open scale's
// held domain over two layers), light and dark; unobserved positions cut (EC-118), uncolored ones in no run, every copy
// colored (EC-65), data the store does not hold or refuses past its bound pending, never an empty run (EC-169).

import { DataStore } from '@molstar/seqx-core/data';
import { partOf, positionCount } from '@molstar/seqx-core/query';
import { DarkTokens, LightTokens, layerColors, openModel, viewTheme } from '@molstar/seqx-core/render';
import type { Channel, ThemeDatum } from '@molstar/seqx-core/theme';
import { describe, expect, it } from 'vitest';
import { settle } from '../../render/_test/helpers.js';
import { colorRuns } from '../colors.js';
import { unitRows } from '../rows.js';
import { chainId, Complex, EntStart, Observed, p53Spec, Up, UpLength, ValueAt } from './p53.js';

const unit = (k: number) => `${Complex}/${chainId(k)}`;
const Size = { width: 900, height: 420 };
const Domains = {
  start: [0, 63, 99, 125, 320],
  end: [61, 94, 120, 140, 360],
  category: ['activation', 'linker', 'binding', 'linker', 'binding'],
};

/** The view drawn and settled once a file (the cvd test draws it again last but one; nothing after reads colors). */
let kept: ReturnType<typeof draw> | undefined;
const drawn = () => (kept ??= draw());

/** The view drawn and settled: its renderer's colors, and each track's drawn color channel. */
async function draw() {
  const s = await settle(p53Spec({ view: true }));
  const f = await s.frame({ size: Size, settled: true });
  const layer = (t: number) => s.view.sections[0].tracks[t].layers[0];
  const channel = (t: number) => (f.layout.sections[0].tracks[t].layers[0].lay as { color: Channel }).color;
  return { ...s, layer, channel, colors: s.renderer.layerColors() };
}

/** The datum of a P04637 position as each layer reads it. */
function datumAt(p: number, kind: 'letters' | 'values' | 'ranges', letters: string, values = ValueAt): ThemeDatum {
  const k = Domains.start.findIndex((s, i) => s <= p && p < Domains.end[i]);
  return {
    value: kind === 'values' ? values(p) : Number.NaN,
    letter: kind === 'letters' ? letters[p] : '',
    start: kind === 'ranges' ? Domains.start[k] : p,
    end: kind === 'ranges' ? Domains.end[k] : p + 1,
    field: (n) =>
      n === 'value' ? values(p) : n === 'letter' ? letters[p] : n === 'category' ? Domains.category[k] : undefined,
  };
}

describe('colorRuns: the drawn track’s colors', () => {
  it('each run’s color is the drawn channel’s at its first residue, for every layer kind (EC-65, EC-118)', async () => {
    const d = await drawn();
    const letters = d.model.systems.get(Up)!.letters!;
    const kinds = ['letters', 'ranges', 'ranges', 'ranges', 'values', 'values'] as const;
    for (let t = 0; t < kinds.length; t++) {
      const r = colorRuns(d.model, d.store, d.layer(t), undefined, { colors: d.colors });
      expect([t, r.complete, r.issues]).toEqual([t, true, []]);
      expect(r.runs.length, `track ${t}`).toBeGreaterThan(0);
      for (const run of r.runs) {
        expect(run.target.parts.map((p) => p.system)).toEqual([0, 1, 2, 3].map(unit));
        const p = partOf(run.target, unit(0))!.ranges[0] + EntStart;
        const values = t === 5 ? (i: number) => 2 + i / UpLength : ValueAt;
        expect(run.color, `track ${t} at ${p}`).toBe(d.channel(t).at(datumAt(p, kinds[t], letters, values)));
      }
    }
  });

  it('a colourless layer takes its place in the categorical palette; an open scale its held domain, not its extent', async () => {
    const d = await drawn();
    const sites = colorRuns(d.model, d.store, d.layer(1), undefined, { colors: d.colors });
    const more = colorRuns(d.model, d.store, d.layer(2), undefined, { colors: d.colors });
    // one color each, the first and the second of the categorical palette: neither `token:range`
    expect([sites.runs.length, more.runs.length]).toEqual([1, 1]);
    expect(more.runs[0].color).not.toBe(sites.runs[0].color);
    expect(more.runs[0].color).not.toBe(LightTokens.colors['seq-range']);
    // the scores' held domain is [0, 3): the union with `higher`; their extent alone would color 0.9 at the top
    const scores = colorRuns(d.model, d.store, d.layer(4), undefined, { colors: d.colors });
    const q = Observed[1] - 1;
    const top = scores.runs.find((x) => {
      const g = partOf(x.target, unit(0))!.ranges;
      for (let i = 0; i + 1 < g.length; i += 2) if (g[i] <= q && q < g[i + 1]) return true;
      return false;
    });
    const high = colorRuns(d.model, d.store, d.layer(5), undefined, { colors: d.colors });
    expect(top?.color).not.toBe(high.runs[high.runs.length - 1].color);
    expect(d.channel(4).at(datumAt(EntStart + q, 'values', ''))).toBe(top?.color);
  });

  it('light and dark: the dark tokens give a scheme’s dark variant; the reader’s substitution is the renderer’s', async () => {
    const d = await drawn();
    const r = colorRuns(d.model, d.store, d.layer(0), undefined, { colors: d.colors });
    expect(r.runs.some((x) => x.color !== x.dark)).toBe(true);
    const lit = r.runs.reduce((n, x) => n + positionCount(x.target), 0);
    expect(lit).toBe(4 * (Observed[1] - Observed[0]));
    const cvd = await d.frame({ size: Size, settled: true, cvdSafe: true });
    const field = colorRuns(d.model, d.store, d.layer(3), undefined, { colors: d.renderer.layerColors() });
    const k = Domains.category.indexOf('binding');
    const at = (cvd.layout.sections[0].tracks[3].layers[0].lay as { color: Channel }).color;
    expect(field.runs.map((x) => x.color)).toContain(at.at(datumAt(Domains.start[k], 'ranges', '')));
    expect(field.runs.map((x) => x.dark)).toContain(
      d.renderer
        .layerColors()
        .channel(d.layer(3), DarkTokens)!
        .at(datumAt(Domains.start[k], 'ranges', '')),
    );
  });
});

describe('data not held (EC-169)', () => {
  async function unsettled(settings?: { maxTilesPerWant: number; tileBins: number }) {
    const store = new DataStore(settings ? { settings } : {});
    const { model, view } = await openModel(p53Spec({ view: true }), store);
    const colors = layerColors(
      view!,
      viewTheme(view!, undefined, () => undefined),
    );
    return { model, store, view: view!, colors, layer: (t: number) => view!.sections[0].tracks[t].layers[0] };
  }

  it('EC-169: tiles not asked are pending, never an empty run', async () => {
    const u = await unsettled();
    const r = colorRuns(u.model, u.store, u.layer(4), undefined, { colors: u.colors });
    expect([r.complete, r.tooMany, Array.from(partOf(r.pending, Up)!.ranges)]).toEqual([false, false, [0, UpLength]]);
    expect(r.runs).toEqual([]);
  });

  it('EC-169: past the store’s tile bound, too many: pending and said, never complete with no runs', async () => {
    const u = await unsettled({ maxTilesPerWant: 1, tileBins: 16 });
    const r = colorRuns(u.model, u.store, u.layer(4), undefined, {
      colors: u.colors,
      want: (w) => u.store.want([w]),
    });
    await u.store.whenSettled();
    const again = colorRuns(u.model, u.store, u.layer(4), undefined, { colors: u.colors });
    for (const x of [r, again]) {
      expect([x.complete, x.tooMany, x.runs.length]).toEqual([false, true, 0]);
      expect(positionCount(x.pending)).toBe(UpLength);
    }
  });

  it('a layer with no colors to give is said', async () => {
    const u = await unsettled();
    const fake = { ...u.layer(0), type: 'ruler', fallbacks: [] };
    expect(colorRuns(u.model, u.store, fake, undefined, { colors: u.colors }).issues.map((i) => i.code)).toEqual([
      'seq.structure.layer',
    ]);
  });
});

describe('unitRows', () => {
  it('each item reaching a unit, with its fields; items past the structure left out', async () => {
    const d = await drawn();
    const r = unitRows(d.model, d.store, '/data/domains');
    expect([r.complete, r.tooMany]).toEqual([true, false]);
    expect(r.rows.map((x) => [x.key, x.fields.category])).toEqual([
      ['DBD', 'binding'],
      ['L1', 'linker'],
    ]);
    expect(r.rows[0].target.parts.map((p) => p.system)).toEqual([0, 1, 2, 3].map(unit));
    expect(unitRows(d.model, d.store, '/data/scores').issues.map((i) => i.code)).toEqual(['seq.structure.rows']);
  });
});
