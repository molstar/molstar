/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Plasmid ring drawn against its data: the stand-in's anchors (every feature keeps its length under every rotation;
// bla and the primer site cross the rotated origin), the ring's record (its window over the origin, its detail in two
// pieces), each feature's arc where its positions are (across the origin too), the rotated linear view's two pieces.
// The Complex's blocks are the Complex's tests'.

import { existsSync, readFileSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { createSeqContext } from '@molstar/seqx-core/context';
import type { Node } from '@molstar/core/repr/tree';
import { DataStore, SeqRuntime } from '@molstar/seqx-core/data';
import {
  attachSources,
  bundledTextMetrics,
  hitView,
  LightTokens,
  layoutView,
  recordDrawLists,
  type RingPlacement,
} from '@molstar/seqx-core/render';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import { PlasmidFeatures, PlasmidLength, plasmidLetters, RotatedOrigin, rotate } from '../plasmid/data.js';
import { ComplexViews, PlasmidViews } from '../views.js';

const update = process.env.SEQ_RECORDS === 'update';

function record(name: string, text: string) {
  const file = fileURLToPath(new URL(`./records/${name}.txt`, import.meta.url));
  if (update) writeFileSync(file, text);
  expect(existsSync(file), `no ${file}: run with SEQ_RECORDS=update`).toBe(true);
  expect(text).toBe(readFileSync(file, 'utf8'));
}

const context = () => createSeqContext({ registry: [{ kind: 'format-runtime', runtime: SeqRuntime }] });
const view = (id: string) => [...PlasmidViews, ...ComplexViews].find((v) => v.id === id)!;

/** A view laid out once its (inline) data has arrived. */
async function laid(spec: Node, viewports = {}, size = { width: 720, height: 900 }) {
  const model = resolveSpec(spec);
  const v = resolveViews(model)[0];
  const store = new DataStore();
  await attachSources(model, store);
  const input = { model, view: v, size, viewports, tokens: LightTokens, text: bundledTextMetrics(), data: store };
  for (let i = 0; i < 3; i++) {
    store.want(layoutView(input).wants);
    await store.whenSettled();
  }
  return layoutView(input);
}

describe('the Plasmid stand-in', () => {
  it('every feature keeps its length under every rotation; bla crosses the rotated origin', () => {
    expect(plasmidLetters().length).toBe(PlasmidLength);
    for (let origin = 0; origin < PlasmidLength; origin += 1) {
      for (const f of PlasmidFeatures) {
        const s = rotate(f.start, origin);
        const e = rotate(f.end, origin);
        expect((e - s + PlasmidLength) % PlasmidLength).toBe(f.end - f.start);
      }
    }
    const bla = PlasmidFeatures.find((f) => f.id === 'bla')!;
    expect(rotate(bla.start, RotatedOrigin)).toBeGreaterThan(rotate(bla.end, RotatedOrigin));
  });

  it('the ring draws each feature where its positions are, across the origin one feature', async () => {
    const layout = await laid(view('ring').build());
    const ring = layout.sections[0];
    const p = ring.placement as RingPlacement;
    const at = (track: number, pos: number) => {
      const keys = new Set<string>();
      const k = track;
      for (let y = p.strips.bandY[k] + 1; y < p.strips.bandY[k] + p.strips.bandH[k]; y += 2) {
        p.warp.map(p.frameX(pos + 0.5), y);
        const h = hitView(layout, p.warp.px, p.warp.py);
        if (h?.element) keys.add(String(h.element.key));
      }
      return [...keys].sort();
    };
    const bla = PlasmidFeatures.find((f) => f.id === 'bla')!;
    // bla on the minus ring both sides of the rotated origin
    expect(at(1, rotate(bla.start + 10, RotatedOrigin))).toContain('bla');
    expect(at(1, rotate(bla.end - 10, RotatedOrigin))).toContain('bla');
    expect(at(1, 0)).toContain('bla');
    // the primer site on the plus ring across the origin
    expect(at(0, PlasmidLength - 5)).toContain('primer');
    expect(at(0, 5)).toContain('primer');
    // one label a feature, placed outside
    const labels = Array.from({ length: p.labels.m }, (_, k) => p.labels.text[p.labels.placed[k]]);
    expect(labels.filter((t) => t === 'bla (AmpR)').length).toBe(1);
    expect(p.title).toBe('pUC19 (stand-in)');
    expect(p.subtitle).toBe('2,686 bp');
  });

  it('the detail follows the ring’s window over the origin, in two pieces; the ring’s record', async () => {
    const viewports = { '/views/main/ring': { start: 2600, end: 2780 } };
    const layout = await laid(view('ring').build(), viewports);
    expect(layout.sections[1].axis!.pieces).toEqual([2600, 2686, 0, 94]);
    const ctx = context();
    const r = await ctx.renderSvg(view('ring').build(), { size: { width: 520, height: 720 }, viewports, window: true });
    expect(r.issues).toEqual([]);
    expect(r.settled).toBe(true);
    record('plasmid-ring', recordDrawLists(r.lists));
    ctx.retire();
  });

  it('the rotated linear view: two pieces meeting at the seam, no break, the features whole', async () => {
    const layout = await laid(view('rotated-linear').build());
    const s = layout.sections[0];
    expect(s.axis!.pieces).toEqual([2000, 2686, 0, 700]);
    expect(s.breaks).toEqual([]);
    expect(s.axis!.total).toBe(686 + 700);
  });
});
