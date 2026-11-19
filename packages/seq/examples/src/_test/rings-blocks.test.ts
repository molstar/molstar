/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Plasmid ring and the Complex's blocks drawn against their data: the stand-in's anchors (every feature keeps its
// length under every rotation; bla and the primer site cross the rotated origin), the ring's record (its window over
// the origin, its detail in two pieces), each feature's arc where its positions are (across the origin too), the
// rotated linear view's two pieces, and the Complex's blocks (copies merged ×2, a feature on its own chain only).

import { existsSync, readFileSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { createSeqContext } from '@molstar/seqx-core/context';
import type { Node } from '@molstar/core/repr/tree';
import { DataStore, SeqRuntime } from '@molstar/seqx-core/data';
import {
  attachSources,
  type BlockPlacement,
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
import { ComplexLengths } from '../complex/blocks.js';
import { complexContacts, complexPae } from '../complex/pae.js';
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

describe('the Complex stand-in’s blocks', () => {
  it('copies merged ×2, the ligands after the polymers; each feature on its own chain', async () => {
    const layout = await laid(view('blocks').build());
    const sec = layout.sections[0];
    expect(sec.section.blocks!.map((b) => b.label)).toEqual(['CDK2 (A) ×2', 'cyclin A (B) ×2', 'ATP ×2']);
    expect(sec.section.axis.ranges).toEqual([0, ComplexLengths.cdk2 + ComplexLengths.cyclin + 1]);
    const p = sec.placement as BlockPlacement;
    const keysAt = (strip: number, pos: number) => {
      const a = p.strips.axes[strip];
      const keys = new Set<string>();
      for (let y = p.strips.bandY[strip * 2]; y < p.strips.bandY[strip * 2] + p.strips.bandH[strip * 2]; y += 2) {
        const h = hitView(layout, a.xOfUnit(pos - p.strips.first[strip] + 0.5), y);
        if (h?.element) keys.add(String(h.element.key));
      }
      return [...keys];
    };
    expect(keysAt(0, 50)).toEqual(['psta']);
    expect(keysAt(1, ComplexLengths.cdk2 + 50)).toEqual(['box1']);
    expect(keysAt(1, ComplexLengths.cdk2 + 135)).toEqual([]);
  });

  it('split: a block a chain and a ligand; stacked: a row each; on a ring: an arc each', async () => {
    expect((await laid(view('blocks-split').build())).sections[0].section.blocks!.length).toBe(6);
    const stackedLayout = await laid(view('blocks-stacked').build());
    const stacked = stackedLayout.sections[0].placement as BlockPlacement;
    expect([stacked.stacked, stacked.strips.n]).toEqual([true, 3]);
    // each row's Sequence its own block's letters, in its own row
    const st = stacked.strips;
    for (let i = 0; i < st.n; i++) {
      const y = st.bandY[i * st.bands + 1] + 2;
      expect(hitView(stackedLayout, st.axes[i].xOfUnit(0.5), y)).toMatchObject({ block: i, system: st.axes[i].system });
    }
    const ring = (await laid(view('blocks-ring').build())).sections[0].placement as RingPlacement;
    expect(ring.strips.n).toBe(3);
  });
});

describe('the Complex’s synthetic PAE and contacts', () => {
  it('reads low across the interface and high elsewhere; every contact lies in the interface', () => {
    const v = complexPae();
    const ny = ComplexLengths.cyclin;
    const mean = (i0: number, i1: number, j0: number, j1: number) => {
      let s = 0;
      for (let i = i0; i < i1; i++) for (let j = j0; j < j1; j++) s += v[i * ny + j];
      return s / ((i1 - i0) * (j1 - j0));
    };
    expect(mean(44, 56, 60, 95)).toBeLessThan(8);
    expect(mean(150, 165, 95, 125)).toBeLessThan(8);
    expect(mean(200, 290, 150, 250)).toBeGreaterThan(22);
    const c = complexContacts();
    for (let k = 0; k < c.a_start.length; k++) expect(v[c.a_start[k] * ny + c.b_start[k]]).toBeLessThan(8);
  });
});
