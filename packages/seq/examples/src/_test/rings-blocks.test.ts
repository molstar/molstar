/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Complex's blocks drawn against their data (the Plasmid's are in `plasmid.test.ts`): copies merged ×2, a feature
// on its own chain only, split, stacked and on a ring.

import type { Node } from '@molstar/core/repr/tree';
import { DataStore } from '@molstar/seqx-core/data';
import {
  attachSources,
  type BlockPlacement,
  bundledTextMetrics,
  hitView,
  LightTokens,
  layoutView,
  type RingPlacement,
} from '@molstar/seqx-core/render';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import { ComplexLengths } from '../complex/blocks.js';
import { complexContacts, complexPae } from '../complex/pae.js';
import { ComplexViews } from '../views.js';

const view = (id: string) => ComplexViews.find((v) => v.id === id)!;

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
