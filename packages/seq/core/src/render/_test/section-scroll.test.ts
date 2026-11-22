/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A tracks section scrolls past its cap (_spec/README.md "Render", Section scroll): its box its `height`, else its
// content up to `section.maxPx`; the tracks below its ruler scroll by the viewport's `top` under the fixed ruler, with
// the one rail at its right edge; a hidden band is no hit; an export shows every track.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { hitView, SectionRail, TracksNav } from '../index.js';
import { renderSvg } from '../render-svg.js';
import { settle } from './helpers.js';

const Size = { width: 800, height: 2000 };

/** A ruler and `n` gene tracks of one feature each; `height` the section's own. */
export function geneTracks(n = 30, height?: number): Node {
  const b = createSeqSpec({ title: 'tracks past the cap' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(100) });
  b.views.view('main', (v) => {
    const sec = v.section('s', { of: s, ...(height !== undefined ? { height } : {}) });
    sec.track('ruler', { label: '' }).layer('r', { type: 'ruler' });
    for (let i = 0; i < n; i++) {
      const d = b.data.ranges(`g${i}`, { target: s, rows: { start: [1 + i], end: [200 + i], id: [`gene ${i}`] } });
      sec.track(`t${i}`, { label: `Gene ${i}` }).layer('l', { type: 'ranges', data: d });
    }
  });
  return b.build();
}

async function laid(spec: Node, config = {}, top?: number) {
  const s = await settle(spec, [], config);
  const f0 = await s.frame({ size: Size, settled: true });
  const start = f0.layout.sections[0].axis!.u0;
  const end = f0.layout.sections[0].axis!.u1;
  const f = top === undefined ? f0 : await s.frame({ size: Size, settled: true, viewport: { start, end, top } });
  return { f, sec: f.layout.sections[0] };
}

describe('a tracks section past its cap', () => {
  it('takes `section.maxPx` (600) by default, its content scrolling inside with the one rail', async () => {
    const { sec } = await laid(geneTracks());
    const p = sec.placement;
    expect(p.contentHeight).toBeGreaterThan(600);
    expect(p.box.height).toBe(600);
    expect(p.clip).toBe(true);
    expect(sec.height).toBe(600);
    expect(sec.sectionRail?.hit.navigator).toBe(SectionRail);
    expect(sec.sectionRail!.end - sec.sectionRail!.start).toBeCloseTo(p.box.y + 600 - p.scrollY!, 6);
  });

  it('a host’s `section.maxPx` (300) and a spec `height` (200) are its box', async () => {
    expect((await laid(geneTracks(), { section: { maxPx: 300 } })).sec.placement.box.height).toBe(300);
    expect((await laid(geneTracks(30, 200))).sec.placement.box.height).toBe(200);
    // under the cap: its content, no rail
    const few = (await laid(geneTracks(3))).sec;
    expect(few.placement.clip).toBe(false);
    expect(few.sectionRail).toBeUndefined();
  });

  it('scrolls the tracks below the ruler by `top`, the ruler fixed; a band hidden under it is no hit', async () => {
    const a = (await laid(geneTracks())).sec;
    const { f, sec: b } = await laid(geneTracks(), {}, 200);
    expect(b.placement.top).toBe(200);
    expect(b.tracks[0].box.y).toBe(a.tracks[0].box.y);
    expect(b.tracks[5].box.y).toBeCloseTo(a.tracks[5].box.y - 200, 6);
    expect(b.sectionRail!.start).toBe(200);
    // clamped at the content's end
    const c = (await laid(geneTracks(), {}, 1e6)).sec;
    expect(c.placement.top).toBeCloseTo(c.placement.contentHeight - 600, 6);
    // a point on the ruler is the ruler's, though a scrolled band lies under it
    const r = b.tracks[0].box;
    const h = hitView(f.layout, r.x + 40, r.y + r.height / 2);
    expect(h?.layer ?? '').not.toMatch(/\/t\d+\//);
  });

  it('the navigation clamps `top` to its content', async () => {
    const { sec } = await laid(geneTracks());
    const c = { total: 400, width: 600, height: 600, cyclic: false, maxScale: 40, placement: sec.placement };
    const v = { start: 0, end: 400 };
    expect(TracksNav.pan(c, v, 0, 1e6).top).toBeCloseTo(sec.placement.contentHeight - 600, 6);
    expect(TracksNav.pan(c, v, 0, -50).top).toBe(0);
  });

  it('an export shows every track (the figure), no rail', async () => {
    const r = await renderSvg(geneTracks(), { size: { width: 800, height: 'auto' } } as never);
    expect(r.svg).toContain('>Gene 29<');
    expect(r.svg).toContain('>Gene 0<');
  });
});
