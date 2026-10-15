/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Reveal (seqspec:interaction#focus): `ctx.reveal` and `seqv:reveal` bring a target into view without selecting,
// focusing or zooming: text scrolls its lines in, tracks pan at their cell (centered; a target wider than the view its
// start in), nothing where it is shown already. The `revealHighlight` behavior (seqspec:interaction#highlight) reveals a
// followed channel's highlight once it rests `dwellMs`; never in the section the pointer is in, during a press, while
// the viewer's own hover is held, nor for a channel it does not follow.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import type { SeqSet } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { Hit, KernelRenderers, type LayerRenderer, RevealHighlightBehavior, SeqCommands } from '../src/index.js';
import { driven } from './input-fixtures.js';
import { seqLoci } from './overlay-fixtures.js';

const S = 'v.main';
/** 2,000 residues, letters and a region over [1500, 1510), as `layout`. */
function doc(layout: 'tracks' | 'text'): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', value: 'ACDEFGHIKL'.repeat(200) });
  b.data.annotation({ ref: 'r', geometry: 'ranges', target: a, value: { start: [1500], stop: [1510] } });
  const v = b.view({ ref: 'v', label: 'V', layout, height: 'box' });
  v.block({ sequence: a });
  v.track({ ref: 'v.letters' }).letters();
  return b.getDocument();
}
async function open(layout: 'tracks' | 'text') {
  const t = await driven(doc(layout), {
    size: { [S]: [600, 300] },
    renderers: KernelRenderers as LayerRenderer[],
  });
  await t.ctx.graph.update();
  const set = t.ctx.graph.assembled.value?.set as SeqSet;
  const view = () => {
    const s = t.vp(S).live.current;
    return { cell: s.cell, x: s.scrollX, y: s.scrollY };
  };
  return { ...t, set, view };
}
/** Nothing selected and no focus changed by a reveal. */
const untouched = (t: Awaited<ReturnType<typeof open>>, focus: unknown) => {
  expect(t.ctx.selection.current.value.seq).toBeNull();
  expect(t.ctx.focus.current.value).toBe(focus);
};

describe('ctx.reveal', () => {
  it('tracks: pans at the cell to center the target; nothing where it is shown; selects and focuses nothing', async () => {
    const t = await open('tracks');
    t.vp(S).setView(12, 0, 'api');
    const focus = t.ctx.focus.current.value;
    expect(t.ctx.reveal(seqLoci(t.set, 'a', 1500, 1510), { animate: false })).toBe(true);
    const v = t.view();
    expect(v.cell).toBe(12);
    const ax = t.ctx.layout.axis(S);
    const mid = ((ax?.x(0, 1500 as never) ?? 0) + (ax?.x(0, 1510 as never) ?? 0)) / 2;
    const l = t.layout(S);
    expect(Math.abs(mid - (l.dataLeft + l.dataRight) / 2)).toBeLessThan(12 + 1);
    // Shown: no move.
    expect(t.ctx.reveal(seqLoci(t.set, 'a', 1502, 1506), { animate: false })).toBe(true);
    expect(t.view()).toEqual(v);
    untouched(t, focus);
  });

  it('tracks: a target wider than the view brings its start in, `pad` px clear', async () => {
    const t = await open('tracks');
    t.vp(S).setView(12, 0, 'api');
    t.ctx.reveal(seqLoci(t.set, 'a', 1000, 1900), { animate: false, pad: 20 });
    const ax = t.ctx.layout.axis(S);
    expect((ax?.x(0, 1000 as never) ?? 0) - t.layout(S).dataLeft).toBeCloseTo(20, 0);
    expect(t.view().cell).toBe(12);
  });

  it('text: scrolls the lines in, and not when they are shown; animated over the viewport animation', async () => {
    const t = await open('text');
    const focus = t.ctx.focus.current.value;
    expect(t.view().y).toBe(0);
    t.ctx.reveal(seqLoci(t.set, 'a', 1500, 1510));
    t.advance(1000);
    const y = t.view().y;
    expect(y).toBeGreaterThan(0);
    t.ctx.reveal(seqLoci(t.set, 'a', 1500, 1510));
    t.advance(1000);
    expect(t.view().y).toBe(y);
    untouched(t, focus);
  });

  it('seqv:reveal: the command, the same move', async () => {
    const t = await open('tracks');
    t.vp(S).setView(12, 0, 'api');
    await t.ctx.commands.dispatch(SeqCommands.Reveal, { target: { refs: ['r'] }, animate: false });
    expect(t.view().x).toBeGreaterThan(0);
    expect(t.view().cell).toBe(12);
    expect(t.ctx.selection.current.value.seq).toBeNull();
  });
});

describe('the revealHighlight behavior', () => {
  async function lit() {
    const t = await open('tracks');
    t.vp(S).setView(12, 0, 'api');
    t.ctx.behaviors.add(RevealHighlightBehavior, { animate: false });
    return t;
  }

  it('reveals a followed channel once it rests `dwellMs`; a newer one restarts the wait', async () => {
    const t = await lit();
    t.ctx.highlight.set('links:3d', seqLoci(t.set, 'a', 800, 805));
    t.advance(150);
    // A newer highlight before the dwell: the wait starts over, for it.
    t.ctx.highlight.set('links:3d', seqLoci(t.set, 'a', 1500, 1505));
    t.advance(150);
    expect(t.view().x).toBe(0);
    t.advance(60);
    const ax = t.ctx.layout.axis(S);
    const l = t.layout(S);
    const x = ax?.x(0, 1502 as never) ?? 0;
    expect(x > l.dataLeft && x < l.dataRight).toBe(true);
  });

  it('follows an external hover of an owner in `hovers` (a linked 3D view), not of another owner', async () => {
    const t = await lit();
    t.ctx.hover.setLoci(seqLoci(t.set, 'a', 1500, 1505), 'other');
    t.advance(300);
    expect(t.view().x).toBe(0);
    t.ctx.behaviors.get<{ hovers: readonly string[] }>(RevealHighlightBehavior.name)?.setParams({ hovers: ['3d'] });
    t.ctx.hover.setLoci(seqLoci(t.set, 'a', 1500, 1505), '3d');
    t.advance(150);
    // Moved on before the dwell: the wait is for the newer target.
    t.ctx.hover.setLoci(seqLoci(t.set, 'a', 1600, 1605), '3d');
    t.advance(150);
    expect(t.view().x).toBe(0);
    t.advance(60);
    const ax = t.ctx.layout.axis(S);
    const l = t.layout(S);
    const x = ax?.x(0, 1602 as never) ?? 0;
    expect(x > l.dataLeft && x < l.dataRight).toBe(true);
    // Cleared during its dwell: nothing.
    t.vp(S).setView(12, 0, 'api');
    t.ctx.hover.setLoci(seqLoci(t.set, 'a', 1700, 1705), '3d');
    t.advance(100);
    t.ctx.hover.setLoci(null, '3d');
    t.advance(300);
    expect(t.view().x).toBe(0);
  });

  it('not a channel it does not follow, not under the pointer, not during a press, not with its own hover', async () => {
    const t = await lit();
    t.ctx.highlight.set('host', seqLoci(t.set, 'a', 1500, 1505));
    t.advance(300);
    expect(t.view().x).toBe(0);
    // The pointer inside the section.
    t.pointer('move', S, 300, 40);
    t.ctx.highlight.set('links:3d', seqLoci(t.set, 'a', 1500, 1505));
    t.advance(300);
    expect(t.view().x).toBe(0);
    // Out of the section, and this viewer's own hover held (a link lighting it back): nothing.
    t.ctx.input.leave(t.ctx.layout.sectionIndex(S));
    expect(t.ctx.pointer.current.inside).toBe(false);
    const hit = Hit.create();
    t.ctx.layout.hitTest(S, 300, 40, hit);
    t.ctx.hover.setFromHit(hit, 'keyboard');
    t.ctx.highlight.set('links:3d', seqLoci(t.set, 'a', 1600, 1605));
    t.advance(300);
    expect(t.view().x).toBe(0);
    // A press: nothing either.
    t.ctx.hover.clear('keyboard');
    t.pointer('down', S, 300, 40);
    t.ctx.input.leave(t.ctx.layout.sectionIndex(S));
    t.ctx.highlight.set('links:3d', seqLoci(t.set, 'a', 1700, 1705));
    t.advance(300);
    expect(t.view().x).toBe(0);
    t.pointer('up', S, 300, 40);
    t.ctx.input.leave(t.ctx.layout.sectionIndex(S));
    // Then, all clear: revealed.
    t.ctx.highlight.set('links:3d', seqLoci(t.set, 'a', 1600, 1605));
    t.advance(300);
    expect(t.view().x).toBeGreaterThan(0);
  });
});
