/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An ESPript-like figure as a step, for the heap guard and its quick twin (text-rows.heap.ts, text-rows.quick.ts;
// seqspec:viewer#rows-in-text, seqspec:viewer#frames, seqspec:viewer#performance), over text-rows-fixtures.ts's globins:
// `draw`, the text scrolled down and back a step each, its lines redrawn (the rows' cells and letters, their labels and
// numbers, a box on a row, the frames' outlines and titles, the helices' cylinders and labels, glyph markers); `hover`,
// the pointer over the rows' cells, a gap and a frame's title through the real dispatcher, an overlay frame each;
// `tracks`, the same document in tracks (`msa`) panned left and right a step each (a copy, the strip and its frames);
// `zoom`, its rows zoomed in and out a step each (a rescale, no build: the frames follow). Each call loads the viewer
// afresh (`vi.resetModules`).

import { vi } from 'vitest';

export type TextRowsKind = 'draw' | 'hover' | 'tracks' | 'zoom';
const S = 'figure.main';

export async function textRows(kind: TextRowsKind) {
  vi.resetModules();
  const { driven } = await import('./input-fixtures.js');
  const { COPY_BACK, recording } = await import('./composer-fixtures.js');
  const { figureDoc } = await import('./text-rows-fixtures.js');
  const { Hit, HoverKind, Inputs, KernelRenderers, PointerType } = await import('../src/index.js');
  const { Dirty, SectionComposer } = await import('../src/draw/index.js');
  const linear = kind === 'tracks' || kind === 'zoom';
  const doc = figureDoc(linear ? { style: { preset: 'msa' } } : {});
  const t = await driven(doc, {
    size: { [S]: linear ? [900, 400] : [900, 260] },
    renderers: KernelRenderers as never,
    labels: false,
  });
  await t.ctx.graph.update();
  // The frames alone: the hover label's strings are one allocation per hovered element by design (IN2's guard).
  t.ctx.behaviors.get('hoverLabel')?.remove();
  const l = t.layout(S);
  const tr = l.textRows;
  if (!linear && (tr.layer < 0 || l.lines.count < 2 || l.scrolledHeight <= l.height))
    throw new Error('the rows are not text, or the text does not scroll');
  if (linear && l.rows.layer < 0) throw new Error('no rows');
  const composer = new SectionComposer(t.ctx, S);
  const be = kind === 'hover' || kind === 'tracks' ? recording(COPY_BACK) : recording();
  const vp = t.vp(S);
  if (kind === 'tracks') vp.setView(12, 0, 'api');
  composer.frame(be, Dirty.all);
  // Hover: points over the first line's rows (cells and gaps) and the frames' title band.
  const pts = new Float64Array(64);
  const ax = t.ctx.layout.axis(S)?.axis;
  if (kind === 'hover') {
    const cell = ax?.fixedCell ?? 0;
    let band = -1;
    let rows = -1;
    for (let e = 0; e < l.entryCount; e++) {
      if (l.entryRow[e] !== 0 || !(l.entryHeight[e] > 0)) continue;
      if (l.entryTrack[e] === tr.track) rows = e;
      else if (l.entryTrack[e] === tr.track - 1) band = e;
    }
    if (rows < 0 || band < 0) throw new Error('no rows or title band on the first line');
    for (let i = 0; i < 32; i++) {
      // Columns within the first group (no group gap to step over), rows by turn; every eighth on the title band.
      pts[2 * i] = l.dataLeft + l.lines.left[0] + ((i * 3) % 10) * cell + cell / 2 + 0.25;
      pts[2 * i + 1] =
        i % 8 === 7
          ? l.scrolledTop + l.entryTop[band] + l.entryHeight[band] / 2
          : l.scrolledTop + l.entryTop[rows] + ((i % tr.count) + 0.5) * tr.height;
    }
    const h = Hit.create();
    let residues = 0;
    let items = 0;
    for (let i = 0; i < 32; i++) {
      if (!t.ctx.layout.hitTest(S, pts[2 * i], pts[2 * i + 1], h)) continue;
      if (h.kind === HoverKind.aligned) residues++;
      else if (h.kind === HoverKind.item) items++;
    }
    if (residues < 20 || items < 2) throw new Error(`the points hit ${residues} residues and ${items} items`);
  }
  const p = Inputs.pointer();
  p.section = t.ctx.layout.sectionIndex(S);
  p.type = PointerType.move;
  p.button = -1;
  p.buttons = 0;
  let hovered = 0;
  let paths = 0;
  let texts = 0;
  const step = (i: number) => {
    be.reset();
    if (kind === 'hover') {
      p.x = pts[2 * (i & 31)];
      p.y = pts[2 * (i & 31) + 1];
      t.ctx.input.pointer(p);
      if (t.ctx.hover.live.current.index >= 0) hovered++;
      composer.frame(be, Dirty.overlay);
      return;
    }
    if (kind === 'draw') vp.panBy(0, i & 1 ? -11 : 11, 'pointer');
    else if (kind === 'tracks') vp.panBy(i & 1 ? -7 : 7, 0, 'pointer');
    else vp.setRowScale(i & 1 ? 1 : 1.25, 0);
    composer.frame(be, Dirty.tracks);
    paths += be.pathItems.n;
    texts += be.textItems.n;
  };
  return {
    step,
    /** What the steps drew: paths (frames, helices) and texts over the steps run, the hovers. */
    state: () => ({ hovered, paths, texts }),
    /** Counts anew (a guard counts over its measured steps alone). */
    resetCounts: () => {
      hovered = paths = texts = 0;
    },
  };
}
