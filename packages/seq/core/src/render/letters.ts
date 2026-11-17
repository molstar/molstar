/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf } from '@molstar/core/repr/tree';
import type { LettersTable } from '@molstar/seqx-core/data';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { type DensityStep, densityStep } from './config.js';
import { NoColor, type PassName, type PassWriter } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import { type LayerSource, layerSource, placeholderSpans, readTiles } from './layer-data.js';
import type { Box, Emphasis, Hit, LayerType, LayoutContext, Overlay, StoreReader } from './layer-type.js';
import { drawSpans, outlineCell } from './marks.js';
import { inkOn, letterToken, resolveColor, type SeqTokens, withAlpha } from './tokens.js';

// The `letters` layer type (_spec/README.md "Render", Layer types): a section's letters, a cell a position. By density
// step: letters on their colored cells, cells only, a tick a position, then one bar a piece in view (presence needs no
// data). Letters are centered in their cells, one size a row, ink by the fill.

export interface LettersPrep {
  readonly src?: LayerSource;
  readonly alphabet: 'protein' | 'dna' | 'rna';
  /** A color the spec gives every cell, as written (resolved per token table at layout). */
  readonly color?: unknown;
}

export interface LettersLay {
  key: string;
  system: string;
  box: Box;
  axis: AxisTransform;
  step: DensityStep;
  /** Cells in view: position, letter code (0: not loaded), fill, ink, snapped x extents. */
  n: number;
  pos: Int32Array;
  code: Uint16Array;
  fill: Uint32Array;
  ink: Uint32Array;
  x0: Float32Array;
  x1: Float32Array;
  font: string;
  fontSize: number;
  /** `x0, x1, cut` triples: the bars of the coarse steps, the pending and failed placeholders. */
  bars: number[];
  pending: number[];
  failed: number[];
  tokens: SeqTokens;
}

function grow<T extends Int32Array | Uint16Array | Uint32Array | Float32Array>(a: T, n: number): T {
  return n <= a.length ? a : new (a.constructor as new (n: number) => T)(Math.max(n, 2 * a.length));
}

/** The letter code at a position of a table, or 0. */
function codeAt(t: LettersTable, p: number): number {
  const i = p - t.start;
  if (i < 0 || i >= t.letters.length) return 0;
  if (typeof t.letters === 'string') return t.letters.charCodeAt(i);
  return t.alphabet ? t.alphabet.charCodeAt(t.letters[i]) || 0 : 0;
}

export const LettersType: LayerType<LettersPrep, LettersLay> = {
  name: 'letters',

  prepare(layer: ResolvedLayer, model: StaticModel): LettersPrep {
    const src = layerSource(model, layer);
    const system = src ? model.systems.get(src.dataSystem) : undefined;
    return {
      ...(src ? { src } : {}),
      alphabet: system?.alphabet ?? 'protein',
      ...(layer.theme.color !== undefined ? { color: layer.theme.color } : {}),
    };
  },

  layout(prep, ctx: LayoutContext, data: StoreReader, _overlay: Overlay, prev?: LettersLay): LettersLay {
    const { axis, box, config, tokens } = ctx;
    const lay: LettersLay = prev ?? {
      key: ctx.key,
      system: axis.system,
      box,
      axis,
      step: 'cells',
      n: 0,
      pos: new Int32Array(64),
      code: new Uint16Array(64),
      fill: new Uint32Array(64),
      ink: new Uint32Array(64),
      x0: new Float32Array(64),
      x1: new Float32Array(64),
      font: '',
      fontSize: 0,
      bars: [],
      pending: [],
      failed: [],
      tokens,
    };
    lay.key = ctx.key;
    lay.system = axis.system;
    lay.box = box;
    lay.axis = axis;
    lay.tokens = tokens;
    lay.n = 0;
    lay.bars.length = 0;
    lay.pending.length = 0;
    lay.failed.length = 0;
    const src = prep.src;
    // letters are the section's own: through a path they would be another system's cells
    if (!src || src.path) {
      if (src?.path)
        ctx.report(
          issueOf(RenderIssues, 'seq.render.layer.letters', 'Letters are drawn on their own system only', {
            address: ctx.key,
            hint: RenderIssues.codes['seq.render.layer.letters'].hint,
          }),
        );
      return lay;
    }
    const step = densityStep(axis.devicePx, config.density);
    lay.step = step;
    const vis = axis.visible([]);
    if (step === 'bars' || step === 'density') {
      for (let i = 0; i + 1 < vis.length; i += 2) axis.spans(vis[i], vis[i + 1], lay.bars);
      return lay;
    }

    let fixed: number | undefined;
    if (prep.color !== undefined) {
      fixed = resolveColor(prep.color, tokens);
      if (fixed === undefined) reportColor(ctx);
    }
    const read = readTiles(ctx, data, src, vis, 1);
    const tables = read.entries.map((e) => e.tables.letters).filter((t): t is LettersTable => !!t);
    tables.sort((a, b) => a.start - b.start);
    let count = 0;
    for (let i = 0; i + 1 < vis.length; i += 2) count += vis[i + 1] - vis[i];
    lay.pos = grow(lay.pos, count);
    lay.code = grow(lay.code, count);
    lay.fill = grow(lay.fill, count);
    lay.ink = grow(lay.ink, count);
    lay.x0 = grow(lay.x0, count);
    lay.x1 = grow(lay.x1, count);
    let t = 0;
    for (let i = 0; i + 1 < vis.length; i += 2) {
      for (let p = vis[i]; p < vis[i + 1]; p++) {
        while (t < tables.length && tables[t].start + tables[t].letters.length <= p) t++;
        let code = 0;
        for (let k = t; k < tables.length && tables[k].start <= p && !code; k++) code = codeAt(tables[k], p);
        if (!code) continue;
        const j = lay.n++;
        const ch = String.fromCharCode(code);
        const fill = fixed ?? tokens.colors[letterToken(ch, prep.alphabet)];
        lay.pos[j] = p;
        lay.code[j] = code;
        lay.fill[j] = fill;
        lay.ink[j] = inkOn(fill, tokens);
        lay.x0[j] = axis.snap(axis.x(p));
        lay.x1[j] = axis.snap(axis.x(p + 1));
      }
    }
    placeholderSpans(src, axis, read.pending, lay.pending);
    placeholderSpans(src, axis, read.failed, lay.failed);

    // one size a row: the largest of the scale whose advance fits the cell and whose height fits the row
    lay.fontSize = 0;
    if (step === 'letters') {
      const cell = axis.scale;
      for (const s of config.letterSizes)
        if (s * config.monoAdvance <= cell - 2 && s <= box.height * 0.75) lay.fontSize = s;
      lay.font = `${lay.fontSize}px mono`;
    }
    return lay;
  },

  draw(lay, pass: PassName, out: PassWriter, emphasis?: Emphasis) {
    const { box, tokens } = lay;
    if (pass === 'data') {
      drawSpans(out, lay.bars, box.y, box.height, tokens.colors['seq-sequence']);
      drawSpans(out, lay.pending, box.y, box.height, tokens.colors['seq-pending']);
      drawSpans(out, lay.failed, box.y, box.height, withAlpha(tokens.colors['seq-conflict'], 0.25));
      const dots = lay.step === 'dots';
      const tick = Math.max(2, Math.round(box.height * 0.3));
      for (let j = 0; j < lay.n; j++) {
        const w = Math.max(lay.x1[j] - lay.x0[j], 1 / lay.axis.dpr);
        if (dots) out.rect(lay.x0[j], box.y + (box.height - tick) / 2, w, tick, lay.fill[j]);
        else out.rect(lay.x0[j], box.y, w, box.height, lay.fill[j]);
      }
      if (lay.step === 'letters' && lay.fontSize > 0 && lay.n > 0) {
        const y = box.y + box.height / 2 + lay.fontSize * 0.35;
        out.beginGlyphs(lay.font, y, NoColor);
        for (let j = 0; j < lay.n; j++) out.letter(lay.code[j], (lay.x0[j] + lay.x1[j]) / 2, lay.ink[j]);
        out.endGlyphs();
      }
      return;
    }
    if (pass === 'hover') {
      const h = emphasis?.hover;
      if (!h || h.layer !== lay.key || h.position === undefined) return;
      outlineCell(out, lay.axis, h.position, box, tokens.colors['seq-hover'], 1.5);
    }
  },

  hit(lay, x, y): Hit | undefined {
    const { box } = lay;
    if (y < box.y || y >= box.y + box.height) return undefined;
    const position = lay.axis.positionAt(x);
    if (position < 0) return undefined;
    return { layer: lay.key, position, system: lay.system };
  },
};

export function reportColor(ctx: LayoutContext): void {
  ctx.report(
    issueOf(RenderIssues, 'seq.render.color', 'The layer’s color is neither a hex color nor a token', {
      address: ctx.key,
      hint: RenderIssues.codes['seq.render.color'].hint,
    }),
  );
}
