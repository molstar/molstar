/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf } from '@molstar/core/repr/tree';
import { DefaultMaxBlockValues } from '@molstar/seqx-core/data';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import { type Channel, glyphAgainst } from '@molstar/seqx-core/theme/index';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { channelForms, type LayerChannels, PlainDatum } from './channels.js';
import { DefaultRenderConfig, type RenderConfig } from './config.js';
import { drawPass, grown, NoColor, type PassName, type PassWriter, passDraws } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import { type LayerSource, layerSource } from './layer-data.js';
import {
  type Box,
  type Emphasis,
  type Hit,
  type LayerType,
  type LayoutContext,
  type Overlay,
  sectionId,
  type StoreReader,
} from './layer-type.js';
import { type Alphabet, LettersColorFallback } from './letters.js';
import { hoverCell } from './marks.js';
import { Pyramid, pyramidLevelFor } from './pyramid.js';
import { readHeld, rowsSource } from './rows.js';
import { type HeldCells, logoHeights, type Profile, profileOf } from './rows-model.js';
import type { SeqTokens } from './tokens.js';

// The `logo` layer type (_spec/README.md "Render", The logo): a conservation logo over an alignment's columns, each
// column's letters stacked by `logo_scale` (`bits`: log2(k) − H, less the small-sample correction with `correction`;
// `probability`), lowest first, each letter placed and stretched with the glyph runs' per-letter y and vertical scale.
// Its profile is computed from the held cells (`gaps`), kept by their version. Below 2 device px a column, one bar a
// device column, the tallest stack there (the shared reduction's one-axis max), in its top letter's color. A letter's
// color is its class's, saturated to `seq-logo-chroma` and held at 3:1 against the track (the `ink` rule with the surface as the fill).

const ColorFallback = LettersColorFallback;
const Env: Record<Alphabet, { alphabet: Alphabet }> = {
  protein: { alphabet: 'protein' },
  dna: { alphabet: 'dna' },
  rna: { alphabet: 'rna' },
};

/** A letter's cap height as a share of its font size: what a placed letter's vertical scale stretches. */
const CapHeight = 0.72;
/**
 * The least glyph height (px) a stacked letter is drawn as a letter: below it the bottom of a stack is an unreadable
 * pile, so the letter is a flat block of its color, its column wide and its height tall (the defaults audit).
 */
export const MinLogoGlyphPx = 9;

export interface LogoPrep {
  readonly src?: LayerSource;
  readonly alphabet: Alphabet;
  readonly channels: LayerChannels;
  readonly scale: 'bits' | 'probability';
  readonly correction: 'none' | 'small-sample';
  readonly gaps: 'ignore' | 'count';
}

export class LogoLay {
  key = '';
  section = '';
  system = '';
  box: Box = { x: 0, y: 0, width: 0, height: 0 };
  axis!: AxisTransform;
  tokens!: SeqTokens;
  hoverWidth = 1.5;
  hoverOffset = 1;
  /** What the hover's squircle reads (its rule). */
  config: RenderConfig = DefaultRenderConfig;
  profile: Profile | undefined;
  profileOf: HeldCells | undefined;
  channel: Channel | undefined;
  /** Placed letters: x, baseline, vertical scale, color, code; `n` of them. */
  n = 0;
  lx = new Float32Array(64);
  ly = new Float32Array(64);
  vs = new Float32Array(64);
  color = new Uint32Array(64);
  code = new Uint8Array(64);
  /** Placed letters under `MinLogoGlyphPx` tall: drawn as a flat block `cellW` wide (kept for the hit); `flatN` of them. */
  flat = new Uint8Array(64);
  flatN = 0;
  cellW = 0;
  font = '20px mono';
  fontSize = 20;
  /** Bars (below 2 device px a column): x, y, w, h, color; `m` of them. */
  m = 0;
  bars = new Float32Array(64);
  barColor = new Uint32Array(64);
  /** Columns in view and each column's stack (its letters' tops), for the hit's letter. */
  pending = false;
  readonly heights = new Float32Array(32);
  readonly datum = new PlainDatum();
  totals = new Float32Array(0);
  tops = new Uint8Array(0);
  bars1d: Pyramid | undefined;
  pxPerUnit = 0;
  /** A color's glyph color against the track (saturated, then the `ink` rule's 3:1 with the surface as the fill), by the token table's version. */
  glyphs = new Map<number, number>();
  glyphsVersion = -1;
}

export const LogoType: LayerType<LogoPrep, LogoLay> = {
  name: 'logo',

  prepare(layer: ResolvedLayer, model: StaticModel): LogoPrep {
    const src = layerSource(model, layer);
    const system = model.systems.get(layer.data?.system ?? '');
    const p = layer.params;
    return {
      ...(src ? { src } : {}),
      alphabet: system?.alphabet ?? 'protein',
      channels: channelForms(layer.theme),
      scale: p.logo_scale === 'probability' ? 'probability' : 'bits',
      correction: p.correction === 'small-sample' ? 'small-sample' : 'none',
      gaps: p.gaps === 'count' ? 'count' : 'ignore',
    };
  },

  layout(prep, ctx: LayoutContext, data: StoreReader, _overlay: Overlay, prev?: LogoLay): LogoLay {
    const { axis, box, tokens, config } = ctx;
    const lay = prev ?? new LogoLay();
    lay.key = ctx.key;
    lay.section = sectionId(ctx.section);
    lay.system = axis.system;
    lay.box = box;
    lay.axis = axis;
    lay.tokens = tokens;
    lay.hoverWidth = config.outline.hover;
    lay.hoverOffset = config.emphasis.hoverOffset;
    lay.config = config;
    lay.n = 0;
    lay.flatN = 0;
    lay.m = 0;
    lay.pending = false;
    if (!prep.src) return lay;
    const at = rowsSource(prep.src, data);
    if (at.rows === undefined || at.columns === undefined || at.rows * at.columns > DefaultMaxBlockValues) {
      // a lazy alignment's profile is a source's (served `matrix`), not read yet: the fallback's bars are none
      if (at.rows !== undefined)
        ctx.report(
          issueOf(RenderIssues, 'seq.render.logo.profile', 'The logo has no profile of a lazy alignment', {
            address: ctx.key,
            hint: RenderIssues.codes['seq.render.logo.profile'].hint,
          }),
        );
      lay.pending = at.rows === undefined;
      return lay;
    }
    const held = readHeld(ctx, data, at.src, at.rows, at.columns);
    if (!held) {
      lay.pending = true;
      return lay;
    }
    if (lay.profileOf !== held) {
      lay.profileOf = held;
      lay.profile = profileOf(held, prep.gaps);
      const C = held.columns;
      lay.totals = new Float32Array(C);
      lay.tops = new Uint8Array(C);
      const hs = new Float32Array(lay.profile.letters.length);
      for (let j = 0; j < C; j++) {
        logoHeights(lay.profile, j, prep.scale, prep.correction, alphabetSize(prep.alphabet), hs);
        let t = 0;
        let top = 0;
        let best = -1;
        for (let k = 0; k < hs.length; k++) {
          t += hs[k];
          if (hs[k] > best) {
            best = hs[k];
            top = lay.profile.letters[k];
          }
        }
        lay.totals[j] = t;
        lay.tops[j] = top;
      }
      lay.bars1d ??= new Pyramid({ factor: 2, pairing: 'independent', reducer: 'max' });
      lay.bars1d.setValues(1, C, lay.totals);
    }
    const profile = lay.profile!;
    const color = prep.channels.get('color', ColorFallback[prep.alphabet], ctx, Env[prep.alphabet]);
    lay.channel = color;
    const max = prep.scale === 'bits' ? Math.log2(alphabetSize(prep.alphabet)) : 1;
    const px = box.height / max;
    lay.pxPerUnit = px;
    const vis = axis.visible([]);
    if (!vis.length) return lay;
    const c0 = Math.max(0, vis[0]);
    const c1 = Math.min(held.columns, vis[vis.length - 1]);
    const bottom = box.y + box.height;
    const missing = tokens.colors['seq-missing'];
    // a letter's color is its class's, its oklab chroma raised to `seq-logo-chroma` (hue kept, a grey stays grey) and moved
    // in lightness to 3:1 against the track (the `ink` rule, the surface as the fill): once a color and table, never a letter a frame
    if (lay.glyphsVersion !== tokens.version) {
      lay.glyphs.clear();
      lay.glyphsVersion = tokens.version;
    }
    const surface = tokens.colors['seq-track'];
    const fillOf = (code: number, column: number) => {
      lay.datum.letter = String.fromCharCode(code);
      lay.datum.start = column;
      lay.datum.end = column + 1;
      const c = color.constant ?? color.at(lay.datum);
      const raw = typeof c === 'number' ? c : missing;
      let g = lay.glyphs.get(raw);
      if (g === undefined) lay.glyphs.set(raw, (g = glyphAgainst(raw, surface, tokens.logoChroma)));
      return g;
    };
    if (axis.devicePx < 2) {
      // one bar a device column: the tallest stack there, in its top letter's color
      const level = lay.bars1d!.level(0, pyramidLevelFor(1 / axis.devicePx, 2));
      const f = 2 ** level.columnLevel;
      for (let b = Math.floor(c0 / f); b * f < c1; b++) {
        const v = level.values[b];
        if (!(v > 0)) continue;
        let top = b * f;
        for (let j = b * f; j < Math.min(c1, (b + 1) * f); j++) if (lay.totals[j] === v) top = j;
        const x0 = axis.x(b * f);
        const x1 = axis.x(Math.min((b + 1) * f, held.columns));
        const k = lay.m++;
        if (4 * k + 4 > lay.bars.length) {
          lay.bars = grown(lay.bars, 4 * k + 4);
          lay.barColor = grown(lay.barColor, k + 1);
        }
        lay.bars[4 * k] = Math.min(x0, x1);
        lay.bars[4 * k + 1] = bottom - v * px;
        lay.bars[4 * k + 2] = Math.max(Math.abs(x1 - x0), 1 / axis.dpr);
        lay.bars[4 * k + 3] = v * px;
        lay.barColor[k] = fillOf(lay.tops[top], top);
      }
      return lay;
    }
    lay.cellW = axis.scale;
    lay.fontSize = Math.max(8, Math.min(64, Math.floor(axis.scale / config.monoAdvance)));
    lay.font = `${lay.fontSize}px mono`;
    const K = profile.letters.length;
    const hs = lay.heights.length >= K ? lay.heights : new Float32Array(K);
    const order = new Int32Array(K);
    for (let j = c0; j < c1; j++) {
      logoHeights(profile, j, prep.scale, prep.correction, alphabetSize(prep.alphabet), hs);
      for (let k = 0; k < K; k++) order[k] = k;
      // lowest first: the tallest letter on top
      order.sort((a, b) => hs[a] - hs[b]);
      if (!axis.cell(j)) continue;
      const x = (axis.sx0 + axis.sx1) / 2;
      let y = bottom;
      for (let q = 0; q < K; q++) {
        const k = order[q];
        const h = hs[k] * px;
        if (h < 0.5) continue;
        const i = lay.n++;
        if (i >= lay.lx.length) {
          lay.lx = grown(lay.lx, i + 1);
          lay.ly = grown(lay.ly, i + 1);
          lay.vs = grown(lay.vs, i + 1);
          lay.color = grown(lay.color, i + 1);
          lay.code = grown(lay.code, i + 1);
          lay.flat = grown(lay.flat, i + 1);
        }
        lay.lx[i] = x;
        lay.ly[i] = y;
        lay.vs[i] = h / (lay.fontSize * CapHeight);
        lay.flat[i] = h < MinLogoGlyphPx ? 1 : 0;
        lay.flatN += lay.flat[i];
        lay.code[i] = profile.letters[k];
        lay.color[i] = fillOf(profile.letters[k], j);
        y -= h;
      }
    }
    return lay;
  },

  draw(lay, pass: PassName, out: PassWriter, emphasis?: Emphasis) {
    drawPass(LogoDraws, lay, out, emphasis);
  },

  hit(lay, at): Hit | undefined {
    const x = at[0];
    const y = at[1];
    const { box } = lay;
    if (y < box.y || y >= box.y + box.height) return undefined;
    const position = lay.axis.positionAt(x);
    if (position < 0) return undefined;
    // the letter under the point: the stack's letter whose slot holds y
    let state: string | undefined;
    for (let i = 0; i < lay.n; i++) {
      if (!lay.axis.cell(position) || lay.lx[i] < lay.axis.sx0 || lay.lx[i] >= lay.axis.sx1) continue;
      const h = lay.vs[i] * lay.fontSize * CapHeight;
      if (y <= lay.ly[i] && y > lay.ly[i] - h) state = String.fromCharCode(lay.code[i]);
    }
    return { layer: lay.key, section: lay.section, position, system: lay.system, ...(state ? { state } : {}) };
  },
};

function alphabetSize(a: Alphabet): number {
  return a === 'protein' ? 20 : 4;
}

function drawLogo(lay: LogoLay, out: PassWriter) {
  const g = out.args;
  for (let k = 0; k < lay.m; k++) {
    g[0] = lay.bars[4 * k];
    g[1] = lay.bars[4 * k + 1];
    g[2] = lay.bars[4 * k + 2];
    g[3] = lay.bars[4 * k + 3];
    g[4] = lay.barColor[k];
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
  if (lay.pending) {
    g[0] = lay.box.x;
    g[1] = lay.box.y;
    g[2] = lay.box.width;
    g[3] = lay.box.height;
    g[4] = lay.tokens.colors['seq-pending'];
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
  if (!lay.n) return;
  // the stack's small letters, flat blocks of their color (never a pile of glyphs)
  for (let i = 0; i < lay.n; i++) {
    if (!lay.flat[i]) continue;
    const h = lay.vs[i] * lay.fontSize * CapHeight;
    g[0] = lay.lx[i] - lay.cellW / 2;
    g[1] = lay.ly[i] - h;
    g[2] = lay.cellW;
    g[3] = h;
    g[4] = lay.color[i];
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
  if (lay.flatN === lay.n) return;
  g[1] = lay.box.y + lay.box.height;
  out.beginGlyphsArgs(lay.font, NoColor);
  for (let i = 0; i < lay.n; i++) {
    if (lay.flat[i]) continue;
    g[0] = lay.lx[i];
    g[1] = lay.ly[i];
    g[2] = 0;
    g[3] = lay.vs[i];
    g[4] = lay.color[i];
    out.placedLetterArgs(lay.code[i]);
  }
  out.endGlyphs();
}

function drawHover(lay: LogoLay, out: PassWriter, emphasis: Emphasis | undefined) {
  hoverCell(lay, out, emphasis);
}

const LogoDraws = passDraws<LogoLay, Emphasis | undefined>({ data: drawLogo, hover: drawHover });
