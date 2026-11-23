/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Emphasis, LayerType, Location, PassName, PassWriter } from '@molstar/seqx-core/extension';
import { PairsReader } from '@molstar/seqx-core/extension';
import { DefaultRondoConfig, type RondoConfig } from './config.js';
import { Kind, layoutRondo, RondoLay, type RondoPrep, stemAt } from './lay.js';
import { polygonAt } from './shapes.js';
import { announceRondo, describeRondo, summaryRondo } from './words.js';

// `rondo:ribbons` (_spec/README.md): a `scope: 'section'` layer over core's pairs reader, on a ring (through its warp)
// and a row of blocks. Its data pass draws the polygons the layout made; its hover a ribbon lit and outlined, a stem's
// or a block's ribbons lit over a wash of the rest; its selection the ribbons with an end in it in the accent; its
// highlight channels the ribbons touching their Location. A reader without it draws `pairs` (its `fallback`).

/** No colour: an absent stroke (the draw lists' `NoColor`). */
const NoColor = 0;
/** A polygon's points copied for the writer (`pathArgs` reads from the start): grown, kept. */
let Scratch = new Float64Array(256);
/** `0xRRGGBBAA` at an opacity. */
const alpha = (c: number, a: number) => ((c & 0xffffff00) | Math.round(Math.max(0, Math.min(1, a)) * 255)) >>> 0;

/**
 * A polygon's paint, set by the caller before `drawPoly`: fill, stroke, stroke width (a colour past 2^31 passed to a
 * call that is not inlined is boxed: the emphasis passes paint through here, nothing allocated once warm).
 */
const Paint = new Float64Array(3);

function drawPoly(lay: RondoLay, out: PassWriter, i: number) {
  const P = lay.polys;
  const n = P.len[i];
  if (n < 3) return;
  if (Scratch.length < 2 * n) Scratch = new Float64Array(2 * Math.max(2 * n, Scratch.length));
  const s = 2 * P.start[i];
  for (let k = 0; k < 2 * n; k++) Scratch[k] = P.xy[s + k];
  const a = out.args;
  a[4] = Paint[0];
  a[5] = Paint[1];
  a[6] = Paint[2];
  a[7] = 0;
  out.pathArgs(Scratch, n, true);
}

/** Sets the paint: a fill, no stroke. */
function fill(c: number) {
  Paint[0] = c;
  Paint[1] = NoColor;
  Paint[2] = 0;
}

/** Whether a Location covers any of a stem's positions on its system. */
function touches(lay: RondoLay, loc: Location | undefined, stem: number): boolean {
  if (!loc) return false;
  const sys = lay.systems[lay.stemSystem[stem]];
  const s = lay.stemStart[stem];
  const e = lay.stemEnd[stem];
  for (let p = 0; p < loc.parts.length; p++) {
    const part = loc.parts[p];
    if (part.system !== sys) continue;
    const r = part.ranges;
    for (let k = 0; k + 1 < r.length; k += 2) if (r[k] < e && r[k + 1] > s) return true;
  }
  return false;
}

/** The ribbons with a stem `lit(stem)` at full opacity over a wash of the hole (`rondo.fade`). */
function washed(lay: RondoLay, out: PassWriter, fade: number, lit: (stem: number) => boolean, color = -1) {
  if (fade > 0 && lay.hole >= 0) {
    fill(alpha(lay.wash, fade));
    drawPoly(lay, out, lay.hole);
  }
  const P = lay.polys;
  for (let i = 0; i < P.n; i++) {
    const kind = lay.polyKind[i];
    const ref = lay.polyRef[i];
    let on = false;
    if (kind === Kind.ribbon) on = lit(lay.ribA[ref]) || lit(lay.ribB[ref]);
    else if (kind === Kind.trumpet || kind === Kind.fan) on = lit(ref);
    else if (kind === Kind.stub) on = lit(lay.stubStem[ref]);
    if (!on) continue;
    fill(color >= 0 ? alpha(color, 0.75) : alpha(lay.polyColor[i], 1));
    drawPoly(lay, out, i);
  }
}

/** The hover's state, read by the lit test (no closure made a frame); cleared after each emphasis draw. */
const Lit = { lay: undefined as unknown as RondoLay, stem: -1, system: '' };

/** Whether the module holds a lay between draws (a test: none is kept). */
export const emphasisHeld = (): boolean => Lit.lay !== undefined || Loc.loc !== undefined;
const litStem = (s: number) => s === Lit.stem;
const litSystem = (s: number) => Lit.lay.systems[Lit.lay.stemSystem[s]] === Lit.system;
const Loc = { loc: undefined as Location | undefined };
const litLocation = (s: number) => touches(Lit.lay, Loc.loc, s);

/** A polygon at full opacity. */
function lit(lay: RondoLay, out: PassWriter, i: number) {
  if (i < 0) return;
  fill(alpha(lay.polyColor[i], 1));
  drawPoly(lay, out, i);
}

/** A polygon at full opacity outlined in `seq-hover`. */
function outlined(lay: RondoLay, out: PassWriter, i: number) {
  Paint[0] = alpha(lay.polyColor[i], 1);
  Paint[1] = lay.hover;
  Paint[2] = 1;
  drawPoly(lay, out, i);
}

function drawHover(lay: RondoLay, out: PassWriter, e: Emphasis, cfg: RondoConfig) {
  const h = e.hover;
  if (!h) return;
  Lit.lay = lay;
  if (h.layer === lay.key && h.element) {
    const i = lay.byKey.get(h.element.key as string);
    if (i === undefined) return;
    const kind = lay.polyKind[i];
    if (kind === Kind.trumpet || kind === Kind.fan) {
      // a stem: its ribbons lit, the rest of the hole washed (d3-chord's fade)
      Lit.stem = lay.polyRef[i];
      washed(lay, out, cfg.fade, litStem);
      outlined(lay, out, i);
      return;
    }
    // a ribbon or a stub: at full opacity, outlined, its two stems' mouths lit
    const a = kind === Kind.ribbon ? lay.ribA[lay.polyRef[i]] : lay.stubStem[lay.polyRef[i]];
    lit(lay, out, lay.stemPoly[a]);
    if (kind === Kind.ribbon) lit(lay, out, lay.stemPoly[lay.ribB[lay.polyRef[i]]]);
    outlined(lay, out, i);
    return;
  }
  // another layer's hover on the section, the same wash over what it touches: an item (a band) the ribbons with an end
  // in its extent; a place (the hue ring, the bare axis) its position's stem; a block's arc or label (no position) every
  // ribbon touching the block
  if (h.section !== lay.sectionId || h.layer === lay.key || !h.system) return;
  if (h.extent) {
    Loc.loc = h.extent;
    washed(lay, out, cfg.fade, litLocation);
    return;
  }
  if (h.position !== undefined) {
    const s = stemAt(lay, h.system, h.position);
    if (s < 0) return;
    Lit.stem = s;
    washed(lay, out, cfg.fade, litStem);
    return;
  }
  let any = false;
  for (let s = 0; s < lay.ns && !any; s++) any = lay.stemPlaced[s] === 1 && lay.systems[lay.stemSystem[s]] === h.system;
  if (!any) return;
  Lit.system = h.system;
  washed(lay, out, cfg.fade, litSystem);
}

/** The emphasis passes: the hover, the selection, the highlights. */
function emphasized(lay: RondoLay, pass: PassName, out: PassWriter, emphasis: Emphasis, cfg: RondoConfig) {
  if (pass === 'hover') drawHover(lay, out, emphasis, cfg);
  else if (pass === 'selection' && emphasis.selection) {
    Loc.loc = emphasis.selection;
    washed(lay, out, 0, litLocation, lay.accent);
  } else if (pass === 'highlight' && emphasis.highlights)
    for (let h = 0; h < emphasis.highlights.length; h++) {
      Loc.loc = emphasis.highlights[h].location;
      washed(lay, out, h === 0 ? cfg.fade : 0, litLocation, emphasis.highlights[h].color);
    }
}

export interface RibbonsOptions {
  readonly config?: RondoConfig;
}

/** The ribbons' layer type with its settings (`rondo.*`): `rondo` composes it with the defaults. */
export function ribbonsType(config: RondoConfig = DefaultRondoConfig): LayerType<RondoPrep, RondoLay> {
  return {
    name: 'rondo:ribbons',
    shape: 'pairs',
    fallback: 'pairs',
    scope: 'section',
    ownHover: true,
    prepare(layer, model) {
      const p = layer.params;
      const theme = JSON.stringify(layer.theme.color ?? '');
      return {
        source: PairsReader.source(layer, model),
        model,
        stems: p.stems === 'none' ? 'none' : 'level',
        bundle: typeof p.bundle === 'number' ? Math.max(0, Math.min(1, p.bundle)) : 0.7,
        stubs: p.stubs === 'none' ? 'none' : 'partner',
        ...(typeof p.opens === 'string' ? { opens: p.opens } : {}),
        hues: theme.includes('rondo:hues'),
        config,
      };
    },
    layout(prep, ctx, data, _overlay, prev) {
      const lay = prev ?? new RondoLay();
      layoutRondo(lay, prep, ctx, data);
      if (lay.empty === 'layout')
        ctx.report({
          code: 'rondo.ribbons.layout',
          severity: 'warning',
          message: 'Ribbons are drawn on rings and rows of blocks: not on this layout',
          address: ctx.key,
        } as never);
      return lay;
    },
    draw(lay, pass: PassName, out, emphasis) {
      const cfg = lay.prep?.config ?? config;
      if (pass === 'data') {
        for (let j = 0; j < lay.polys.n; j++) {
          const i = lay.order[j];
          if (i === lay.hole) continue;
          fill(lay.polyColor[i]);
          drawPoly(lay, out, i);
        }
        return;
      }
      if (pass === 'labels') {
        const a = out.args;
        for (let s = 0; s < lay.nb; s++) {
          const t = lay.stubLabel[s];
          if (!t) continue;
          a[0] = lay.stubLabelX[s];
          a[1] = lay.stubLabelY[s] + 4;
          a[2] = lay.stubLabelW[s];
          a[3] = Number.POSITIVE_INFINITY;
          a[4] = lay.ink;
          out.textArgs(t, 1, lay.font);
        }
        return;
      }
      if (!emphasis) return;
      Lit.lay = lay;
      try {
        emphasized(lay, pass, out, emphasis, cfg);
      } finally {
        // no lay or Location kept past the draw (a viewer disposed is not held by this module)
        Lit.lay = undefined as unknown as RondoLay;
        Loc.loc = undefined;
      }
    },
    hit(lay, at) {
      const i = polygonAt(lay, at[0], at[1]);
      return i < 0 ? undefined : lay.hits[i];
    },
    filtered: (lay) => lay.read?.filtered ?? 0,
    describe: describeRondo,
    summary: summaryRondo,
    announce: announceRondo,
  };
}

/** `rondo:ribbons` with the default settings. */
export const RibbonsType = ribbonsType();
