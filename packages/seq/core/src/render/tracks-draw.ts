/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { AxisTransform } from './axis.js';
import {
  Anchor,
  type DrawLists,
  drawPass,
  NoColor,
  passDraws,
  Passes,
  type PassName,
  type PassWriter,
} from './draw-lists.js';
import type { Box, Emphasis } from './layer-type.js';
import type { SeqTokens } from './tokens.js';
import type { RenderConfig } from './config.js';
import { WarpKind } from './layout/strips.js';
import { TranslateWarp } from './layout/warp.js';
import { drawRail } from './rail.js';
import { wellClipArgs, wellCorners } from './wells.js';
import type { LayerLayout, SectionLayout, ViewLayout } from './tracks-types.js';
import { isExtensionName } from '@molstar/seqx-core/lang';
import { drawExtensionHover, drawSectionEmphasis, drawTints, emphasisPart, Part } from './tracks-emphasis.js';
import { hasHeaderActions } from './headers-draw.js';

// The view's draw (_spec/README.md "Render", Tracks layout): its passes, the background (the page and the wells), each
// section's chrome, rails and strips through their warps, the layers' marks clipped to their wells, and the breaks.

/** How far from a break's notch a point hits it, CSS px (fixed, as other hit boxes are). */
export const BreakSlop = 3;
/** A break's notch above the band: its width and height, CSS px. */
export const Notch = { width: 5, height: 3 };

/** Whether a layer type's name is an extension's, kept by name: the test allocates, a warm hover frame does not. */
const Extensions = new Map<string, boolean>();
function extensionType(name: string): boolean {
  let x = Extensions.get(name);
  if (x === undefined) Extensions.set(name, (x = isExtensionName(name)));
  return x;
}

/** Draws the layout's passes (all, or those named, e.g. the emphasis alone) into the lists, cleared first. */
export function drawView(
  layout: ViewLayout,
  lists: DrawLists,
  emphasis?: Emphasis,
  passes: readonly PassName[] = Passes,
): void {
  lists.reset(passes);
  lists.width = layout.width;
  lists.height = layout.height;
  // indexed loops: an iterator a loop is an allocation where the compiler keeps it (nothing allocated once warm)
  for (let pi = 0; pi < passes.length; pi++) {
    const pass = passes[pi];
    drawPass(ViewDraws, layout, lists.writer(pass), emphasis);
  }
}

/** The view's own draw of each pass, through `drawPass` like a layer type's. */
const ViewDraws = passDraws<ViewLayout, Emphasis | undefined>({
  background: drawBackground,
  data: drawTracks,
  window: drawTracks,
  highlight: drawTracks,
  selection: drawTracks,
  hover: drawTracks,
  cursor: drawTracks,
  labels: drawTracks,
  tip: drawTracks,
});

/**
 * The background pass: the page, each track's well (a layout that fills its bands: a squircle by the squircle rule,
 * `trackShape`), and with the tints under the data (`emphasis.selectionPlace: 'under'`) the highlight channels' tints
 * and then the selection's above them, a band a section or a box a track (`emphasis.selectionShape`: inside the well's
 * clip), a strip at a time through its warp.
 */
function drawBackground(layout: ViewLayout, out: PassWriter, emphasis: Emphasis | undefined) {
  const { tokens, config } = layout;
  out.fill(0, 0, layout.width, layout.height, tokens.colors['msui-bg']);
  const under = config.emphasis.selectionPlace === 'under' && emphasis !== undefined;
  const column = config.emphasis.selectionShape === 'column';
  for (let si = 0; si < layout.sections.length; si++) {
    const s = layout.sections[si];
    if (s.layout.fillBands) {
      for (let ti = 0; ti < s.tracks.length; ti++) {
        const t = s.tracks[ti];
        const { box, well } = t;
        if (box.y >= layout.height) continue;
        // the well through `args` (nothing boxed): its squircle by its height
        const g = out.args;
        g[0] = well.x;
        g[1] = well.y;
        g[2] = well.width;
        g[3] = well.height;
        wellCorners(g, 7, config);
        g[4] = tokens.colors['seq-track'];
        g[5] = NoColor;
        g[6] = 0;
        out.rectArgs();
        if (under && !column && s.axis) {
          const strips = s.placement.strips;
          if (strips.n) {
            // in the page's frame: the first strip's (a layout that fills its bands translates it by nothing)
            wellClipArgs(out.args, s, 0, ti, box);
            out.args[0] += strips.dx[0];
            out.args[1] += strips.dy[0];
            out.clipSquircleArgs();
            emphasisPart(Part.tints);
          }
          drawTints(out, s, s.axis, box, emphasis!, tokens, config);
          if (strips.n) {
            out.unclip();
            emphasisPart(Part.outlines);
            drawTints(out, s, s.axis, box, emphasis!, tokens, config);
            emphasisPart(Part.all);
          }
        }
      }
      if (under && column && s.axis && s.band) drawTints(out, s, s.axis, s.band, emphasis!, tokens, config);
      continue;
    }
    if (!under) continue;
    const strips = s.placement.strips;
    const pb = s.placement.box;
    if (s.placement.clip === true) clipTo(out, pb);
    for (let i = 0; i < strips.n; i++) {
      const w = openStrip(s, i, out);
      const a = strips.axes[i];
      const col = s.columns[i];
      if (column && col) drawTints(w, s, a, col, emphasis!, tokens, config);
      else if (!column)
        for (let j = 0; j < strips.bands; j++) {
          const box = s.boxes[i * strips.bands + j];
          if (box.height > 0) drawTints(w, s, a, box, emphasis!, tokens, config);
        }
      closeStrip(s, i, out);
    }
    if (s.placement.clip === true) out.unclip();
  }
}

/** The strips' translation (one at a time: drawing is not reentrant). */
const Translate = new TranslateWarp();

/** The writer a strip draws into for `out`'s pass: its polar warp's, else `out` with the translation marked. */
/** A clip to a box through the writer's `args` (its fractional numbers passed to a call would be boxed). */
function clipTo(out: PassWriter, b: Box): void {
  const g = out.args;
  g[0] = b.x;
  g[1] = b.y;
  g[2] = b.width;
  g[3] = b.height;
  out.clipArgs();
}

function openStrip(s: SectionLayout, i: number, out: PassWriter): PassWriter {
  const strips = s.placement.strips;
  if (strips.warp[i] === WarpKind.polar) return s.placement.warpOf!(i).writer(out);
  Translate.open(out);
  return out;
}

/** What a strip drew, onto the page through its warp. */
function closeStrip(s: SectionLayout, i: number, out: PassWriter): void {
  const strips = s.placement.strips;
  if (strips.warp[i] === WarpKind.polar) s.placement.warpOf!(i).close(out);
  else Translate.close(out, strips.dx, strips.dy, i);
}

/**
 * Any other pass (the writer's): per section its layout's chrome, then per strip (through its warp) its emphasis marks
 * of the pass over its band column (`column`) or each band (`tracks`), each band's chrome (its label) and its layers'
 * marks, and its breaks.
 */
function drawTracks(layout: ViewLayout, out: PassWriter, emphasis: Emphasis | undefined) {
  const { config } = layout;
  const pass = out.pass.name;
  const column = config.emphasis.selectionShape === 'column';
  for (let si = 0; si < layout.sections.length; si++) {
    const s = layout.sections[si];
    // the tip is a rows layer's (a hovered cell read whole) and the headers' actions: a section with neither has nothing
    // to walk, as the window pass
    if (pass === 'tip' && !hasRows(s) && s.layout.name !== 'rings' && !hasHeaderActions(s)) continue;
    s.layout.chrome(s, pass, out, emphasis, -1, -1);
    if (pass === 'window' && s.rails)
      for (let i = 0; i < s.rails.length; i++) drawRail(s.rails[i], out, emphasis, s.tokens, s.config);
    if (pass === 'window' && s.rowsRails)
      for (let i = 0; i < s.rowsRails.length; i++) drawRail(s.rowsRails[i], out, emphasis, s.tokens, s.config);
    if (pass === 'window' && s.lanesRails)
      for (let i = 0; i < s.lanesRails.length; i++) drawRail(s.lanesRails[i], out, emphasis, s.tokens, s.config);
    const strips = s.placement.strips;
    if (!strips.n) {
      // pending: its tracks' chrome alone
      for (let ti = 0; ti < s.tracks.length; ti++)
        if (s.tracks[ti].box.y < layout.height) s.layout.chrome(s, pass, out, emphasis, 0, ti);
      continue;
    }
    drawStrips(layout, s, out, emphasis, pass, column);
    // what crosses strips, on the page; a section that scrolls inside its box is clipped to it
    const sl = s.sectionLays;
    if (!sl.length) continue;
    const clipped = s.placement.clip === true;
    if (clipped) clipTo(out, s.placement.box);
    for (let li = 0; li < sl.length; li++) sl[li].type.draw(sl[li].lay, pass, out, emphasis);
    if (clipped) out.unclip();
  }
  // the bridges' ribbons, on the page between their sides (_spec/README.md "Render", Bridges)
  layout.bridges?.draw(pass, out, emphasis);
}

/** Whether a section has a rows layer (kept by the layout: its layers do not change within it). */
const RowsIn = new WeakMap<SectionLayout, boolean>();
function hasRows(s: SectionLayout): boolean {
  let has = RowsIn.get(s);
  if (has === undefined) {
    has = false;
    for (const band of s.lays) for (const l of band ?? []) if (l.type.name === 'rows') has = true;
    RowsIn.set(s, has);
  }
  return has;
}

/**
 * Whether the hover is a letters layer's lifted through a codon map on this section (a residue over its codon): only
 * that lights the hover's own section, on its own-system letters tracks (axis.md, An item split by a path). An item
 * through a path lights its own outline and connector, nothing else.
 */
function liftedLetters(s: SectionLayout, emphasis: Emphasis): boolean {
  const h = emphasis.hover;
  if (!h || !h.extent || !h.extent.parts.length || h.extent.parts[0].system === h.system) return false;
  for (const band of s.lays)
    if (band)
      for (let i = 0; i < band.length; i++)
        if (band[i].key === h.layer)
          return band[i].type.name === 'letters' && (band[i].lay as { through?: string }).through === 'codon';
  return false;
}

/** Whether a track's box holds a letters layer on the section's own system (not the hovered one): the DNA track. */
function ownLettersIn(lays: readonly LayerLayout[] | undefined, emphasis: Emphasis): boolean {
  if (lays)
    for (let i = 0; i < lays.length; i++)
      if (
        lays[i].type.name === 'letters' &&
        lays[i].key !== emphasis.hover?.layer &&
        (lays[i].lay as { through?: string }).through === 'own'
      )
        return true;
  return false;
}

/** A section's strips of a pass, each through its warp (`drawTracks`). */
function drawStrips(
  layout: ViewLayout,
  s: SectionLayout,
  out: PassWriter,
  emphasis: Emphasis | undefined,
  pass: PassName,
  column: boolean,
) {
  const { tokens } = layout;
  const lift = pass === 'hover' && emphasis !== undefined && liftedLetters(s, emphasis);
  const strips = s.placement.strips;
  const nb = strips.bands;
  // a section that scrolls inside its box: each strip clipped to it, in the strip's frame (its warp moves the clip)
  const clipped = s.placement.clip === true;
  const pb = s.placement.box;
  {
    for (let i = 0; i < strips.n; i++) {
      const w = openStrip(s, i, out);
      const a = strips.axes[i];
      const col = s.columns[i];
      const dy = strips.dy[i];
      if (clipped) {
        const g = w.args;
        g[0] = pb.x - strips.dx[i];
        g[1] = pb.y - dy;
        g[2] = pb.width;
        g[3] = pb.height;
        w.clipArgs();
      }
      if (column && emphasis && col) drawSectionEmphasis(w, s, a, col, emphasis, layout, true);
      // the lifted hover is no column: its codon over each own-system letters track alone
      if (column && lift && emphasis)
        for (let j = 0; j < nb; j++) {
          const k = i * nb + j;
          if (s.boxes[k].height > 0 && ownLettersIn(s.lays[k], emphasis))
            drawSectionEmphasis(w, s, a, s.boxes[k], emphasis, layout, false);
        }
      for (let j = 0; j < nb; j++) {
        const k = i * nb + j;
        const box = s.boxes[k];
        if (box.y + dy >= layout.height) continue;
        if (box.height <= 0 && !s.layout.fillBands) continue;
        s.layout.chrome(s, pass, w, emphasis, i, j);
        if (pass !== 'labels' && !column && emphasis) {
          // a box a track: its tints inside its well (the well's shape), the outlines around marks unclipped
          const tint =
            s.layout.fillBands &&
            ((pass === 'highlight' && emphasis.highlights !== undefined) ||
              (pass === 'selection' &&
                (emphasis.selection !== undefined || emphasis.projected?.selection !== undefined)));
          if (tint) {
            // the tints inside the well's clip, then the outlines over them unclipped
            wellClipArgs(w.args, s, i, j, box);
            w.clipSquircleArgs();
            emphasisPart(Part.tints);
            drawSectionEmphasis(w, s, a, box, emphasis, layout, !(lift && ownLettersIn(s.lays[k], emphasis)));
            w.unclip();
            emphasisPart(Part.outlines);
            drawSectionEmphasis(w, s, a, box, emphasis, layout, !(lift && ownLettersIn(s.lays[k], emphasis)));
            emphasisPart(Part.all);
          } else drawSectionEmphasis(w, s, a, box, emphasis, layout, !(lift && ownLettersIn(s.lays[k], emphasis)));
        }
        const lays = s.lays[k];
        // the window pass is the layouts' own (their chrome): no layer draws there, so a hover's redraw walks none
        if (!lays || !lays.length || pass === 'window') continue;
        // the data clipped to its well: one clip a track (its squircle in a layout that fills its bands)
        const clip = pass === 'data';
        if (clip) {
          wellClipArgs(w.args, s, i, j, box);
          w.clipSquircleArgs();
        }
        for (let li = 0; li < lays.length; li++) {
          const l = lays[li];
          if (l.zoomIn && pass === 'data') {
            drawZoomIn(w, box, tokens, layout.config);
            continue;
          }
          l.type.draw(l.lay, pass, w, emphasis);
          // an extension's layer type is given the hover's look by the view (Emphasis), unless it draws its own
          if (pass === 'hover' && emphasis?.hover?.layer === l.key && !l.type.ownHover && extensionType(l.type.name))
            drawExtensionHover(w, l.type, l.lay, emphasis.hover, a, layout.tokens, layout.config);
        }
        if (clip) w.unclip();
      }
      const breaks = s.stripBreaks[i];
      if (pass === 'data' && col && breaks.length) drawBreaks(w, breaks, a, col, tokens);
      if (clipped) w.unclip();
      closeStrip(s, i, out);
    }
  }
}

/**
 * A strip's breaks (_spec/README.md "Render", Tracks layout): a hairline a device pixel wide in the break ink across
 * the band column, over the cells either side (they touch), and a notch above it. Zero units wide: the hairline has no
 * width of its own on the axis.
 */
function drawBreaks(out: PassWriter, b: readonly number[], axis: AxisTransform, band: Box, tokens: SeqTokens) {
  // through `args` (a number passed to a call that is not inlined is boxed: nothing allocated once warm)
  const g = out.args;
  const dpr = axis.dpr;
  const px = 1 / dpr;
  for (let i = 0; i + 2 < b.length; i += 3) {
    const x = Math.round(b[i] * dpr) / dpr;
    for (let k = 0; k < 2; k++) {
      g[0] = k === 0 ? x - px / 2 : x - Notch.width / 2;
      g[1] = k === 0 ? band.y : band.y - Notch.height;
      g[2] = k === 0 ? px : Notch.width;
      g[3] = k === 0 ? band.height : Notch.height;
      g[4] = tokens.colors['seq-break'];
      g[5] = NoColor;
      g[6] = 0;
      g[7] = 0;
      g[8] = 0;
      out.rectArgs();
    }
  }
}

/** The text a layer refused by the store's bounds says over its band (_spec/README.md "Render", Data). */
export const ZoomInText = 'zoom in to read this track';

/**
 * A refused layer's band: the pending hatch (stripes of the pending ink, as a band read nothing shows) with the text
 * "zoom in to read this track" in the middle, the label's font and the muted ink: distinct from pending, which settles
 * by itself.
 */
function drawZoomIn(w: PassWriter, box: Box, tokens: SeqTokens, config: RenderConfig): void {
  const g = w.args;
  const pending = tokens.colors['seq-pending'];
  // the hatch: a stripe every 8 px across the band, 4 px wide
  for (let x = box.x; x < box.x + box.width; x += 8) {
    g[0] = x;
    g[1] = box.y;
    g[2] = Math.min(4, box.x + box.width - x);
    g[3] = box.height;
    g[4] = pending;
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    w.rectArgs();
  }
  g[0] = box.x + box.width / 2;
  g[1] = box.y + box.height / 2 + 4;
  g[2] = 0;
  g[3] = box.width;
  g[4] = tokens.colors['msui-text-muted'];
  w.textArgs(ZoomInText, Anchor.middle, config.labelFont);
}
