/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hashFnv32a } from '@molstar/core/data/util/hash-functions';
import {
  Anchor,
  type DrawLists,
  NoColor,
  PathCommandPoints,
  type PassName,
  Passes,
  type PassVisitor,
  walkPass,
} from './draw-lists.js';
import { colorHex } from './tokens.js';

// The recording backend (_spec/README.md "Render", Backends): one line per primitive, pass by pass, numbers to two
// decimals, colors `#rrggbbaa`: the format of the committed records, read as diffs.
//
//   pass data
//   clip 136 8 656 20
//   rect 136 8 15.5 20 #cfe0f5ff
//   rect 140 8 16 20 none stroke #2d271eff 1.5
//   rect 140 30 80 10 #7a8fa6ff squircle 3.33 left
//   glyphs "10px mono" 22.5 per M 143.75 #2d271eff, E 159.25 #2d271eff
//   glyphs "10px mono" placed #2d271eff M 143.75 22.5 90 1, E 159.25 22.5 92 1   (x, y, degrees, vertical scale)
//   text "Domains" 8 30 start "12px sans" #6b655cff width 47.34 max 112
//   unclip

/** A number to at most two decimals, no trailing zeros, no `-0`. */
export function fmt(n: number): string {
  if (n === Infinity) return 'inf';
  if (n === -Infinity) return '-inf';
  const r = Math.round(n * 100) / 100;
  return Object.is(r, -0) ? '0' : String(r);
}

const AnchorNames = ['start', 'middle', 'end'];
const CommandNames = ['M', 'L', 'C', 'Z'];

/** A fill: `#rrggbbaa`, or `none`. */
function paint(c: number): string {
  return c === NoColor ? 'none' : colorHex(c);
}

function quote(s: string): string {
  return JSON.stringify(s);
}

/** A squircle's radius and its round ends (`left`, `right`, both unsaid) and cut edges (`cut-top`, `cut-bottom`). */
function roundOf(radius: number, round: number): string {
  if (!(radius > 0)) return '';
  let s = ` squircle ${fmt(radius)}`;
  if ((round & 3) === 1) s += ' left';
  else if ((round & 3) === 2) s += ' right';
  if (round & 4) s += ' cut-top';
  if (round & 8) s += ' cut-bottom';
  return s;
}

/** The draw lists as lines, pass by pass (every pass, or those named); an empty pass is its header alone. */
export function recordDrawLists(lists: DrawLists, passes: readonly PassName[] = Passes): string {
  const lines: string[] = [`size ${fmt(lists.width)} ${fmt(lists.height)}`];
  const v: PassVisitor = {
    rect(l, i) {
      let s = `rect ${fmt(l.x[i])} ${fmt(l.y[i])} ${fmt(l.w[i])} ${fmt(l.h[i])} ${paint(l.fill[i])}`;
      if (l.stroke[i] !== NoColor) s += ` stroke ${colorHex(l.stroke[i])} ${fmt(l.strokeWidth[i])}`;
      s += roundOf(l.radius[i], l.round[i]);
      lines.push(s);
    },
    path(l, i) {
      const parts: string[] = ['path'];
      let pt = l.pointStart[i];
      for (let c = l.cmdStart[i], end = c + l.cmdCount[i]; c < end; c++) {
        const cmd = l.commands[c];
        parts.push(CommandNames[cmd]);
        for (let k = 0; k < PathCommandPoints[cmd]; k++, pt++)
          parts.push(fmt(l.points[2 * pt]), fmt(l.points[2 * pt + 1]));
      }
      parts.push(paint(l.fill[i]));
      if (l.stroke[i] !== NoColor) parts.push('stroke', colorHex(l.stroke[i]), fmt(l.strokeWidth[i]));
      if (l.dash[i] > 0) parts.push('dash', fmt(l.dash[i]));
      lines.push(parts.join(' '));
    },
    glyphs(l, i) {
      const one = l.color[i] !== NoColor;
      const placed = l.placed[i] === 1;
      const letters: string[] = [];
      for (let j = l.start[i], end = j + l.count[i]; j < end; j++) {
        let s = `${String.fromCharCode(l.codes[j])} ${fmt(l.x[j])}`;
        // a placed letter: its own baseline, angle (degrees) and vertical scale
        if (placed) s += ` ${fmt(l.ly[j])} ${fmt((l.angle[j] * 180) / Math.PI)} ${fmt(l.vs[j])}`;
        letters.push(one ? s : `${s} ${colorHex(l.colors[j])}`);
      }
      const color = one ? colorHex(l.color[i]) : 'per';
      const y = placed ? 'placed' : fmt(l.y[i]);
      lines.push(`glyphs ${quote(lists.string(l.font[i]))} ${y} ${color} ${letters.join(', ')}`);
    },
    text(l, i) {
      let s = `text ${quote(lists.string(l.str[i]))} ${fmt(l.x[i])} ${fmt(l.y[i])} ${AnchorNames[l.anchor[i] ?? Anchor.start]}`;
      s += ` ${quote(lists.string(l.font[i]))} ${colorHex(l.color[i])} width ${fmt(l.width[i])}`;
      if (l.maxWidth[i] !== Infinity) s += ` max ${fmt(l.maxWidth[i])}`;
      lines.push(s);
    },
    image(l, i) {
      const img = l.images[i];
      const hash = (hashFnv32a(img.data) >>> 0).toString(16).padStart(8, '0');
      lines.push(
        `image ${fmt(l.x[i])} ${fmt(l.y[i])} ${fmt(l.w[i])} ${fmt(l.h[i])} ${img.width}x${img.height} ${hash}`,
      );
    },
    clip(l, i) {
      lines.push(`clip ${fmt(l.x[i])} ${fmt(l.y[i])} ${fmt(l.w[i])} ${fmt(l.h[i])}${roundOf(l.radius[i], l.round[i])}`);
    },
    unclip() {
      lines.push('unclip');
    },
  };
  for (const name of passes) {
    lines.push(`pass ${name}`);
    walkPass(lists.passes[name], v);
  }
  return `${lines.join('\n')}\n`;
}
