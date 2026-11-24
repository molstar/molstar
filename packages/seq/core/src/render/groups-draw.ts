/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Anchor, type PassWriter } from './draw-lists.js';
import type { Hit } from './layer-type.js';
import { ChevronPx, CountGapPx } from './groups.js';
import type { SectionLayout } from './tracks-types.js';

// A section's group headers in the label column (_spec/README.md "Render", Track groups): the chevron (▾ open, ▸
// folded), the label and a folded group's quiet count, drawn from what the layout measured (a draw measures nothing);
// the header's hit, `group-header`, its group's address in `layer`.

const AnchorStart: number = Anchor.start;

/** Every group header of a section, in the labels pass; nothing drawn across the data. */
export function drawGroupHeaders(out: PassWriter, s: SectionLayout): void {
  const headers = s.groups;
  if (!headers) return;
  const g = out.args;
  const { config, tokens } = s;
  for (let i = 0; i < headers.length; i++) {
    const h = headers[i];
    if (h.box.y >= s.pageHeight) continue;
    const ink = tokens.colors['msui-text-muted'];
    // the chevron is chrome: a figure draws none, its label at the column's start
    const x = h.chevron ? h.box.x + ChevronPx : h.box.x;
    if (h.chevron) {
      g[0] = h.box.x;
      g[1] = h.baseline;
      g[2] = ChevronPx;
      g[3] = ChevronPx;
      g[4] = ink;
      out.textArgs(h.folded ? '▸' : '▾', AnchorStart, config.labelFont);
    }
    g[0] = x;
    g[1] = h.baseline;
    g[2] = h.labelWidth;
    g[3] = h.labelWidth;
    g[4] = tokens.colors['msui-text'];
    if (h.label) out.textArgs(h.label, AnchorStart, config.labelFont);
    if (!h.count || !h.countWidth) continue;
    g[0] = x + h.labelWidth + CountGapPx;
    g[1] = h.baseline;
    g[2] = h.countWidth;
    g[3] = h.countWidth;
    g[4] = ink;
    out.textArgs(h.count, AnchorStart, config.labelFont);
  }
}

/** The group header under a point (its line, or a folded row's label column): `group-header`, its address in `layer`. */
export function hitGroupHeader(s: SectionLayout, x: number, y: number): Hit | undefined {
  const headers = s.groups!;
  for (let i = 0; i < headers.length; i++) {
    const b = headers[i].box;
    if (x >= b.x && x < b.x + b.width && y >= b.y && y < b.y + b.height) return headers[i].hit;
  }
  return undefined;
}
