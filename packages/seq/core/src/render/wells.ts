/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A track's well (_spec/README.md "Render", Shapes and Tracks layout): its data's band grown by its padding, a squircle
// by the squircle rule over the well's height (`trackShape: 'square'`: a box), filled in the background and the clip of
// its data. Its numbers go through `args` (a number passed to a call that is not inlined is boxed: nothing allocated
// once warm).

import type { RenderConfig } from './config.js';
import { Round } from './draw-lists.js';
import type { Box } from './layer-type.js';
import type { SectionLayout } from './tracks.js';

/**
 * The well's squircle radius and round ends into `g[r]`, `g[r + 1]` for the well `g[2]` × `g[3]` px: the squircle
 * rule (a third of its height, capped; none below `squircleMin` in either dimension), round at both ends; 0 for a
 * square well (`trackShape: 'square'`).
 */
export function wellCorners(g: Float64Array, r: number, config: RenderConfig): void {
  const w = g[2];
  const h = g[3];
  if (config.trackShape !== 'squircle' || w < config.squircleMin || h < config.squircleMin) {
    g[r] = 0;
    g[r + 1] = 0;
    return;
  }
  g[r] = Math.min(h * config.squircleRatio, config.squircleMax, w / 2, h / 2);
  g[r + 1] = Round.both;
}

/**
 * The clip of band `j` on strip `i` into `g` (x, y, w, h, radius, round ends: `PassWriter.clipSquircleArgs`), in the
 * strip's frame: the band's well (its padding round it) over `box`'s x. In a layout that fills its bands (tracks, a
 * grid) it is the well's squircle (its radius by the track's well), round where `box` reaches the track's well's ends
 * (a block within the row: square where it stops short of them); a grid's cell (a pair of blocks) is square too on
 * an edge inside the plot (`cutTop`, `cutBottom`), so only the plot's outer corners round. Elsewhere (text, rings: their
 * wells are boxes) the well's box.
 */
export function wellClipArgs(g: Float64Array, s: SectionLayout, i: number, j: number, box: Box): void {
  const strips = s.placement.strips;
  const k = i * strips.bands + j;
  const pt = strips.padT[k];
  g[0] = box.x;
  g[1] = strips.bandY[k] - pt;
  g[3] = strips.bandH[k] + pt + strips.padB[k];
  g[4] = 0;
  g[5] = 0;
  const well = s.layout.fillBands ? s.tracks[j]?.well : undefined;
  if (well) {
    // the radius by the track's well, the clip's x the box's
    const h = g[3];
    g[2] = well.width;
    g[3] = well.height;
    wellCorners(g, 4, s.config);
    g[3] = h;
    if (g[4] > 0) {
      let round =
        (box.x <= well.x + 0.5 ? Round.left : 0) | (box.x + box.width >= well.x + well.width - 0.5 ? Round.right : 0);
      // a grid's cells share one well (the plot's): an edge inside it is cut
      if (s.layout.axesY) {
        if (g[1] > well.y + 0.5) round |= Round.cutTop;
        if (g[1] + h < well.y + well.height - 0.5) round |= Round.cutBottom;
      }
      // no round end: a box (a cell inside the plot)
      if ((round & 3) === 0) g[4] = 0;
      g[5] = (round & 3) === 0 ? 0 : round;
    }
  }
  g[2] = box.width;
}
