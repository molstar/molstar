/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What a place on data read by region is, in words (_spec/README.md "Hover card", By region): a pairs bin, a region
// grid's cell, a matrix cell, each by the core's hover texts (`@molstar/seqx-core/render`: `pairsBinHover`,
// `gridCellHover`, `matrixCellHover`); the card lays them out as its rows, the announcer says their sentence.

import type { Key } from '@molstar/seqx-core/query';
import {
  type GridLay,
  gridCellHover,
  type LayerLayout,
  type MatrixLay,
  matrixCellHover,
  type PairsLay,
  pairsBinHover,
  type RegionHover,
  type SystemName,
} from '@molstar/seqx-core/render';

/**
 * The words of a place on a layer read by region: its element (a pairs bin's key) or its cell (`i` on x, `j` on its y
 * system: a grid's second axis, a matrix's state; NaN when none); undefined for anything else.
 */
export function regionHoverOf(
  l: LayerLayout | undefined,
  element: Key | undefined,
  i: number,
  j: number,
  name: SystemName,
): RegionHover | undefined {
  if (!l) return undefined;
  switch (l.type.name) {
    case 'pairs':
      return element ? pairsBinHover(l.lay as PairsLay, element.key, name) : undefined;
    case 'grid':
      return Number.isNaN(i) || Number.isNaN(j) ? undefined : gridCellHover(l.lay as GridLay, i, j, name);
    case 'matrix':
      return Number.isNaN(i) || Number.isNaN(j) ? undefined : matrixCellHover(l.lay as MatrixLay, i, j, name);
    default:
      return undefined;
  }
}
