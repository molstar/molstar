/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { WindowPart } from '@molstar/seqx-core/render';

// The cursor over a window's part (_spec/README.md "Composition", Navigators): the canvas's rails and the navigators
// read it, so it is here, not with the navigators (a viewer that lists none draws rails).

export type PartCursor =
  | 'ew-resize'
  | 'ns-resize'
  | 'nwse-resize'
  | 'move'
  | 'grab'
  | 'grabbing'
  | 'pointer'
  | 'default';

/**
 * The cursor that says the gesture over a window's part: a grip resizes along the axis it moves (`ew-resize` along
 * positions, `ns-resize` along rows), the body grabs (`grabbing` while dragged), outside it a click centres (`pointer`;
 * a press outside keeps the default while it runs).
 */
export function partCursor(part: WindowPart, axis: 'positions' | 'rows', dragging: boolean): PartCursor {
  if (part === 'outside') return dragging ? 'default' : 'pointer';
  if (part === 'body') return dragging ? 'grabbing' : 'grab';
  return axis === 'rows' ? 'ns-resize' : 'ew-resize';
}
