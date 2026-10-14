/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Navigators (seqspec:viewer#extension-points): a navigator only draws; its gestures are the `navigator` behavior's
// (windows.ts).
// An entry is its parts; `navigatorFrom` makes one from a partial set, and the canvas engine draws any entry from its
// parts (canvas/navigators.ts) unless it brings its own draws. Presets name a navigator or give parts.

import type { NavView } from "./draw/nav.ts";
import type { RenderProps } from "./layout.ts";
import type { SeqViewerRegistry } from "./registry.ts";

export type NavigatorAxis = "overview" | "miniature" | "labels" | "hairline";
export interface NavigatorParts {
  /** The row, and the axis drawing centered in it. */
  readonly height: number;
  readonly thickness: number;
  /** Highlight channels' and the selection's extents as marks on the axis. */
  readonly highlightMark: boolean;
  /** What the axis shows, the first that applies. */
  readonly axis: readonly NavigatorAxis[];
  /** Track refs for the miniature; null: all but letters, ruler and rows. */
  readonly tracks: readonly string[] | null;
  /** A ground per block, its corner radius; a lane per block when stacked. */
  readonly ground: boolean;
  readonly radius: number;
  readonly lanes: boolean;
  /** The window's marks, combined: an outlined box, the outside faded by `fade`, a segment `segment` px thick. */
  readonly box: boolean;
  readonly fade: number;
  readonly segment: number;
  readonly ticks: boolean;
  /** Stronger marks under the pointer; the hovered position's line. */
  readonly hoverGrow: boolean;
  readonly hoverMark: boolean;
  /** Px near a window end that grab it to zoom. */
  readonly edge: number;
  /** No row of its own: drawn on the pinned ruler's bottom border (`rulerRail`); without a ruler, `rail`. */
  readonly onRuler: boolean;
}
export const DefaultNavigatorParts: NavigatorParts = {
  height: 20,
  thickness: 12,
  highlightMark: true,
  axis: ["overview", "miniature"],
  tracks: null,
  ground: true,
  radius: 2,
  lanes: false,
  box: false,
  fade: 0.62,
  segment: 0,
  ticks: true,
  hoverGrow: true,
  hoverMark: true,
  edge: 5,
  onRuler: false,
};

export interface NavigatorEntry {
  readonly name: string;
  readonly parts: NavigatorParts;
  /** 0: no row. */
  height(props: RenderProps): number;
  readonly edge: number;
  /** The whole axis, and the view's window on the overlay. */
  draw?(n: NavView): void;
  drawWindow?(n: NavView): void;
}

export function navigatorFrom(parts: Partial<NavigatorParts>, name = "custom"): NavigatorEntry {
  const p: NavigatorParts = { ...DefaultNavigatorParts, ...parts };
  return { name, parts: p, edge: p.edge, height: () => (p.onRuler ? 0 : p.height) };
}

const rail = { thickness: 2, axis: ["hairline"], ground: false, fade: 0, segment: 3, ticks: false } as const;
export const DefaultNavigators: readonly NavigatorEntry[] = [
  navigatorFrom({ ...rail, height: 0, onRuler: true }, "rulerRail"),
  navigatorFrom({ height: 36, thickness: 28, lanes: true, box: true, fade: 0 }, "band"),
  navigatorFrom({ height: 20, thickness: 12 }, "strip"),
  navigatorFrom({ ...rail, height: 12 }, "rail"),
  navigatorFrom({ height: 0, axis: [], highlightMark: false, hoverMark: false }, "none"),
];

const made = new WeakMap<object, NavigatorEntry>();
/** A preset's navigator: a registered name (unknown: none), or parts over the defaults (one entry per parts). */
export function navigatorOf(
  registry: SeqViewerRegistry,
  nav: string | Partial<NavigatorParts> | undefined,
): NavigatorEntry | null {
  if (nav === undefined) return null;
  if (typeof nav === "string") return registry.get("navigators", nav) ?? null;
  let e = made.get(nav);
  if (!e) {
    e = navigatorFrom(nav);
    made.set(nav, e);
  }
  return e;
}
