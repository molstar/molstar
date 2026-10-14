/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Element descriptions (seqspec:libraries#nightingale): plain data that `mount` turns into Nightingale's elements,
// compared by tests without a DOM; and the supported set to resolve with.

import { KernelRegistry } from "@molstar/seq-core/format";

export interface ElementDesc {
  tag: string;
  attributes: Record<string, string>;
  /**
   * Data arrays and text, set as properties; `svs` is mount's metadata, never set on the element: a manager's axis
   * (`AxisMeta`), a linegraph's progressive data to fetch (`StreamedValues`).
   */
  properties: Record<string, unknown>;
  /** Imperative setters (`setHeatmapData`), called once the element is defined. */
  calls: { method: string; args: unknown[] }[];
  children: ElementDesc[];
}

export function el(
  tag: string,
  attributes: Record<string, string> = {},
  properties: Record<string, unknown> = {},
  children: ElementDesc[] = [],
  calls: ElementDesc["calls"] = [],
): ElementDesc {
  return { tag, attributes, properties, calls, children };
}

/** The kernel layer types Nightingale draws (`emphasis` is drawn as its fallback, `points`). */
const NightingaleTypes = [
  "letters",
  "rows",
  "ranges",
  "points",
  "variants",
  "secondary",
  "values",
  "color",
  "ruler",
  "pairs",
  "logo",
  "matrix",
] as const;

/**
 * The supported set to resolve with, the same for both layouts: its layer types and the kernel's themes; no `char` or
 * `path` glyphs and no extension types or themes, so those resolve to their fallbacks.
 */
export const NightingaleSupported: ReadonlySet<string> = new Set([
  ...NightingaleTypes,
  ...(KernelRegistry.themes ?? []).map((t) => t.name),
]);
