/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The built-in presets (seqspec:viewer#extension-points): a section's geometry, render props over the context's, its
// navigator and per-behavior overrides (`forSection`). Picked per section by the layout (`ui.style`, the section's
// `style` through `Presets.Alias`, else `Presets.Default` by layout and axis).

import type { RenderProps } from "./layout.ts";
import type { NavigatorParts } from "./navigators.ts";

export interface PresetEntry {
  readonly name: string;
  readonly geometry: "linear" | "wrapped" | "circular";
  /** Over `ctx.renderProps`, for its sections; the vertical overview is `render.verticalOverview`. */
  readonly render?: Partial<RenderProps>;
  /** A registered navigator's name, or parts (`navigatorFrom`); none when left out. */
  readonly navigator?: string | Partial<NavigatorParts>;
  /** Per behavior: off (`false`), or params over its own (seqspec:viewer#behaviors, `forSection`). */
  readonly behaviors?: Readonly<Record<string, false | object>>;
}

/**
 * Text (seqspec:interaction#selection): selected as text is, a plain drag on the letters selecting positions line by
 * line, with the crosshair over them; the rows scroll by the wheel, the scrollbar, the keyboard, a selection drag past
 * their edges, and a drag on the text's inset and its right strip (`TEXT_STRIP`: a touch reader's scroll), never across;
 * no zoom.
 */
const textBehaviors = {
  // A double click in text is a click (nothing takes it).
  select: { dragRegions: ["data", "column"], repeatClicks: "click" },
  // Its default regions: on the data `select` claims a drag first (priority), so the inset and the strip pan.
  dragPan: { x: false },
  cursor: {
    byRegion: {
      data: "crosshair",
      column: "crosshair",
      gutter: "grab",
      header: "pointer",
      navigator: "pointer",
      vnav: "pointer",
    },
  },
  wheel: { zoom: [] },
} as const;
/** Px right of the text that a drag scrolls (the text's inset is the least): a touch reader's grip on the rows. */
const TEXT_STRIP = 24;
export const DefaultPresets: readonly PresetEntry[] = [
  // No row borders (seqspec:viewer#render-props): a document opts in (`style.row_borders`), a host by its render props.
  { name: "tracks", geometry: "linear", navigator: "rulerRail" },
  {
    name: "msa",
    geometry: "linear",
    // A header with room for a tree and the labels; gutters for each row's first and last residue numbers.
    render: {
      headerWidth: 200,
      rowHeight: 18,
      gutterLeft: 40,
      gutterRight: 40,
      verticalOverview: 72,
      verticalNav: "both",
    },
    navigator: "rulerRail",
  },
  { name: "strip", geometry: "linear", render: { headerWidth: 0 }, navigator: "strip", behaviors: { select: false } },
  {
    name: "wrapped",
    geometry: "wrapped",
    render: { headerWidth: 0, gutterRight: TEXT_STRIP },
    navigator: "none",
    behaviors: textBehaviors,
  },
  {
    name: "compact",
    geometry: "wrapped",
    render: {
      headerWidth: 0,
      gutterRight: TEXT_STRIP,
      wrappedGutter: 8,
      wrappedGroup: 10,
      wrappedGroupGap: 4,
      wrappedGap: 0,
      wrappedTrackGap: 0,
      wrappedLaneHeight: 4,
    },
    navigator: "none",
    behaviors: textBehaviors,
  },
  // Rings (seqspec:viewer#the-render-backend): the whole axis always, no header, gutters or navigator; no zoom or pan.
  {
    name: "circular",
    geometry: "circular",
    render: { headerWidth: 0, gutterLeft: 0, gutterRight: 0, rowBorders: false },
    navigator: "none",
    behaviors: { wheel: { zoom: [] }, dragPan: false },
  },
];
