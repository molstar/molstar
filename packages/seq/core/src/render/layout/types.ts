/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The section layout interface (_spec/README.md "Render", Section layouts): a layout cuts a section's display units
// into strips (each one `AxisTransform` over a stack of bands, a track a band) and places each on the page with a warp;
// layer types lay out and draw a strip at a time in its own frame, hits come back through the placement's inverse, and
// the viewer's navigation is the layout's arithmetic (`SectionNav`). Internal: no application implements one yet.

import type { Issue } from '@molstar/core/repr/tree';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedSection, ResolvedTrack } from '@molstar/seqx-core/view';
import type { AxisTransform, Viewport } from '../axis.js';
import type { RenderConfig } from '../config.js';
import type { PassName, PassWriter } from '../draw-lists.js';
import type { Box, Emphasis, Hit, LabelSink } from '../layer-type.js';
import type { TextMetrics } from '../text.js';
import type { SeqTokens } from '../tokens.js';
import type { SectionLayout } from '../tracks.js';
import type { Strips } from './strips.js';
import type { PolarWarp } from './warp.js';

/**
 * What `Placement.at` finds under a page point: a band (`data`), a ring's outside label, a text's gutter (row labels
 * and numbers), a window, a gap between groups or lines (no position for a hover; the nearest for a drag), nothing.
 */
export type PointKind = 'data' | 'label' | 'gutter' | 'window' | 'gap' | 'empty';

/** A window's parts: its body (dragged: pans) and its two ends (dragged: resize about the other end). */
export type WindowPart = 'body' | 'start' | 'end';

/**
 * A point placed: its strip and band, the point, its kind, the display unit there. Reused. From `at` the point is in
 * the strip's frame (what a layer's hit reads); from `place` it is on the page (where the unit is drawn: the middle of
 * the strip's first band with height), so `at(place(o))` gives `o` back.
 */
export class PlacePoint {
  strip = -1;
  band = -1;
  /** The point: in the strip's frame from `at`, on the page from `place`. */
  x = 0;
  y = 0;
  kind: PointKind = 'empty';
  /** The display unit (fractional) there, −1 for none. */
  offset = -1;
}

/** A section laid out: its box, its whole content's height, its strips, and the inverse of its warps. */
export interface Placement {
  /** The section's region on the page (its visible box: a text section scrolls inside it). */
  readonly box: Box;
  /** All of its content (text: every line), for scrolling and export. */
  readonly contentHeight: number;
  readonly strips: Strips;
  /** What is under a page point, into a reused record; false when the point is outside the section. */
  at(x: number, y: number, out: PlacePoint): boolean;
  /** The display unit under a point, the nearest when between strips or in a gutter (drags, zoom anchors). */
  offsetAt(x: number, y: number): number;
  /**
   * Where a display unit is drawn: its strip, its first band with height and the page point there (the cursor, reveal,
   * section-scope layers); `at` of that point gives the unit back.
   */
  place(offset: number, out: PlacePoint): boolean;
  /** A polar strip's warp (`strips.warp` polar). */
  warpOf?(strip: number): PolarWarp;
  /** Whether its strips are clipped to its box (a section that scrolls inside it). */
  readonly clip?: boolean;
  /** What a follower of the section follows: the display units it shows (text: its lines in view). */
  readonly window?: Viewport;
}

/** What a layout reads of the view and the section: never a layer's params. */
export interface PlaceContext {
  readonly width: number;
  readonly height: number;
  readonly dpr: number;
  readonly config: RenderConfig;
  readonly tokens: SeqTokens;
  readonly text: TextMetrics;
  /** The view’s model (its systems’ numberings: a text’s numbers; a ring’s title: its system’s name and alphabet). */
  readonly model?: StaticModel;
  /** The section's region: its x and width (the data box's, past the label column), its top, the page below it. */
  readonly box: Box;
  /** The section's transform over its window (its pieces, its direction, the viewport or a leader's image). */
  readonly axis: AxisTransform;
  /** The section's bands: its tracks as shown (templates expanded), with their heights. */
  readonly bands: readonly ResolvedTrack[];
  readonly heights: readonly number[];
  /**
   * Lays out the layers of a strip's bands (each in its box in the strip's frame) and records what each band used
   * (`strips.used`): what a layout that collapses empty bands reads before it fixes the strip's height.
   */
  layStrip(placement: Placement, strip: number, labels?: LabelSink): void;
  /** Says what the layout could not place (a text past `text.maxLines`). */
  report(issue: Issue): void;
  /** Lays out all of the section's content (an export: a text's every line), not only what is in view. */
  readonly figure?: boolean;
  /** Draws the windows an export leaves out by default (`export.window`: a ring's). */
  readonly windows?: boolean;
  /**
   * A grid section's y (_spec/README.md "Render", The grid layout): its axis, the pieces it shows and its window in
   * display units (x's, a leader's image, its own `viewport.y`, or all of it); undefined while its length is unknown.
   */
  readonly y?: {
    readonly axis: {
      readonly system: string;
      readonly dir: 1 | -1;
      readonly length?: number;
      readonly cyclic: boolean;
    };
    readonly pieces: readonly number[];
    readonly viewport?: Viewport;
  };
}

/** The pure arithmetic of a section's viewport (_spec/README.md "Render", Section layouts, Navigation). */
export interface NavContext {
  /** Display units in the section. */
  readonly total: number;
  /** The data box's width and visible height, px. */
  readonly width: number;
  readonly height: number;
  readonly cyclic: boolean;
  /** The most px a position (the viewer's `viewport.maxPxPerPosition`), and the least (text: the letters' threshold). */
  readonly maxScale: number;
  readonly minScale?: number;
  /** The section's placement, when the layout needs it (text: lines). */
  readonly placement?: Placement;
  /** A grid's y: its display units, when y has a window of its own (`viewport.y`). */
  readonly totalY?: number;
}

/** A cursor's step: along the reading order, a line up or down, to the ends. */
export type StepDir = 'left' | 'right' | 'up' | 'down' | 'home' | 'end';

/** Gestures to viewports, pure: what the viewer calls instead of doing arithmetic of its own. */
export interface SectionNav {
  /** Moves the view by a scroll's deltas, px: positive shows what is right of and below what is in view. */
  pan(c: NavContext, v: Viewport, dx: number, dy: number): Viewport;
  /** Zooms by `factor` (above 1: in) about display unit `at`, kept where it is. */
  zoom(c: NavContext, v: Viewport, factor: number, at: number): Viewport;
  fit(c: NavContext, v: Viewport): Viewport;
  /**
   * The least move that shows display units `[start, end)` (the cursor, a linked selection), unchanged when it does: a
   * window pans to the nearer edge, a text scrolls the line into view.
   */
  reveal(c: NavContext, v: Viewport, start: number, end: number): Viewport;
  /**
   * The viewport a focus on `[start, end)` gives (a double-click, an address): a window takes the range; a text brings
   * its first line to the top, its letter size kept. Absent: the range itself, clamped.
   */
  focus?(c: NavContext, v: Viewport, start: number, end: number): Viewport;
  /** The display unit a cursor at `offset` steps to, or −1 for none. */
  step(c: NavContext, offset: number, dir: StepDir): number;
  clamp(c: NavContext, v: Viewport): Viewport;
}

/** A section's layout: tracks, text, rings. `Kept` is what it keeps between frames (its buffers). */
export interface SectionLayoutType<Kept = unknown> {
  readonly name: string;
  /** Whether bands are filled with the track ink in the background pass (tracks: a track a box). */
  readonly fillBands: boolean;
  /** A track's height by its layers' types when the spec gives none (the config's `heights` when absent). */
  heights?(config: RenderConfig): RenderConfig['heights'];
  /** Cuts the section into strips for its box and viewport, into prev's buffers. */
  place(section: ResolvedSection, ctx: PlaceContext, viewport: Viewport | undefined, prev?: Kept): Placement & Kept;
  /** Whether two viewports place alike: the part of the section's key this layout reads. */
  sameKey(a: Viewport | undefined, b: Viewport | undefined): boolean;
  readonly nav: SectionNav;
  /**
   * Its window may pass a cyclic axis's origin and draws there (tracks: the system twice; a ring's): the viewer's
   * navigation then wraps it (`NavContext.cyclic`). Absent, a window stays within its axis (text: an anchor).
   */
  readonly wraps?: boolean;
  /**
   * What a follower of a section in this layout follows when it is not the transform's window (a ring: its window,
   * past the total when it wraps the origin; text: its lines in view). The layout then places over the whole axis.
   */
  followWindow?(
    section: ResolvedSection,
    axis: AxisTransform,
    viewport: Viewport | undefined,
    ctx: {
      readonly width: number;
      readonly config: RenderConfig;
      readonly text: TextMetrics;
      readonly heights: readonly number[];
    },
  ): Viewport;
  /**
   * Sets the window a placement draws (a ring's: what its follower shows), its viewport's or its leader's range, on a
   * kept placement too: the window is out of the section's key (_spec/README.md "Render", Rings layout). `active`:
   * something follows it (or it follows); `shown`: drawn at all (`rings.window`, an export's `windows`).
   */
  window?(placement: Placement, viewport: Viewport | undefined, active: boolean, shown: boolean): void;
  /**
   * The layout's own marks of a pass (label column, gutters, numbers, a ring's guides and labels, a window): once a
   * section with `band` −1 before its strips, then once a band before its layers (`strip`, `band`), in the strip's frame.
   */
  chrome(
    section: SectionLayout,
    pass: PassName,
    out: PassWriter,
    emphasis: Emphasis | undefined,
    strip: number,
    band: number,
  ): void;
  /** What its own marks hit (a filter chip, a gutter's row, a window), before the bands' layers. */
  hitChrome?(section: SectionLayout, x: number, y: number): Hit | undefined;
  /**
   * Its strips carry a second transform (`strip.y`, a grid's): every layer is laid out a strip at a time with it
   * (`ctx.axisY`), a `section` layer too (pairs as dots), and its bands overlay the strip.
   */
  readonly axesY?: boolean;
}
