/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Issue } from '@molstar/core/repr/tree';
import type { Want } from '@molstar/seqx-core/data';
import type { OverlayNode, StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer, ResolvedSection, ResolvedTrack, ResolvedView } from '@molstar/seqx-core/view';
import { AxisTransform, type Viewport } from './axis.js';
import type { RenderConfig } from './config.js';
import type { Box, LayerType, Overlay, StoreReader } from './layer-type.js';
import type { ThemeContext } from './channels.js';
import type { TextMetrics } from './text.js';
import type { SeqTokens } from './tokens.js';
import { LabelLane, type LayerContent, type TrackSizing, type TrackSizingHook } from './sizing.js';
import type { Placement, SectionLayoutType } from './layout/types.js';
import { Rail } from './rail.js';
import type { ImageCache } from './images.js';
import type { SectionCache } from './tracks-cache.js';

// The tracks layout's records (_spec/README.md "Render", Tracks layout): what it is given (`LayoutInput`), what it
// gives (the view's, sections', tracks' and layers' layouts) and what it keeps between frames (`LayoutState`).

export interface LayerLayout {
  readonly key: string;
  readonly layer: ResolvedLayer;
  readonly type: LayerType<any, any>;
  readonly lay: unknown;
}

export interface TrackLayout {
  readonly track: ResolvedTrack;
  /** The data box: its band, inside its well by its padding (what its layers lay out, hit and emphasize in). */
  readonly box: Box;
  /** Its well: the data box grown by its padding (_spec/README.md "Render", Tracks layout), filled in the background. */
  readonly well: Box;
  /**
   * The label as drawn (cut to the label column with an ellipsis, never a hard clip).
   * @internal
   */
  readonly labelText: string;
  /** The label's measured width and its baseline (drawing never measures). */
  readonly labelWidth: number;
  readonly labelY: number;
  readonly layers: readonly LayerLayout[];
  /**
   * The "n filtered" chip: the rows a layer's filter left out (the first layer with some), its text and its box at the
   * label column's right end; a click on it clears the filter (`target: 'filter'`).
   */
  readonly filter?: {
    readonly layer: string;
    readonly count: number;
    readonly text: string;
    readonly x: number;
    readonly width: number;
    /** Its line below the label's (px), on a track tall enough for two lines: the label then whole; else 0, beside it. */
    readonly dy: number;
  };
}

export interface SectionLayout {
  readonly section: ResolvedSection;
  /** The section's id (its address's last step): what a hit names it by. */
  readonly id: string;
  /** Its layout (_spec/README.md "Render", Section layouts) and where it placed the section. */
  readonly layout: SectionLayoutType<any>;
  readonly placement: Placement;
  /** Undefined while the section's length is unknown (pending); else its transform over its window (a follower reads it). */
  readonly axis?: AxisTransform;
  /** Its tracks on the page and the gaps between them, the data's x: the column a `column` emphasis spans. */
  readonly band?: Box;
  readonly y: number;
  readonly height: number;
  /** Its tracks as bands: on the page in tracks; in text, those of its first line. */
  readonly tracks: readonly TrackLayout[];
  /** The breaks in view, `[x, hiddenStart, hiddenEnd, …]` (x unsnapped): a hairline each across the band. */
  readonly breaks: readonly number[];
  /** Its ruler rail, when it draws one (_spec/README.md "Render", Navigators). */
  readonly rail?: Rail;
  /** Its rows layers' rails (upright, at the right edge of each band), when they draw them. */
  readonly rowsRails?: readonly Rail[];
  /** The direction mark in the first track's label row (`axis.directionMark`): its text, x and measured width. */
  readonly mark?: { readonly text: string; readonly x: number; readonly width: number };
  /** The layers' lays a strip and band (`strip * bands + band`), in the strip's frame. */
  readonly lays: readonly (readonly LayerLayout[])[];
  /** Each strip's breaks in view (its frame), as `breaks`. */
  readonly stripBreaks: readonly (readonly number[])[];
  /** Its `section` layers' lays: laid out once a section, drawn in page coordinates. */
  readonly sectionLays: readonly LayerLayout[];
  /** Each band's box a strip (`strip * bands + band`, its frame) and each strip's band column (none when empty). */
  readonly boxes: readonly Box[];
  readonly columns: readonly (Box | undefined)[];
  /** What its chrome reads of the view. */
  readonly pageHeight: number;
  readonly tokens: SeqTokens;
  readonly config: RenderConfig;
}

export interface ViewLayout {
  readonly width: number;
  readonly height: number;
  readonly dpr: number;
  readonly tokens: SeqTokens;
  readonly config: RenderConfig;
  readonly sections: readonly SectionLayout[];
  /** What the layers need from the store now. */
  readonly wants: readonly Want[];
  /** The issues of this layout (a color it cannot read, a layer it cannot draw), each once. */
  readonly issues: readonly Issue[];
  /** The dynamic nodes this layout shows, with what they name (with `nodes.holdShown`): its owner holds them. */
  readonly shown: readonly OverlayNode[];
}

export interface LayoutInput {
  readonly model: StaticModel;
  readonly view: ResolvedView;
  readonly size: { readonly width: number; readonly height: number };
  /** The window of display units of every section that does not follow another; all of each by default. */
  readonly viewport?: Viewport;
  /** A section's own window (by its address), before `viewport`: a caller placing sections one by one. */
  readonly viewports?: Readonly<Record<string, Viewport>>;
  readonly dpr?: number;
  readonly tokens: SeqTokens;
  readonly text: TextMetrics;
  readonly config?: RenderConfig;
  readonly data: StoreReader;
  readonly overlay?: Overlay;
  readonly settled?: boolean;
  /** The layer types drawn, by name (the composition's registry): none when absent (every layer is said, not drawn). */
  readonly types?: Readonly<Record<string, LayerType<any, any>>>;
  /** The schemes and measures themes read (the built-ins by default) and the reader's color-blind substitution. */
  readonly theme?: Partial<Pick<ThemeContext, 'schemes' | 'measures' | 'cvdSafe'>>;
  /** The section layouts by name (the composition's registry): tracks for every section when absent. */
  readonly layouts?: Readonly<Record<string, SectionLayoutType<any>>>;
  /** Lays out every line of a text section, its box its content (an export: the figure). */
  readonly figure?: boolean;
  /** Draws the windows a figure leaves out by default (a ring's in an export). */
  readonly windows?: boolean;
  /** The clades a rows layer folds by its track's address, where its band gives none (an export's `folds`). */
  readonly folds?: Readonly<Record<string, readonly string[]>>;
  /**
   * The track sizings a section's `heights` may name besides the built-ins (`content`, `uniform`, `fixed`): the
   * composition's `track-sizing` registry entries, by name (_spec/README.md "Render", Track heights).
   */
  readonly sizings?: Readonly<Record<string, TrackSizing>>;
  /** The host's `trackSizing` port: the last word on each section's tracks' heights. */
  readonly trackSizing?: TrackSizingHook;
  /**
   * A navigator's picture: its tracks small, no label lane (its labels are the view's).
   * @internal
   */
  readonly picture?: boolean;
  /** The shared image cache (_spec/README.md "Render", The shared reduction): the renderer's. */
  readonly images?: ImageCache;
  /**
   * The sections that draw a ruler rail, and the rows tracks that draw a rows rail, by address (_spec/README.md
   * "Render", Navigators): part of each section's key by identity (a slim rail takes a row, a rows rail a column at the
   * right edge). Which can move is set on the rails each frame by the renderer, out of the key: idle wakes with no
   * relayout.
   */
  readonly rails?: Readonly<Record<string, unknown>>;
  /**
   * Another layout's held domains (_spec/README.md "Render", Navigators: a navigator's picture colors by its view's,
   * never by a domain of its own over the whole axis): looked up by `map(key)`; this layout then holds none.
   */
  readonly domains?: HeldDomains;
}

/** A layout's held domains as another reads them: its state, and a key here to a key there. */
export interface HeldDomains {
  readonly state: LayoutState;
  map(key: string): string;
}

/** What a layout keeps between frames: each layer's prepared form (per resolved layer) and its last layout. */
export class LayoutState {
  readonly preps = new WeakMap<ResolvedLayer, unknown>();
  readonly lays = new Map<string, { type: LayerType<any, any>; lay: unknown }>();
  /**
   * Open scales' domains (_spec/README.md "Themes", Palettes and scales): the union over the layers naming each, of their
   * values in view, computed when the view settles and held while it moves.
   */
  domains = new Map<string, readonly [number, number]>();
  /** Bumped when `domains` changes: part of a section's reuse. */
  domainsVersion = 0;
  /** Each section's last layout and what it was laid out for (_spec/README.md "Render", The key: per section). */
  readonly sections = new Map<string, SectionCache>();
  /** The row columns' fit, kept while what it read is the same. */
  columnsFit: unknown;
  /** The template instances shown last (their tracks' addresses): kept while a zoom's new level is pending. */
  readonly instancesShown = new Set<string>();
  /**
   * Track labels cut to the label column, by label and width (the metrics and font are the state's renderer's), for
   * the metrics' version (`elidedVersion`: a web font arriving cuts them again).
   */
  readonly elided = new Map<string, Map<number, string>>();
  elidedVersion = -1;
  /** Each section's tracks' heights and contents while their key is (`bandHeights`), and as a leader's reading. */
  readonly sizes = new Map<string, SizedHeights>();
  readonly leadSizes = new Map<string, SizedHeights>();
  /** The bands' label lanes, by band index (one a band, reset before its layers are laid out). */
  readonly labelLanes: LabelLane[] = [];
}

/** A section's tracks' heights with what they were sized for (`layoutView`'s `bandHeights`). */
export interface SizedHeights {
  readonly tracks: readonly ResolvedTrack[];
  readonly fit: number;
  readonly config: RenderConfig;
  readonly layout: SectionLayoutType<any>;
  readonly version: number;
  readonly sizings: LayoutInput['sizings'];
  readonly hook: LayoutInput['trackSizing'];
  readonly heights: readonly number[];
  /** Each track's content (one lane block a track), undefined for a track not sized by its content. */
  readonly contents: readonly (LayerContent | undefined)[];
  readonly issue?: Issue;
}

/** A layer's key: its address, else its track's and its index. */
export function layerKey(track: ResolvedTrack, layer: ResolvedLayer): string {
  return layer.address ?? `${track.address}:${layer.index}`;
}

/** A layer's type: the first of its fallbacks drawn here (a loop: no array made a layer a frame). */
export function firstType(
  layer: ResolvedLayer,
  types: Readonly<Record<string, LayerType<any, any>>>,
): LayerType<any, any> | undefined {
  const f = layer.fallbacks;
  for (let i = 0; i < f.length; i++) {
    const t = types[f[i]];
    if (t !== undefined) return t;
  }
  return undefined;
}
