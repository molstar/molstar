/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer's configuration items of the first slice (_spec/README.md "Configuration"): one default object of the
// shape the models read. Where a behavior could go more than one way, it is an item here with the recommended default.

/** What the models and behaviors read; read-only after composition. */
export interface ViewerConfig {
  readonly foundations: { readonly preset: 'comfortable' | 'compact' };
  readonly density: {
    /** The least pixels a position needs for its letter to be drawn. */
    readonly letterMinPx: number;
    /** The least pixels a position needs for a dot. */
    readonly dotMinPx: number;
  };
  readonly bindings: { readonly preset: 'default' | 'genome-browser' };
  readonly store: {
    /** The loaded bytes the data store keeps (64 MB). */
    readonly cacheBytes: number;
    /** Requests in flight at most (6). */
    readonly inFlight: number;
    /** Bins per tile (1,024). */
    readonly tileBins: number;
  };
  readonly hover: {
    /** Project the hover onto every linked section (true); false: only sections of the hovered system. */
    readonly projectLinked: boolean;
    /** Close the hover before projecting: a nucleotide lights its whole codon and its residue (true, EC-30, EC-144). */
    readonly closeCodons: boolean;
  };
  readonly viewport: {
    /** The most pixels a position may take when zoomed in: the zoom's upper bound. */
    readonly maxPxPerPosition: number;
    /** Sections on equal axes share one viewport (true, EC-113); false: a viewport per section. */
    readonly groupEqualAxes: boolean;
    /**
     * When a section's width changes: keep its zoom about the center (`keep-zoom`, the default: a section hidden for
     * lack of width comes back as it was, EC-112) or keep its range (`keep-range`). A fitted viewport stays fitted.
     */
    readonly onResize: 'keep-zoom' | 'keep-range';
  };
  readonly follow: {
    /**
     * A move asked of a follower (followers never lead, EC-114): `ignore` it (the default) or `detach` the follower,
     * which then has its own viewport until it is attached again (`attach`, or a fit).
     */
    readonly followerMoves: 'ignore' | 'detach';
    /**
     * A follower's range over a cyclic origin: kept as the wrap, two pieces (`wrap`, the default), or the hull of
     * the image (`hull`).
     */
    readonly cyclicArc: 'wrap' | 'hull';
  };
  readonly navigate: {
    /** The zoom factor of one wheel notch (100 px of delta) with the zoom modifier, and of a key zoom. */
    readonly wheelZoomStep: number;
    /** The part of the viewport an arrow key pans. */
    readonly keyPanFraction: number;
    /** Wheel deltas within one frame are applied once, at the frame (true); false: at each event. */
    readonly coalesceWheel: boolean;
  };
}

export const DefaultViewerConfig: ViewerConfig = {
  foundations: { preset: 'comfortable' },
  density: { letterMinPx: 7, dotMinPx: 2 },
  bindings: { preset: 'default' },
  store: { cacheBytes: 64 * 1024 * 1024, inFlight: 6, tileBins: 1024 },
  hover: { projectLinked: true, closeCodons: true },
  viewport: { maxPxPerPosition: 32, groupEqualAxes: true, onResize: 'keep-zoom' },
  follow: { followerMoves: 'ignore', cyclicArc: 'wrap' },
  navigate: { wheelZoomStep: 1.25, keyPanFraction: 0.1, coalesceWheel: true },
};

type DeepPartial<T> = { readonly [K in keyof T]?: T[K] extends object ? DeepPartial<T[K]> : T[K] };

/** The defaults with the given items over them, one level of groups deep. */
export function viewerConfig(over: DeepPartial<ViewerConfig> = {}): ViewerConfig {
  const out: Record<string, unknown> = {};
  for (const [group, items] of Object.entries(DefaultViewerConfig))
    out[group] = { ...items, ...(over as Record<string, object | undefined>)[group] };
  return out as unknown as ViewerConfig;
}
