/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Issue } from '@molstar/core/repr/tree';
import type { Coverage, Want } from '@molstar/seqx-core/data';
import type { ResolvedSection } from '@molstar/seqx-core/view';
import { AxisTransform, type Viewport } from './axis.js';
import type { RenderConfig } from './config.js';
import type { Overlay, StoreReader } from './layer-type.js';
import type { ThemeContext } from './channels.js';
import type { PlaceContext, SectionLayoutType } from './layout/types.js';
import type { LayoutInput, LayoutState, SectionLayout } from './tracks-types.js';

// A section's reuse (_spec/README.md "Render", The key: per section): its reads of the store, what its last layout gave
// the frame, and whether the view's parts, its place, axis and window are the same so its kept layout is drawn again.

/**
 * A section's reads of the store (_spec/README.md "Render", The key: per section): the entries and coverages its
 * layout read with their versions, so a section whose reads hold is not laid out again. Its arrays are kept from layout
 * to layout (a count, not a new array: a layout allocates a bounded amount).
 */
/** A store that bounds no layer: every load fits. */
const Fits = { fits: true } as const;

class SectionReads implements StoreReader {
  private keys: string[] = [];
  private versions: number[] = [];
  private n = 0;
  private sources: string[] = [];
  private coverages: (Coverage | undefined)[] = [];
  private m = 0;
  private refusalsAt: number | undefined;
  private data!: StoreReader;
  begin(data: StoreReader): void {
    this.data = data;
    this.n = 0;
    this.m = 0;
    this.refusalsAt = data.refusals;
  }
  tilesOf(w: Want) {
    return this.data.tilesOf?.(w) ?? [];
  }
  /** A layer's refusal: a section is laid out again when any changes (`refusals`). */
  layerRefused(layer: string) {
    return this.data.layerRefused?.(layer);
  }
  layerLoad(wants: readonly Want[]) {
    return this.data.layerLoad?.(wants) ?? Fits;
  }
  get refusals() {
    return this.data.refusals;
  }
  coverage(source: string) {
    const c = this.data.coverage(source);
    this.sources[this.m] = source;
    this.coverages[this.m++] = c;
    return c;
  }
  keysOf(w: Want) {
    return this.data.keysOf(w);
  }
  tooMany(w: Want) {
    return this.data.tooMany?.(w) ?? false;
  }
  get settings() {
    return this.data.settings;
  }
  get(key: string) {
    const e = this.data.get(key);
    this.keys[this.n] = key;
    this.versions[this.n++] = e?.version ?? -1;
    return e;
  }
  /** Whether they read the same again (read through `data`, so a renderer's reader records them as read). */
  hold(data: StoreReader): boolean {
    for (let i = 0; i < this.n; i++) if ((data.get(this.keys[i])?.version ?? -1) !== this.versions[i]) return false;
    for (let i = 0; i < this.m; i++) if (data.coverage(this.sources[i]) !== this.coverages[i]) return false;
    if (data.refusals !== this.refusalsAt) return false;
    return true;
  }
}

/** Where a section's layout reports to this frame: the view's wants, issues, held domains, lays used, nodes shown. */
export interface FrameSinks {
  readonly wants: Want[];
  readonly report: (i: Issue) => void;
  readonly theme: ThemeContext;
  readonly used: Set<string>;
  readonly showing: Set<string>;
}

/**
 * What a section kept from its last layout, what it was laid out for, and what that layout gave the frame (wants,
 * issues, extents, lays used, nodes shown), replayed when it is reused. One a section, kept across frames: its arrays,
 * its reads, its theme and its sinks are made once (a layout allocates a bounded amount).
 */
export class SectionCache {
  layout!: SectionLayout;
  section!: ResolvedSection;
  stamp!: ViewStamp;
  y = 0;
  axis: AxisTransform | undefined;
  gridY: PlaceContext['y'];
  viewport: Viewport | undefined;
  readonly reads = new SectionReads();
  private wants: Want[] = [];
  private nw = 0;
  private issues: Issue[] = [];
  private ni = 0;
  private used: string[] = [];
  private nu = 0;
  private showing: string[] = [];
  private ns = 0;
  private extentKeys: string[] = [];
  private extentValues: number[] = [];
  private ne = 0;
  private sinks!: FrameSinks;

  /** Starts a layout of the section for a frame. */
  begin(data: StoreReader, sinks: FrameSinks): void {
    this.reads.begin(data);
    this.sinks = sinks;
    this.nw = 0;
    this.ni = 0;
    this.nu = 0;
    this.ns = 0;
    this.ne = 0;
  }

  readonly want = (w: Want): void => {
    this.wants[this.nw++] = w;
    this.sinks.wants.push(w);
  };

  readonly report = (i: Issue): void => {
    this.issues[this.ni++] = i;
    this.sinks.report(i);
  };

  addUsed(key: string): void {
    this.used[this.nu++] = key;
    this.sinks.used.add(key);
  }

  addShowing(a: string, b: string): void {
    this.showing[this.ns++] = a;
    this.showing[this.ns++] = b;
    this.sinks.showing.add(a).add(b);
  }

  /** The themes a section's layers read: the frame's, with the extents it saw recorded. */
  readonly theme: ThemeContext = (() => {
    const self = this;
    return {
      get schemes() {
        return self.sinks.theme.schemes;
      },
      get measures() {
        return self.sinks.theme.measures;
      },
      get cvdSafe() {
        return self.sinks.theme.cvdSafe;
      },
      scale: (a) => self.sinks.theme.scale(a),
      domain: (k) => self.sinks.theme.domain(k),
      extent: (k, lo, hi) => {
        self.extentKeys[self.ne] = k;
        self.extentValues[2 * self.ne] = lo;
        self.extentValues[2 * self.ne + 1] = hi;
        self.ne++;
        self.sinks.theme.extent(k, lo, hi);
      },
    } as ThemeContext;
  })();

  /** Gives this frame what the kept layout gave its own. */
  replay(sinks: FrameSinks): void {
    for (let i = 0; i < this.nw; i++) sinks.wants.push(this.wants[i]);
    for (let i = 0; i < this.ni; i++) sinks.report(this.issues[i]);
    for (let i = 0; i < this.nu; i++) sinks.used.add(this.used[i]);
    for (let i = 0; i < this.ns; i++) sinks.showing.add(this.showing[i]);
    for (let i = 0; i < this.ne; i++)
      sinks.theme.extent(this.extentKeys[i], this.extentValues[2 * i], this.extentValues[2 * i + 1]);
  }
}

/** The view's parts a section's layout reads (size, ratio, tokens, text, overlay, config, types, themes, domains). */
interface ViewStamp {
  readonly model: string;
  readonly width: number;
  readonly height: number;
  readonly dpr: number;
  readonly tokens: number;
  readonly text: number;
  readonly overlay: number;
  readonly settled: boolean;
  readonly config: RenderConfig;
  readonly types: unknown;
  readonly schemes: unknown;
  readonly measures: unknown;
  readonly cvdSafe: boolean;
  readonly domains: number;
  readonly figure: boolean;
  readonly windows: boolean;
  readonly rails: unknown;
  readonly headers: unknown;
}

export function viewStamp(
  input: LayoutInput,
  config: RenderConfig,
  dpr: number,
  overlay: Overlay,
  state: LayoutState,
): ViewStamp {
  return {
    model: input.model.hash,
    width: input.size.width,
    height: input.size.height,
    dpr,
    tokens: input.tokens.version,
    text: input.text.version,
    overlay: overlay.version,
    settled: input.settled ?? false,
    config,
    types: input.types,
    schemes: input.theme?.schemes,
    measures: input.theme?.measures,
    cvdSafe: input.theme?.cvdSafe ?? false,
    domains: input.domains ? input.domains.state.domainsVersion : state.domainsVersion,
    figure: input.figure ?? false,
    windows: input.windows ?? false,
    rails: input.rails,
    headers: input.headers,
  };
}

function sameStamp(a: ViewStamp, b: ViewStamp): boolean {
  return (
    a.model === b.model &&
    a.width === b.width &&
    a.height === b.height &&
    a.dpr === b.dpr &&
    a.tokens === b.tokens &&
    a.text === b.text &&
    a.overlay === b.overlay &&
    a.settled === b.settled &&
    a.config === b.config &&
    a.types === b.types &&
    a.schemes === b.schemes &&
    a.measures === b.measures &&
    a.cvdSafe === b.cvdSafe &&
    a.domains === b.domains &&
    a.figure === b.figure &&
    a.windows === b.windows &&
    a.rails === b.rails &&
    a.headers === b.headers
  );
}

/** Two transforms over the same window of the same pieces, at the same place. */
function sameAxis(a: AxisTransform | undefined, b: AxisTransform | undefined): boolean {
  if (a === b) return true;
  if (!a || !b) return false;
  return (
    a.u0 === b.u0 &&
    a.u1 === b.u1 &&
    a.x0 === b.x0 &&
    a.width === b.width &&
    a.dir === b.dir &&
    samePieces(a.pieces, b.pieces)
  );
}

/** Two pieces arrays alike by value (a follower's are made again each frame). */
function samePieces(a: readonly number[], b: readonly number[]): boolean {
  if (a === b) return true;
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

/** Whether a section's kept layout is this frame's: its view parts, place, window, viewport by its layout, reads. */
/** Two grids' y place alike: the same pieces and window. */
export function sameGridY(a: PlaceContext['y'], b: PlaceContext['y']): boolean {
  if (a === b) return true;
  if (!a || !b || !samePieces(a.pieces, b.pieces) || a.axis.dir !== b.axis.dir) return false;
  const va = a.viewport;
  const vb = b.viewport;
  return va === vb || (!!va && !!vb && va.start === vb.start && va.end === vb.end);
}

export function reusable(
  kept: SectionCache,
  section: ResolvedSection,
  layout: SectionLayoutType<any>,
  stamp: ViewStamp,
  y: number,
  axis: AxisTransform | undefined,
  viewport: Viewport | undefined,
  data: StoreReader,
): boolean {
  return (
    kept.section === section &&
    kept.layout.layout === layout &&
    kept.y === y &&
    axis !== undefined &&
    sameStamp(kept.stamp, stamp) &&
    sameAxis(kept.axis, axis) &&
    layout.sameKey(kept.viewport, viewport) &&
    kept.reads.hold(data)
  );
}
