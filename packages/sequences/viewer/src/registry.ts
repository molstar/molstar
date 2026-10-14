/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Registry entries and registration over both halves (seqspec:viewer#registration): the format's registries and the
// viewer's kinds,
// pure data grouped by kind, registered in a fixed order, conflicts checked before any change, an idempotent undo.

import {
  createRegistries,
  type ResolvedLayer,
  type ResolvedSection,
  register as registerSvs,
  type SvsRegistries,
  type SvsRegistryEntry,
  type ThemeEntry,
} from '@molstar/seq-core/format';
import type { Annotation, Field, Issue, SeqSet } from '@molstar/seq-core/model';
import type { AxisView } from './axis.js';
import type { SeqCommand } from './commands.js';
import type { RingLabelSink } from './draw/ring-labels.js';
import type { DrawView } from './draw/view.js';
import type { NodeKind, Skin } from './graph.js';
import type { Hit } from './hit.js';
import type { HoverState } from './hover.js';
import type { LabelContext, LabelProvider, LabelSink } from './labels.js';
import type { RenderProps } from './layout.js';
import type { NavigatorEntry } from './navigators.js';
import type { PresetEntry } from './presets.js';
import type { LineTable } from './runs.js';
import { isViewerOwn } from './spec.js';
import { type SeqState, state, type WritableSeqState } from './state.js';
import type { StreamLayer } from './streams.js';
import type { LegendSink, Theme, ThemeProvider } from './themes.js';

// The entry kinds. VC1 declares what registration and the graph read (a name, the format entry a kind pairs with,
// `prepare`, `create`, `content`); the phases that run each kind extend these with their full interfaces.

export interface NamedEntry {
  readonly name: string;
}
/** seqspec:viewer#extension-points; its name is its layer type, which a registration must hold. */
export interface LayerRenderer<S = unknown> extends NamedEntry {
  /**
   * The geometries it draws: the supported set per layout (seqspec:model#pieces-and-layout-maps); without `circular`, a
   * section on rings draws its extents as bars and says so (`svs.render.fallback`, seqspec:viewer#the-render-backend).
   */
  readonly geometries: readonly ('linear' | 'wrapped' | 'circular')[];
  /**
   * False: `prepare` never reads residues or alignment content, so a residue edit keeps its state
   * (seqspec:viewer#updates-and-status).
   */
  readonly content?: boolean;
  /** Headless and pure, once per input change (the graph's prepare stage). */
  prepare?(input: LayerInput, rc: PrepareContext): S;
  /** Data held by region: a pure description the stream manager loads from (seqspec:sources#the-viewers-side). */
  stream?(s: S): StreamLayer | null;
  dispose?(state: S): void;
  // The layout's and the hit's members (LY1), headless.
  /** The same pixels for the same positions wherever the view is; CS1 reads it. */
  readonly stable?: boolean | ((s: S, p: unknown) => boolean);
  /** Below it, `tiny` (or bars); CS1 reads it. */
  readonly minHeight?: number;
  /**
   * On rings, its labeled items for the section's labels outside the outermost ring: edge-rate, when the layout or the
   * view changes.
   */
  ringLabels?(s: S, p: unknown, out: RingLabelSink): void;
  /**
   * The names of the size channels it draws (seqspec:format#layer-types); a size theme on another is reported and not
   * drawn.
   */
  readonly sized?: readonly string[];
  /** Rows loaded by view: one layout entry for all rows. */
  readonly lazy?: boolean;
  /**
   * Free sizing, CSS px; `rowScale` is the section's row zoom (1: none), which only a renderer that `zooms` reads. An
   * integer keeps the row zoom allocation-free (a double returned is boxed). `width` is the data area's width and
   * `room` the height the section shows its rows in (a box's; 0: as tall as its content), whole CSS px, so a plot can
   * keep its aspect (`dot_plot`'s `aspect_ratio`); the layout asks again when the section's size changes.
   */
  height?(s: S, p: unknown, props: RenderProps, rowScale?: number, width?: number, room?: number): number;
  /** Its height follows the row zoom (`height`'s `rowScale`): the zoom rescales its entry with the rows. */
  readonly zooms?: boolean;
  /** The least a track holding it is (CSS px): a track `height` hint is raised to it, as a row never under a device px. */
  minTrackHeight?(s: S, p: unknown): number;
  /** Wrapped: band heights per line, from each line's block coordinate system runs (`axis.shownToRuns`). */
  lineHeights?(s: S, p: unknown, axis: AxisView, lines: LineTable, out: Float32Array): void;
  /** An expanded split row's height; default: the track's. */
  readonly splitHeight?: number;
  /**
   * Fills locations for an item under the point, `out.point` (no double crosses the call: V8 boxes a fractional one
   * it passes to a call it does not inline); false: none of its own.
   */
  hit?(v: AxisView, s: S, p: unknown, out: Hit): boolean;
  /**
   * Its glyphs overhang a linear block's end (an end boundary's, a point's on the last position): `hit` is also asked
   * just past the block, with the block and no position (`out.index` −1), so only drawn geometry holds the point.
   * Default false: past the end its layers are not asked.
   */
  readonly hitsPastEnd?: boolean;
  /**
   * Lines for the hovered element of its layers (after the generic providers'), from its prepared state
   * (seqspec:viewer#labels).
   */
  label?(h: HoverState, lc: LabelContext, out: LabelSink, s: S, p: unknown): void;
  /** Its channels' themes, in its type's channel order, and the prepared state, so it lists what the data has. */
  legend?(p: unknown, themes: readonly (Theme | null)[], out: LegendSink, s: S): void;
  /** False: never drawn in a miniature (a map, a navigator); default true. */
  readonly miniature?: boolean;
  /** A rows layer's labels in the header, which the composer draws: a data row's, a collapsed clade's wedge row's. */
  readonly rowModel?: { label(s: S, row: number): string; wedge?(s: S, node: number): string };
  /** Equal bands over the track (a matrix's states), labeled in the header. */
  readonly bands?: { count(s: S): number; label(s: S, band: number): string };
  /** A tree beside rows: the px row labels move right, a node under the point (as `hit`), the rows a node stands for. */
  readonly header?: {
    inset(s: S, p: unknown): number;
    hit(v: AxisView, s: S, p: unknown, out: Hit): boolean;
    span(s: S, node: number, out: Int32Array): void;
    /** Display rows [first, last] of the tree (`v.top`: row 0's), and a hovered node's clade on the overlay. */
    draw?(v: DrawView, s: S, p: unknown, first: number, last: number): void;
    hover?(v: DrawView, s: S, p: unknown, node: number, first: number, last: number): void;
  };
  // The drawing half: the renderer API only (draw/api.ts), no Canvas2D. The geometry is `v.geometry`.
  draw?(v: DrawView, s: S, p: unknown): void;
  /** Display rows [first, last] of a multi-row layer in one call (dense or lazy rows). */
  rows?(v: DrawView, s: S, p: unknown, first: number, last: number): void;
  /** Below `minHeight`; without it the engine draws the layer's placed extents as bars. */
  tiny?(v: DrawView, s: S, p: unknown): void;
  /** Numbers that follow the view, in the gutters. */
  gutter?(v: DrawView, s: S, p: unknown): void;
  /** What is lit, from `v.marks`. */
  overlay?(v: DrawView, s: S, p: unknown): void;
  /**
   * A data hover's row band on its track (the `hoverRow` behavior): its device rows [out[0], out[1]) over the entry
   * (`v.top`, `v.height`), as a matrix lights its hovered band; false: none (the ruler). Without it, the entry's box.
   */
  rowHover?(v: AxisView, s: S, p: unknown, h: HoverState, out: Int32Array): boolean;
  /** Navigators and map rows: the whole axis in miniature into the view's box; false: nothing drawn. */
  overview?(v: DrawView, s: S, p: unknown): boolean;
  /** The engine's, per prepared state and skin: made on the first draw, disposed when either changes. */
  readonly cache?: { create(s: S, skin: Skin): unknown; dispose?(c: unknown): void };
}
/** A data column's renderer by its type, an unknown type is drawn as `text`. */
export interface ColumnRenderer<S = unknown> extends NamedEntry {
  /** Default CSS px. */
  readonly width: number;
  /** The column's field (the data's own or a field provider's, null: none), per alignment row its `rows` annotation row. */
  prepare?(field: Field | null, rows: Int32Array, ann: Annotation | null, set: SeqSet): S;
  /** One row's cell: in `v.columnX` over `v.columnWidth`, over `v.top` and `v.height`; `row` the alignment row. */
  draw?(v: DrawView, s: S, row: number): void;
  /** The cell under `out.point` (as a renderer's `hit`): fills locations beyond the row's own; false: none. */
  hit?(v: AxisView, s: S, row: number, out: Hit): boolean;
}
/** What a renderer prepares from: one layer against its section's axis, never the composed view. */
export interface LayerInput {
  readonly set: SeqSet;
  readonly section: ResolvedSection;
  readonly layer: ResolvedLayer;
  readonly annotations: readonly Annotation[];
  readonly params: Readonly<Record<string, unknown>>;
}
export interface PrepareContext {
  readonly registry: SeqViewerRegistry;
  /** Reports what a prepare leaves out (a dot plot's block pair past its cap), on the layer's node. */
  readonly report?: (issue: Issue) => void;
}
export interface SeqViewerRegistryEntry {
  /** The format's half, as its `register` takes it. */
  readonly svs?: SvsRegistryEntry | readonly SvsRegistryEntry[];
  readonly nodeKinds?: readonly NodeKind[];
  readonly renderers?: readonly LayerRenderer[];
  readonly columns?: readonly ColumnRenderer[];
  readonly themes?: readonly ThemeProvider[];
  readonly navigators?: readonly NavigatorEntry[];
  readonly presets?: readonly PresetEntry[];
  readonly labels?: readonly LabelProvider[];
  /** A command's name is its id (seqspec:viewer#commands). */
  readonly commands?: readonly SeqCommand<any>[];
  readonly actions?: readonly NamedEntry[];
}

export interface ViewerEntries {
  nodeKinds: NodeKind;
  renderers: LayerRenderer;
  columns: ColumnRenderer;
  themes: ThemeProvider;
  navigators: NavigatorEntry;
  presets: PresetEntry;
  labels: LabelProvider;
  commands: SeqCommand<unknown>;
  actions: NamedEntry;
}
export type ViewerKind = keyof ViewerEntries;
export type ViewerEntryOf<K extends ViewerKind> = ViewerEntries[K];

/** The interface's order: within one entry, kinds register in it. */
const KINDS: readonly ViewerKind[] = [
  'nodeKinds',
  'renderers',
  'columns',
  'themes',
  'navigators',
  'presets',
  'labels',
  'commands',
  'actions',
];

const nameOf = (kind: ViewerKind, item: unknown) =>
  kind === 'commands' ? (item as { id: string }).id : (item as { name: string }).name;

/** A renderer's or theme provider's format entry, registered or in this call: none, or one missing a kind it serves. */
function unpaired(kind: ViewerKind, item: unknown, reg: SvsRegistries, svs: readonly SvsRegistryEntry[]): boolean {
  const name = nameOf(kind, item);
  if (kind === 'renderers')
    return !reg.get('layerTypes', name) && !svs.some((s) => s.layerTypes?.some((e) => e.name === name));
  if (kind !== 'themes') return false;
  const e = reg.get('themes', name) ?? svs.flatMap((s) => s.themes ?? []).find((t: ThemeEntry) => t.name === name);
  return !e || (item as ThemeProvider).kinds.some((k) => !e.params[k]);
}

export interface SeqViewerRegistry {
  /** The format's registries, passed to its functions. */
  readonly svs: SvsRegistries;
  /** Bumped once by every registration and undo. */
  readonly version: SeqState<number>;
  get<K extends ViewerKind>(kind: K, name: string): ViewerEntryOf<K> | undefined;
  /** Registration order; cached until a change. */
  list<K extends ViewerKind>(kind: K): readonly ViewerEntryOf<K>[];
}

/** Registration conflicts in either half, listed together; nothing was registered. */
export class SeqRegistrationError extends Error {
  readonly conflicts: readonly string[];
  constructor(conflicts: readonly string[]) {
    super(`register: ${conflicts.length} conflict(s), nothing registered:\n- ${conflicts.join('\n- ')}`);
    this.name = 'SeqRegistrationError';
    this.conflicts = conflicts;
  }
}

/**
 * The viewer's own format entries: the only ones `register` lets name its `seqv` namespace (seqspec:format#extensions),
 * known by their brand (`isViewerOwn`), so the entries of another copy of the viewer in the page are its own too.
 */
class ViewerOwned extends Set<SvsRegistryEntry> {
  override has(e: SvsRegistryEntry): boolean {
    return isViewerOwn(e);
  }
}
const SEQV: ReadonlySet<SvsRegistryEntry> = new ViewerOwned();

/** The format's error lists its conflicts as "- " lines. */
function svsConflicts(e: unknown): string[] {
  const message = e instanceof Error ? e.message : String(e);
  const lines = message.split('\n- ').slice(1);
  return lines.length > 0 ? lines.map((l) => `svs ${l}`) : [`svs: ${message}`];
}

export class ViewerRegistry implements SeqViewerRegistry {
  readonly svs = createRegistries();
  private readonly versionState: WritableSeqState<number> = state(0);
  private readonly maps = new Map(KINDS.map((k) => [k, new Map<string, unknown>()]));
  private readonly cache = new Map<ViewerKind, readonly unknown[]>();
  /** Undos of live registrations, so dispose clears both halves. */
  private readonly undos = new Set<() => void>();

  get version(): SeqState<number> {
    return this.versionState;
  }
  get<K extends ViewerKind>(kind: K, name: string) {
    return this.maps.get(kind)?.get(name) as ViewerEntryOf<K> | undefined;
  }
  list<K extends ViewerKind>(kind: K) {
    let l = this.cache.get(kind);
    if (!l) {
      l = [...(this.maps.get(kind)?.values() ?? [])];
      this.cache.set(kind, l);
    }
    return l as readonly ViewerEntryOf<K>[];
  }

  private changed() {
    this.cache.clear();
    this.versionState.set(this.versionState.value + 1);
  }

  register(entry: SeqViewerRegistryEntry | readonly SeqViewerRegistryEntry[]): () => void {
    const entries: readonly SeqViewerRegistryEntry[] = Array.isArray(entry) ? entry : [entry as SeqViewerRegistryEntry];
    const svs: SvsRegistryEntry[] = [];
    for (const e of entries) if (e.svs) svs.push(...(Array.isArray(e.svs) ? e.svs : [e.svs as SvsRegistryEntry]));

    // The viewer's kinds first, so their conflicts are known before the format's call changes anything.
    const added: [Map<string, unknown>, string, unknown][] = [];
    const pending = new Map<string, unknown>();
    const conflicts: string[] = [];
    for (const e of entries) {
      for (const kind of KINDS) {
        const map = this.maps.get(kind) as Map<string, unknown>;
        for (const item of (e[kind] ?? []) as readonly unknown[]) {
          const name = nameOf(kind, item);
          const held = map.get(name) ?? pending.get(`${kind} ${name}`);
          if (held === item) continue;
          if (held !== undefined) {
            conflicts.push(`${kind}: "${name}" is ${map.has(name) ? 'already registered' : 'listed twice'}`);
            continue;
          }
          if (unpaired(kind, item, this.svs, svs))
            conflicts.push(`${kind}: "${name}" has no format entry in either registration`);
          pending.set(`${kind} ${name}`, item);
          added.push([map, name, item]);
        }
      }
    }

    let undoSvs: (() => void) | null = null;
    if (svs.length > 0) {
      try {
        undoSvs = registerSvs(this.svs, svs, { seqv: SEQV });
      } catch (e) {
        throw new SeqRegistrationError([...svsConflicts(e), ...conflicts]);
      }
    }
    if (conflicts.length > 0) {
      undoSvs?.();
      throw new SeqRegistrationError(conflicts);
    }

    for (const [map, name, item] of added) map.set(name, item);
    this.changed();
    let undone = false;
    const undo = () => {
      if (undone) return;
      undone = true;
      this.undos.delete(undo);
      for (let i = added.length - 1; i >= 0; i--) {
        const [map, name, item] = added[i];
        if (map.get(name) === item) map.delete(name);
      }
      undoSvs?.();
      this.changed();
    };
    this.undos.add(undo);
    return undo;
  }

  /** Completes the version, then undoes every live registration, newest first, with no emission. */
  dispose() {
    this.versionState.complete();
    for (const undo of [...this.undos].reverse()) undo();
  }
}
