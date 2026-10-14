/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Labels (seqspec:viewer#labels): label providers, in priority order, write into a reused sink; one `HoverLabel` per
// hovered element is the only allocation. When to compute and hide it is the `hoverLabel` behavior's. A label is
// computed again only when the element, the registry or the view changed: three plain fields, no work on a move.

import { pathTo, type StreamedData } from "@molstar/seq-core/format";
import {
  type ComposedPath,
  CoordinateSystemKind as K,
  Numbering,
  type SeqIndex,
  type SeqSet,
} from "@molstar/seq-core/model";
import { sameList } from "./axis.ts";
import type { ViewModel } from "./graph.ts";
import { HoverKind } from "./hit.ts";
import type { HoverState } from "./hover.ts";
import { type HeaderLegend, headerLegend } from "./legend.ts";
import type { SeqViewerRegistry } from "./registry.ts";
import { type SeqState, state, type WritableSeqState } from "./state.ts";
import type { StreamView } from "./streams.ts";

export interface LabelSink {
  /** The first title wins (higher priorities write first). */
  title(text: string): void;
  line(key: string, value: string): void;
  replaceTitle(text: string): void;
}
/** A numbering of another molecule that mappings reach from a unit, with the path to it. */
export interface ReachedNumbering {
  readonly n: Numbering;
  readonly path: ComposedPath;
}
export interface LabelContext {
  readonly set: SeqSet;
  readonly view: ViewModel;
  readonly section: number;
  /** The viewer's registry: renderers' label entries, and the format's field providers (`registry.svs`). */
  readonly registry: SeqViewerRegistry;
  format(n: Numbering, i: SeqIndex): string;
  /** Layer k's regions held by view (`StreamView`); null: none (not streamed). */
  stream(layer: number): StreamView | null;
  /**
   * The label shows data still in flight (a streamed row's tile): it is made again when a region arrives, not only
   * when the element, the registry or the view changes.
   */
  waiting(): void;
  /** Per unit, the named numberings of the other molecules its mappings reach; composed once per set. */
  reached(unit: number): readonly ReachedNumbering[];
  /** A lazy molecule's letters at positions (`StreamManager.lettersAt`; ST7's codon); none: nothing held. */
  letters?(data: StreamedData, s: number, at: Int32Array, n: number, out: Uint8Array): number;
}
export interface LabelProvider {
  readonly name: string;
  readonly priority: number;
  /** By kind and region; no allocation. */
  applies(h: HoverState): boolean;
  label(h: HoverState, lc: LabelContext, out: LabelSink): void;
}
export interface HoverLabel {
  readonly title: string;
  readonly lines: readonly { readonly key: string; readonly value: string }[];
  /**
   * A track header's hover (no tree node): the track's legend, swatches and ramps (seqspec:interaction#the-legend);
   * absent: none.
   */
  readonly legend?: HeaderLegend;
}
export interface LabelManager {
  /** Providers in priority order; for the live hover, once per hovered element. */
  compute(h: HoverState): HoverLabel | null;
  readonly current: WritableSeqState<HoverLabel | null>;
  /** Registered, by priority. */
  readonly providers: SeqState<readonly string[]>;
}

export interface LabelHost {
  readonly registry: SeqViewerRegistry;
  readonly hover: { readonly current: HoverState; readonly version: SeqState<number> };
  view(): ViewModel | null;
  /** A layer's regions by its ref (`ctx.streams.get`), and the manager's version (bumped by every arrival). */
  stream(layer: string): StreamView | null;
  streamVersion(): number;
  /** A lazy molecule's letters, asked when not held (`StreamManager.lettersAt`). */
  letters?(data: StreamedData, s: number, at: Int32Array, n: number, out: Uint8Array): number;
  /** A hovered header's full text when the canvas draws it cut (`cutHeaderText`); null: drawn whole, or none. */
  cutText?(h: HoverState): string | null;
}

const name = (set: SeqSet, id: number) => (id >= 0 ? set.strings.get(id) : "");
/**
 * The numberings of the other molecules mappings reach from a unit (an isoform's UniProt numbering), with paths. Only
 * a named numbering is worth a path: composing one is O(the molecules' lengths), and a genome's chromosome reaches
 * hundreds of proteins numbered by ordinal only (40 ms per chromosome's first hover).
 */
function reach(set: SeqSet, unit: number): ReachedNumbering[] {
  const out: ReachedNumbering[] = [];
  const own = set.units.molecule[unit];
  for (const m of set.molecules) {
    if (!m || m.id === own || m.alphabet === null) continue;
    if (!m.numberings.items.some((n) => name(set, n.name) !== "ordinal")) continue;
    const p = pathTo(set, { kind: K.unit, id: unit, row: -1 }, { kind: K.molecule, id: m.id, row: -1 });
    if (!p || p === "ambiguous") continue;
    for (const n of m.numberings.items) if (name(set, n.name) !== "ordinal") out.push({ n, path: p });
  }
  return out;
}

export class Labels implements LabelManager {
  readonly current = state<HoverLabel | null>(null);
  readonly providers = state<readonly string[]>([], sameList);
  private sorted: readonly LabelProvider[] = [];
  private readonly host: LabelHost;
  private titleText: string | null = null;
  private readonly keys: string[] = [];
  private readonly values: string[] = [];
  private readonly counts: number[] = [];
  /** What the live label was computed for: the hover's version, the registry's, the view's. */
  private cachedHover = -1;
  private cachedRegistry = -1;
  private cachedView: ViewModel | null = null;
  private cachedLabel: HoverLabel | null = null;
  /** Whether the label waits on data in flight (`LabelContext.waiting`), and the stream version it was made at. */
  private waits = false;
  private cachedStream = -1;
  /** The reached numberings per unit of one set, dropped with it. */
  private reachedSet: SeqSet | null = null;
  private readonly reachedBy = new Map<number, readonly ReachedNumbering[]>();
  private readonly lc: { set: SeqSet; view: ViewModel; section: number; registry: SeqViewerRegistry } & LabelContext = {
    set: null as unknown as SeqSet,
    view: null as unknown as ViewModel,
    section: -1,
    registry: null as unknown as SeqViewerRegistry,
    format: (n, i) => Numbering.format(n, i, this.lc.set.strings),
    reached: (unit) => this.reachedFrom(unit),
    stream: (layer) => this.host.stream(this.lc.set.strings.get(this.lc.view.resolved.layers.ref[layer])),
    waiting: () => {
      this.waits = true;
    },
    letters: (data, s, at, n, out) => this.host.letters?.(data, s, at, n, out) ?? 0,
  };
  private readonly sink: LabelSink = {
    title: (t) => {
      this.titleText ??= t;
    },
    line: (k, v) => {
      this.keys.push(k);
      this.values.push(v);
    },
    replaceTitle: (t) => {
      this.titleText = t;
    },
  };

  constructor(host: LabelHost) {
    this.host = host;
    this.lc.registry = host.registry;
    host.registry.version.subscribe(() => {
      this.sorted = [...host.registry.list("labels")].sort((a, b) => b.priority - a.priority);
      this.providers.set(this.sorted.map((p) => p.name));
    });
  }

  compute(h: HoverState): HoverLabel | null {
    const live = h === this.host.hover.current;
    const version = this.host.hover.version.value;
    const registry = this.host.registry.version.value;
    const view = this.host.view();
    const streamed = this.host.streamVersion();
    if (
      live &&
      version === this.cachedHover &&
      registry === this.cachedRegistry &&
      view === this.cachedView &&
      !(this.waits && streamed !== this.cachedStream)
    )
      return this.cachedLabel;
    this.waits = false;
    this.titleText = null;
    this.keys.length = 0;
    this.values.length = 0;
    if (view) {
      if (view.set !== this.reachedSet) {
        this.reachedSet = view.set;
        this.reachedBy.clear();
      }
      this.lc.set = view.set;
      this.lc.view = view;
      this.lc.section = h.section;
      for (const p of this.sorted) if (p.applies(h)) p.label(h, this.lc, this.sink);
    }
    // A track's header (not a tree's): the full text where the canvas cuts it, then the track's legend (one object
    // per track and view, made once), titled by the track when the text is whole.
    // Labels on (a provider registered): with none, a hover has no label at all.
    const header = view && this.sorted.length > 0 && h.kind === HoverKind.header && h.treeNode < 0 && h.track >= 0;
    const legend = header ? headerLegend(this.host.registry, view, h.track) : null;
    const full = header ? (this.host.cutText?.(h) ?? null) : null;
    if (full !== null) this.titleText = full;
    else if (legend && view) {
      const t = view.resolved.tracks;
      this.titleText ??= t.label[h.track] >= 0 ? view.set.strings.get(t.label[h.track]) : "";
    }
    const label = this.titleText === null && this.keys.length === 0 && !legend ? null : this.group(legend);
    if (live) {
      this.cachedHover = version;
      this.cachedRegistry = registry;
      this.cachedView = view;
      this.cachedStream = streamed;
      this.cachedLabel = label;
    }
    return label;
  }

  private reachedFrom(unit: number): readonly ReachedNumbering[] {
    const held = this.reachedBy.get(unit);
    if (held) return held;
    const out = reach(this.lc.set, unit);
    this.reachedBy.set(unit, out);
    return out;
  }

  /** Equal lines once, `× n` when repeated (a merged copy group's units): the label's own objects only. */
  private group(legend: HeaderLegend | null): HoverLabel {
    const lines: { key: string; value: string }[] = [];
    const n = this.counts;
    n.length = 0;
    for (let i = 0; i < this.keys.length; i++) {
      let j = 0;
      while (j < lines.length && (lines[j].key !== this.keys[i] || lines[j].value !== this.values[i])) j++;
      if (j < lines.length) n[j]++;
      else {
        lines.push({ key: this.keys[i], value: this.values[i] });
        n.push(1);
      }
    }
    for (let j = 0; j < lines.length; j++) if (n[j] > 1) lines[j].value = `${lines[j].value} × ${n[j]}`;
    return legend ? { title: this.titleText ?? "", lines, legend } : { title: this.titleText ?? "", lines };
  }

  dispose() {
    this.current.complete();
    this.providers.complete();
  }
}
