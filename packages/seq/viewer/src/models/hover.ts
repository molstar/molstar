/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover's model (_spec/README.md "HoverModel"): a method, never a command, run per pointer move. `set` copies the
// hit into one reused record and moves a version, allocating nothing; one owner at a time, a stale owner's clear does
// nothing (EC-141). The hover's Location and its projection onto each linked section are computed lazily, when a
// section is drawn, once per version: closed first (a nucleotide lights its whole codon and its residue, EC-30,
// EC-144), then carried along the view's links (paths found under the copy rule, EC-25).

import { type ReadonlyState, ReactiveModel, State } from '@molstar/core/reactive';
import type { Hit as RenderHit } from '@molstar/seqx-core/render';
import { type ElementKey, EmptyLocation, type Location, location } from '@molstar/seqx-core/query';
import { projectRanges, type Path } from '@molstar/seqx-core/resolve';
import { projectOntoView } from '@molstar/seqx-core/view';
import { DefaultViewerConfig, type ViewerConfig } from '@molstar/seqx-viewer/config';
import type { SpecModel } from './spec.js';
import type { ViewModel } from './view.js';

/** What a hit's target is, for the gesture bindings (_spec/README.md "Bindings"): the granularity comes from the hit. */
export type HitTarget =
  | 'cell'
  | 'item'
  | 'row'
  | 'column'
  | 'track'
  | 'empty'
  | 'axis'
  | 'break'
  | 'filter'
  | 'window'
  | 'tree-node'
  | 'summary-row';

/**
 * What the pointer is over, as the layout's hit test gives it: the layer, an element key and its extent, or a
 * position, on a system. `section` and `target` may be given; otherwise they are read from the layer's address and
 * from what was hit (an element: an item; a position: a cell; neither: empty space).
 */
export interface Hit
  extends Omit<RenderHit, 'layer' | 'section' | 'target'>,
    Partial<Pick<RenderHit, 'layer' | 'section'>> {
  readonly target?: HitTarget;
}

/** The kind of target a hit is. */
export function hitTarget(hit: Hit | undefined): HitTarget {
  if (!hit) return 'empty';
  if (hit.target) return hit.target;
  if (hit.element) return 'item';
  if (hit.position !== undefined) return 'cell';
  return 'empty';
}

/** The hover, as one record reused for every move: read it, never keep it (copy what you need). */
export interface HoverRecord {
  readonly owner: unknown;
  readonly layer: string | undefined;
  readonly section: string | undefined;
  readonly system: string | undefined;
  /** NaN when the hover is no position. */
  readonly position: number;
  readonly element: ElementKey | undefined;
  readonly extent: Location | undefined;
  /** The position on a second axis (a rows cell's row, a grid's y) and its system: NaN and undefined when none. */
  readonly ySystem: string | undefined;
  readonly yPosition: number;
  /** What the hit was (a tree node, a summary row: the clade's hover), as the hit said; undefined when it said nothing. */
  readonly target: HitTarget | undefined;
}

type MutableRecord = { -readonly [K in keyof HoverRecord]: HoverRecord[K] };

const keysEqual = (a: ElementKey | undefined, b: ElementKey | undefined) =>
  a === b || (!!a && !!b && a.owner === b.owner && a.key === b.key);

export class HoverModel extends ReactiveModel {
  private readonly _version = new State(0);
  private readonly _record: MutableRecord = {
    owner: undefined,
    layer: undefined,
    section: undefined,
    system: undefined,
    position: Number.NaN,
    element: undefined,
    extent: undefined,
    ySystem: undefined,
    yPosition: Number.NaN,
    target: undefined,
  };
  /** The version, the view and the model the caches are for. */
  private _cacheKey: [number, unknown, unknown] = [-1, undefined, undefined];
  private _location: Location | undefined;
  private _closed: Location | undefined;
  private _all: Location | undefined;
  private readonly _projected = new Map<string, Location>();
  private readonly _buffer: number[] = [];
  private readonly _view: ViewModel;
  private readonly _spec: SpecModel;
  readonly config: ViewerConfig;
  /** Moves on every change of the hover: what draws it listens to this. */
  readonly version: ReadonlyState<number> = this._version.asReadonly;
  /** The current hover, one reused record. */
  readonly current: HoverRecord = this._record;

  constructor(spec: SpecModel, view: ViewModel, config: ViewerConfig = DefaultViewerConfig) {
    super();
    this._spec = spec;
    this._view = view;
    this.config = config;
  }

  get owner(): unknown {
    return this._record.owner;
  }

  /** Whether anything is hovered. */
  get active(): boolean {
    return this._record.owner !== undefined;
  }

  /**
   * Hovers what `hit` names for `owner`, taking the hover from any other owner; no hit clears it if `owner` holds it.
   * Copies the hit (no object is kept or made); false when nothing changed. An element is the same hover by its key
   * (a layer may make its extent again on each move); an undefined owner is refused (it could never clear it).
   */
  set(hit: Hit | undefined, owner: unknown): boolean {
    if (owner === undefined) return false;
    if (!hit) return this.clear(owner);
    const r = this._record;
    const position = hit.position ?? Number.NaN;
    const yPosition = hit.yPosition ?? Number.NaN;
    const section = hit.section ?? this._sectionOf(hit.layer);
    if (
      r.owner === owner &&
      r.system === hit.system &&
      r.layer === hit.layer &&
      r.section === section &&
      (r.position === position || (Number.isNaN(r.position) && Number.isNaN(position))) &&
      r.ySystem === hit.ySystem &&
      r.target === hit.target &&
      (r.yPosition === yPosition || (Number.isNaN(r.yPosition) && Number.isNaN(yPosition))) &&
      keysEqual(r.element, hit.element) &&
      (hit.element !== undefined || r.extent === hit.extent)
    )
      return false;
    r.owner = owner;
    r.layer = hit.layer;
    r.section = section;
    r.system = hit.system;
    r.position = position;
    r.element = hit.element;
    r.extent = hit.extent;
    r.ySystem = hit.ySystem;
    r.yPosition = yPosition;
    r.target = hit.target;
    this._version.set(this._version.value + 1);
    return true;
  }

  /** Clears the hover when `owner` holds it; a stale owner's clear (a late leave) does nothing (EC-141). */
  clear(owner: unknown): boolean {
    const r = this._record;
    if (r.owner === undefined || r.owner !== owner) return false;
    r.owner = r.layer = r.section = r.system = r.element = r.extent = r.ySystem = r.target = undefined;
    r.position = r.yPosition = Number.NaN;
    this._version.set(this._version.value + 1);
    return true;
  }

  /**
   * The hover's own Location: the element's extent, else its position, and its position on a second axis (a grid's
   * cell lights both: the crosshair's two positions on their sections); undefined when nothing is hovered.
   */
  location(): Location | undefined {
    this._fresh();
    if (this._location === undefined && this._record.owner !== undefined) {
      const r = this._record;
      if (r.extent) this._location = r.extent;
      else {
        const parts: { system: string; ranges: number[] }[] = [];
        if (!Number.isNaN(r.position) && r.system !== undefined)
          parts.push({ system: r.system, ranges: [Math.floor(r.position), Math.floor(r.position) + 1] });
        if (!Number.isNaN(r.yPosition) && r.ySystem !== undefined)
          parts.push({ system: r.ySystem, ranges: [Math.floor(r.yPosition), Math.floor(r.yPosition) + 1] });
        this._location = parts.length ? location({ parts }) : EmptyLocation;
      }
    }
    return this._location;
  }

  /**
   * The hover as a section draws it, computed when asked and kept for this version: the closed Location's parts on the
   * section's system, and the others carried along the view's links (when `hover.projectLinked`). Empty when the
   * section is not linked; undefined when nothing is hovered.
   */
  projected(section: string): Location | undefined {
    const own = this.location();
    if (!own) return undefined;
    const hit = this._projected.get(section);
    if (hit) return hit;
    const view = this._view.resolved.value;
    const target = view?.sections.find((s) => s.address === section);
    if (!view || !target) return EmptyLocation;
    const system = target.axis.system;
    const closed = this._closedOf(own);
    const parts: { system: string; ranges: ArrayLike<number> }[] = [];
    for (const part of closed.parts) {
      if (part.system === system) parts.push({ system, ranges: part.ranges });
      else if (this.config.hover.projectLinked) {
        const path = this._pathOf(part.system, section);
        if (path) parts.push({ system, ranges: projectRanges(path, part.ranges, this._buffer).slice() });
      }
    }
    const out = parts.length ? location({ parts }) : EmptyLocation;
    this._projected.set(section, out);
    return out;
  }

  /**
   * The hover as every section draws it (`Emphasis.projected.hover`): the closed Location carried onto the view's
   * sections' systems, one Location keyed by system, computed when asked and kept for this version; undefined when
   * nothing is hovered or no view is shown.
   */
  projectedAll(): Location | undefined {
    const own = this.location();
    if (!own) return undefined;
    if (this._all) return this._all;
    const view = this._view.resolved.value;
    if (!view) return undefined;
    this._all = projectOntoView(this._spec.model.value, view, own, {
      own: true,
      close: this.config.hover.closeCodons,
      linked: this.config.hover.projectLinked,
    });
    return this._all;
  }

  // ---- internals ----

  /** Drops what was computed for an older hover, view or model. */
  private _fresh(): void {
    const k = this._cacheKey;
    const v = this._version.value;
    const view = this._view.resolved.value;
    const model = this._spec.model.value;
    if (k[0] === v && k[1] === view && k[2] === model) return;
    k[0] = v;
    k[1] = view;
    k[2] = model;
    this._location = undefined;
    this._closed = undefined;
    this._all = undefined;
    this._projected.clear();
  }

  /** The hover closed among the view's systems: once per version. */
  private _closedOf(own: Location): Location {
    if (this._closed) return this._closed;
    const systems = this._view.resolved.value?.sections.map((s) => s.axis.system) ?? [];
    this._closed = this.config.hover.closeCodons ? this._spec.model.value.close(own, systems) : own;
    return this._closed;
  }

  /** The view's link from a section on `system` to `section`; else the model's static path (a data system). */
  private _pathOf(system: string, section: string): Path | undefined {
    const view = this._view.resolved.value;
    if (view)
      for (const l of view.links) {
        if (l.to !== section) continue;
        const from = view.sections.find((s) => s.address === l.from);
        if (from?.axis.system === system) return l.path;
      }
    const to = view?.sections.find((s) => s.address === section);
    if (!to || view?.sections.some((s) => s.axis.system === system)) return undefined;
    return this._spec.model.value.path(system, to.axis.system).path;
  }

  /** The section a layer is in: the view's section whose address starts its address. */
  private _sectionOf(layer: string | undefined): string | undefined {
    if (layer === undefined) return undefined;
    for (const s of this._view.resolved.value?.sections ?? [])
      if (
        layer.startsWith(s.address) &&
        (layer.length === s.address.length || layer.charCodeAt(s.address.length) === 47) // '/'
      )
        return s.address;
    return undefined;
  }
}
