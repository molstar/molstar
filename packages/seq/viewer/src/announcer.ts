/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ReadonlyState, State } from '@molstar/core/reactive';
import { DomHost } from '@molstar/ui/dom';
import {
  itemText,
  type LayerLayout,
  LettersType,
  locationText,
  type SectionLayout,
  shownLetters,
  spanText,
} from '@molstar/seqx-core/render';
import type { SeqViewer } from './viewer.js';

// What the viewer says to a reader who does not see the canvas (_spec/README.md "Input", Announcements): one polite
// live region's text, made here DOM-free and shown by the mount. Three things are said: the keyboard cursor's place
// (its section, system, position, letter and the items there) as it moves, the hover's (when it rests: a pointer
// crossing the picture says nothing) and a selection's change (its size and extent, once). Everything goes through one
// rate limit: a message never follows another within `minGapMs`, the latest one waiting its turn (a held arrow key
// says where it is each `minGapMs`, not each key). Why one text a time: a live region queues what it is given, and a
// reader hears the queue after the pointer has moved on.

export interface AnnouncerOptions {
  /** How long the hover rests before it is said, ms. */
  readonly restMs?: number;
  /** The least time between two messages, ms. */
  readonly minGapMs?: number;
}

export const DefaultAnnouncerOptions = { restMs: 500, minGapMs: 700 } as const;

export class Announcer {
  private readonly _message = new State('');
  /** What the live region says now; a message equal to the last is said again only after the region was cleared. */
  readonly message: ReadonlyState<string> = this._message.asReadonly;
  private readonly _viewer: SeqViewer;
  private readonly _rest: number;
  private readonly _gap: number;
  private _last = Number.NEGATIVE_INFINITY;
  private _pending: string | undefined;
  private _gapTimer: unknown;
  private _restTimer: unknown;
  private readonly _stop: (() => void)[] = [];
  private _selected = '';
  private _disposed = false;

  constructor(viewer: SeqViewer, options: AnnouncerOptions = {}) {
    this._viewer = viewer;
    this._rest = options.restMs ?? DefaultAnnouncerOptions.restMs;
    this._gap = options.minGapMs ?? DefaultAnnouncerOptions.minGapMs;
    this._stop.push(
      viewer.selection.cursor.listen((c) => {
        // the pointer's own cursor is shown, not said: only the keyboard's moves are
        if (c?.keyboard) this.say(this.cursorText());
      }),
      viewer.hover.version.listen(() => this.hovered()),
      viewer.selection.value.listen(() => this.selectionChanged()),
    );
  }

  /** The text of a say, at once when the rate limit allows, else when it does (the latest of those waiting). */
  say(text: string): void {
    if (this._disposed || !text) return;
    const wait = this._last + this._gap - DomHost.now();
    if (wait <= 0) {
      this._pending = undefined;
      this._emit(text);
      return;
    }
    this._pending = text;
    if (this._gapTimer === undefined)
      this._gapTimer = DomHost.setTimeout(() => {
        this._gapTimer = undefined;
        const t = this._pending;
        this._pending = undefined;
        if (t && !this._disposed) this._emit(t);
      }, wait);
  }

  private _emit(text: string): void {
    this._last = DomHost.now();
    this._message.set(text);
  }

  // ---- the hover, when it rests ----

  private hovered(): void {
    DomHost.clearTimeout(this._restTimer);
    this._restTimer = undefined;
    if (!this._viewer.hover.active) return;
    this._restTimer = DomHost.setTimeout(() => {
      this._restTimer = undefined;
      this.say(this.hoverText());
    }, this._rest);
  }

  // ---- the selection, once a change ----

  private selectionChanged(): void {
    const shown = this._viewer.selection.value.value.shown;
    const text = shown.parts.length
      ? `Selected ${locationText(this._viewer.spec.model.value, shown)}`
      : this._selected
        ? 'Selection cleared'
        : '';
    if (text === this._selected && shown.parts.length) return;
    this._selected = shown.parts.length ? text : '';
    this.say(text);
  }

  // ---- the words ----

  private layout() {
    return this._viewer.frame?.layout;
  }

  private systemName(system: string): string {
    return this._viewer.spec.model.value.systems.get(system)?.name ?? system;
  }

  /** The items of a section at a position: its tracks' layers' items under the position's x. */
  private itemsAt(section: SectionLayout, position: number, only?: LayerLayout): string[] {
    const axis = section.axis;
    if (!axis) return [];
    const a = axis.x(position);
    const b = axis.x(position + 1);
    if (!Number.isFinite(a) || !Number.isFinite(b)) return [];
    const x = (a + b) / 2;
    const out: string[] = [];
    for (const t of section.tracks)
      for (const l of t.layers) {
        if (only && l !== only) continue;
        const data = l.layer.data?.system;
        for (const it of l.type.items?.(l.lay) ?? [])
          if (it.x0 <= x && x <= it.x1) out.push(itemText(this.systemName(data ?? section.section.axis.system), it));
      }
    return out;
  }

  /** The letter a letters layer of the section shows at a position (none while it is not loaded). */
  private letterAt(section: SectionLayout, position: number): string {
    const model = this._viewer.spec.model.value;
    for (const t of section.tracks)
      for (const l of t.layers) {
        if (l.type.name !== 'letters') continue;
        const r = shownLetters(
          LettersType.prepare(l.layer, model),
          section.section.axis.dir,
          this._viewer.store.store,
          [position],
        );
        if (r.text.trim()) return r.text.trim();
      }
    return '';
  }

  /** `protein, p53 12, letter L, items: Domain A, Domain B`. */
  placeText(section: SectionLayout, system: string, position: number, only?: LayerLayout): string {
    const parts = [`${section.id}, ${this.systemName(system)} ${spanText(position, position + 1)}`];
    const letter = this.letterAt(section, position);
    if (letter) parts.push(`letter ${letter}`);
    const items = this.itemsAt(section, position, only);
    if (items.length) parts.push(`items: ${items.join('; ')}`);
    return parts.join(', ');
  }

  cursorText(): string {
    const c = this._viewer.selection.cursor.value;
    const section = c && this.layout()?.sections.find((s) => s.section.address === c.section);
    return c && section ? this.placeText(section, c.system, c.position) : '';
  }

  hoverText(): string {
    const h = this._viewer.hover.current;
    const layout = this.layout();
    if (!layout || !this._viewer.hover.active) return '';
    const section = layout.sections.find((s) => s.section.address === h.section);
    if (!section) return '';
    const layer = section.tracks.flatMap((t) => t.layers).find((l) => l.key === h.layer);
    // an item: its own words; a cell: its place
    if (h.element && layer?.type.items) {
      const it = layer.type.items(layer.lay).find((d) => d.key === h.element!.key);
      if (it) return itemText(this.systemName(layer.layer.data?.system ?? section.section.axis.system), it);
    }
    if (h.system !== undefined && !Number.isNaN(h.position))
      return this.placeText(section, h.system, Math.floor(h.position));
    return '';
  }

  dispose(): void {
    this._disposed = true;
    for (const stop of this._stop) stop();
    this._stop.length = 0;
    DomHost.clearTimeout(this._gapTimer);
    DomHost.clearTimeout(this._restTimer);
  }
}
