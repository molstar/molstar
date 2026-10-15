/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The public input streams (seqspec:interaction#streams): what a consumer outside the viewer (a UI, Prism, a host) sees
// of input, fed by the dispatcher after the behaviors ran. Observe-only: a subscriber watches a gesture, never claims
// it. A payload is a fresh object a consumer may keep, in Mol*'s terms for buttons and modifiers; it is built only when
// its stream has a subscriber (`observed`), so the input path costs nothing more unobserved.

import type { ButtonsType, ModifiersKeys } from '@molstar/core/util/input/input-observer';
import { filter, type Observable, Subject } from 'rxjs';
import type { Assembled } from './graph.js';
import { type Hit, HitRegion } from './hit.js';
import type { ClickEvt, DragEvt, KeyInput, PointerInput, WheelEvt } from './input.js';

/** The DOM event behind a payload, structurally (the headless build has no DOM types). */
export interface SeqNativeEvent {
  readonly type: string;
  readonly defaultPrevented: boolean;
  preventDefault(): void;
}
/** The surface's regions (`HitRegion` by name; `none`: the surface, nothing under the point). */
export type SeqRegionKind = keyof typeof HitRegion;
/** The default UIs' parts that report input (their roots' `data-seqv-part`). */
export type SeqChromeKind =
  | 'viewTab'
  | 'sectionHeader'
  | 'sectionSwitch'
  | 'trackAction'
  | 'hiddenChip'
  | 'rowControl'
  | 'legend'
  | 'snapshotPlayer'
  | 'search'
  | 'export'
  | 'issues'
  | 'hoverCard'
  | 'tooltip'
  | 'description'
  | 'annotations';
export type SeqPartKind = SeqRegionKind | SeqChromeKind;
/** Where: a region of a surface or a part of the chrome, with refs (names, never indices). */
export interface SeqPart {
  readonly kind: SeqPartKind;
  readonly section: string | null;
  readonly track: string | null;
  readonly group: string | null;
  /** A chrome node's `data-seqv-key` (a tab's, a legend entry's); null on a surface. */
  readonly key: string | null;
}
/** What, on a surface: the layer and lane, the position (block, model index), the display row, an item. -1: none. */
export interface SeqAt {
  readonly layer: string | null;
  readonly lane: number;
  readonly block: number;
  readonly index: number;
  readonly row: number;
  /** A split track's row (its value's; -2 its label's). */
  readonly split: number;
  readonly item: { readonly annotation: string; readonly item: number; readonly id: string | null } | null;
}
/** Surface px (a chrome part's own, from its root's corner), and the page's client px. NaN: not known (a key's). */
export interface SeqPlace {
  readonly x: number;
  readonly y: number;
  readonly clientX: number;
  readonly clientY: number;
}
export interface SeqInputEvent {
  readonly stream: SeqStreamName;
  readonly origin: 'pointer' | 'keyboard' | 'chrome';
  readonly part: SeqPart;
  /** Null on chrome and on a key. */
  readonly at: SeqAt | null;
  readonly place: SeqPlace;
  readonly buttons: ButtonsType;
  readonly button: ButtonsType.Flag;
  readonly modifiers: Readonly<ModifiersKeys>;
  readonly time: number;
  /** The behavior that took the gesture; '' none. */
  readonly owner: string;
  /**
   * The DOM event, during the call only (a consumer may `preventDefault`); after it, and without a DOM, a stand-in.
   * Headless, the stand-in of a call records `preventDefault` (what `menu` and `key` report as handled).
   */
  readonly native: SeqNativeEvent;
}
export interface SeqClickEvent extends SeqInputEvent {
  /** 2: a double click, 3: a triple. */
  readonly count: number;
}
export interface SeqDragEvent extends SeqInputEvent {
  readonly phase: 'start' | 'move' | 'end' | 'cancel';
  readonly dx: number;
  readonly dy: number;
  readonly totalX: number;
  readonly totalY: number;
  /** Where the press was. */
  readonly start: SeqPart;
}
export interface SeqWheelEvent extends SeqInputEvent {
  readonly deltaX: number;
  readonly deltaY: number;
  /** A behavior acted on it. */
  readonly handled: boolean;
}
export interface SeqKeyEvent extends SeqInputEvent {
  readonly key: string;
  readonly code: string;
  readonly handled: boolean;
}
export interface SeqStreams {
  /** The element under the pointer changed. */
  readonly hover: SeqInputEvent;
  readonly press: SeqInputEvent;
  readonly click: SeqClickEvent;
  /** start, move, end, cancel; moves and the end only of a drag a behavior took. */
  readonly drag: SeqDragEvent;
  /** The browser's context-menu request (seqspec:interaction#the-track-menu). */
  readonly menu: SeqInputEvent;
  readonly wheel: SeqWheelEvent;
  readonly key: SeqKeyEvent;
  readonly leave: SeqInputEvent;
  /** A track header activated (the `headers` behavior; the keyboard's Enter on a header). */
  readonly activate: SeqInputEvent;
}
export type SeqStreamName = keyof SeqStreams;
/** What `on` keeps: each given field must match (a list of kinds: any of them). */
export interface SeqInputScope {
  readonly part?: SeqPartKind | readonly SeqPartKind[];
  readonly section?: string;
  readonly track?: string;
}
/** A UI's chrome input, for `SeqInteraction.chrome`. */
export type SeqChromeStream = 'hover' | 'click' | 'menu' | 'leave';

export type SeqInteraction = { readonly [K in SeqStreamName]: Observable<SeqStreams[K]> } & {
  /** The stream, kept to `scope`; unscoped, the stream itself. */
  on<K extends SeqStreamName>(stream: K, scope?: SeqInputScope): Observable<SeqStreams[K]>;
  /** A stream has a subscriber: its payloads are built. */
  observed(stream: SeqStreamName): boolean;
  /** A UI reports its chrome's input (the default UIs' delegated listener; a host's own parts may). */
  chrome(
    stream: SeqChromeStream,
    part: SeqPart,
    native: SeqNativeEvent,
    place: SeqPlace,
    opts?: { count?: number; buttons?: number; button?: number; modifiers?: number },
  ): void;
};

/** After the call, and where there is no DOM event. */
const NO_NATIVE: SeqNativeEvent = Object.freeze({ type: '', defaultPrevented: false, preventDefault() {} });
/** A headless call's stand-in: records `preventDefault`. Built only with a payload. */
const standIn = (type: string): SeqNativeEvent => {
  const e = {
    type,
    defaultPrevented: false,
    preventDefault() {
      e.defaultPrevented = true;
    },
  };
  return e;
};
const REGIONS = Object.keys(HitRegion) as SeqRegionKind[];
/** The DOM's `button` index as Mol*'s flag: primary, auxiliary (middle), secondary, fourth, fifth. */
const BUTTON = [1, 4, 2, 8, 16];
/** `Modifier` bits (shift 1, ctrl 2, alt 4, meta 8) as Mol*'s `ModifiersKeys`, one frozen object each. */
const MODIFIERS: readonly Readonly<ModifiersKeys>[] = Array.from({ length: 16 }, (_, m) =>
  Object.freeze({ shift: !!(m & 1), control: !!(m & 2), alt: !!(m & 4), meta: !!(m & 8) }),
);
const PHASES = ['start', 'move', 'end', 'cancel'] as const;
const NAMES: readonly SeqStreamName[] = [
  'hover',
  'press',
  'click',
  'drag',
  'menu',
  'wheel',
  'key',
  'leave',
  'activate',
];

export interface InteractionHost {
  sectionRef(section: number): string | null;
  assembled(): Assembled | null;
  now(): number;
}

type Writable<T> = { -readonly [K in keyof T]: T[K] };
type Payload = Writable<SeqInputEvent> & Record<string, unknown>;
type DomLike = {
  readonly timeStamp?: number;
  readonly clientX?: number;
  readonly clientY?: number;
  readonly buttons?: number;
  readonly currentTarget?: unknown;
  readonly target?: unknown;
};
const num = (v: unknown, or: number) => (typeof v === 'number' ? v : or);

/** The streams and their payloads' builders; the dispatcher asks `observed` before it builds. */
export class Interaction {
  readonly subjects: { readonly [K in SeqStreamName]: Subject<SeqStreams[K]> };
  readonly api: SeqInteraction;
  private readonly host: InteractionHost;

  constructor(host: InteractionHost) {
    this.host = host;
    const s = {} as Record<SeqStreamName, Subject<SeqInputEvent>>;
    for (const n of NAMES) s[n] = new Subject();
    this.subjects = s as unknown as Interaction['subjects'];
    const api = {
      on: (stream: SeqStreamName, scope?: SeqInputScope) => this.on(stream, scope),
      observed: (stream: SeqStreamName) => this.subjects[stream].observed,
      chrome: this.chrome,
    } as Record<string, unknown>;
    for (const n of NAMES) api[n] = s[n].asObservable();
    this.api = api as unknown as SeqInteraction;
  }

  private on(stream: SeqStreamName, scope?: SeqInputScope): Observable<SeqInputEvent> {
    const subject = this.subjects[stream] as Subject<SeqInputEvent>;
    if (!scope) return subject.asObservable();
    const { part, section, track } = scope;
    const kinds: readonly SeqPartKind[] | null = part === undefined ? null : typeof part === 'string' ? [part] : part;
    return subject.pipe(
      filter(
        (e) =>
          (kinds === null || kinds.includes(e.part.kind)) &&
          (section === undefined || e.part.section === section) &&
          (track === undefined || e.part.track === track),
      ),
    );
  }

  observed(stream: SeqStreamName): boolean {
    return this.subjects[stream].observed;
  }

  /** The hit's part: its region, and the section's, track's and group's refs. */
  private partOf(hit: Hit | null, section: number): SeqPart {
    const a = this.host.assembled();
    const rv = a?.view?.resolved;
    const str = (id: number) => (a ? a.set.strings.get(id) : null);
    const h = hit && hit.section === section ? hit : null;
    return {
      kind: REGIONS[h?.region ?? 0] ?? 'none',
      section: this.host.sectionRef(section),
      track: h && rv && h.track >= 0 ? str(rv.tracks.ref[h.track]) : null,
      group: h && rv && h.group >= 0 ? str(rv.groups.ref[h.group]) : null,
      key: null,
    };
  }

  private atOf(hit: Hit | null): SeqAt | null {
    if (!hit || hit.region === HitRegion.none) return null;
    const a = this.host.assembled();
    const rv = a?.view?.resolved;
    const ann = a && hit.annotation >= 0 && hit.item >= 0 ? a.set.annotations[hit.annotation] : undefined;
    const key = ann ? ann.itemKey[hit.item] : -1;
    return {
      layer: a && rv && hit.layer >= 0 ? a.set.strings.get(rv.layers.ref[hit.layer]) : null,
      lane: hit.lane,
      block: hit.block,
      index: hit.index,
      row: hit.displayRow,
      split: hit.split,
      item:
        a && ann
          ? { annotation: a.set.strings.get(ann.key), item: hit.item, id: key >= 0 ? a.set.strings.get(key) : null }
          : null,
    };
  }

  /**
   * The payload's common part. `x`, `y`: surface px; the client's from the native (a keyboard's: the surface's corner,
   * the native's current target, plus x and y).
   */
  private base(
    stream: SeqStreamName,
    origin: SeqInputEvent['origin'],
    part: SeqPart,
    at: SeqAt | null,
    x: number,
    y: number,
    buttons: number,
    button: number,
    modifiers: number,
    owner: string,
    native: SeqNativeEvent | null,
  ): Payload {
    const d = native as DomLike | null;
    let clientX = num(d?.clientX, Number.NaN);
    let clientY = num(d?.clientY, Number.NaN);
    if (origin === 'keyboard' || Number.isNaN(clientX)) {
      const el = (d?.currentTarget ?? d?.target) as { getBoundingClientRect?(): { left: number; top: number } } | null;
      const r = el?.getBoundingClientRect?.();
      clientX = r ? r.left + x : Number.NaN;
      clientY = r ? r.top + y : Number.NaN;
    }
    return {
      stream,
      origin,
      part,
      at,
      place: { x, y, clientX, clientY },
      buttons: buttons as ButtonsType,
      button: (button >= 0 ? (BUTTON[button] ?? 0) : 0) as ButtonsType.Flag,
      modifiers: MODIFIERS[modifiers & 15],
      time: num(d?.timeStamp, this.host.now()),
      owner,
      native: native ?? standIn(stream),
    };
  }

  /** Emits and lets go of the native: a kept payload holds the stand-in. Returns whether the default was prevented. */
  private send(stream: SeqStreamName, e: Payload): boolean {
    (this.subjects[stream] as Subject<SeqInputEvent>).next(e);
    const prevented = e.native.defaultPrevented;
    e.native = NO_NATIVE;
    return prevented;
  }

  /** hover, press, menu, activate from a pointer record and its hit (a keyboard's: its cell). */
  pointer(
    stream: SeqStreamName,
    origin: 'pointer' | 'keyboard',
    e: PointerInput | null,
    hit: Hit,
    owner: string,
  ): boolean {
    const section = e ? e.section : hit.section;
    const x = e ? e.x : hit.x;
    const y = e ? e.y : hit.y;
    const p = this.base(
      stream,
      origin,
      this.partOf(hit, section),
      this.atOf(hit.section === section ? hit : null),
      x,
      y,
      e?.buttons ?? 0,
      e?.button ?? -1,
      e?.modifiers ?? 0,
      owner,
      e?.native ?? null,
    );
    return this.send(stream, p);
  }

  /** activate from a key (the keyboard's Enter on its cursor's header). */
  keyed(stream: SeqStreamName, k: KeyInput, hit: Hit, owner: string): boolean {
    const p = this.base(
      stream,
      'keyboard',
      this.partOf(hit, hit.section),
      this.atOf(hit),
      hit.x,
      hit.y,
      0,
      -1,
      k.modifiers,
      owner,
      k.native,
    );
    return this.send(stream, p);
  }

  click(c: ClickEvt, origin: 'pointer' | 'keyboard', owner: string) {
    const e = c.input;
    const p = this.base(
      'click',
      origin,
      this.partOf(c.hit, e.section),
      this.atOf(c.hit),
      e.x,
      e.y,
      e.buttons,
      e.button,
      e.modifiers,
      owner,
      e.native,
    );
    p.count = c.count;
    this.send('click', p);
  }

  drag(d: DragEvt, section: number) {
    const e = d.input;
    const p = this.base(
      'drag',
      'pointer',
      this.partOf(d.hit, section),
      this.atOf(d.hit),
      e.x,
      e.y,
      e.buttons,
      d.button,
      e.modifiers,
      d.owner,
      e.native,
    );
    p.phase = PHASES[d.phase];
    p.dx = d.dx;
    p.dy = d.dy;
    p.totalX = d.totalX;
    p.totalY = d.totalY;
    p.start = this.partOf(d.start, d.start.section);
    this.send('drag', p);
  }

  wheel(w: WheelEvt) {
    const e = w.input;
    const p = this.base(
      'wheel',
      'pointer',
      this.partOf(w.hit, e.section),
      this.atOf(w.hit),
      e.x,
      e.y,
      0,
      -1,
      e.modifiers,
      '',
      e.native,
    );
    p.deltaX = e.deltaX;
    p.deltaY = e.deltaY;
    p.handled = w.handled;
    this.send('wheel', p);
  }

  key(k: KeyInput, handled: boolean) {
    const p = this.base(
      'key',
      'keyboard',
      this.partOf(null, k.section),
      null,
      Number.NaN,
      Number.NaN,
      0,
      -1,
      k.modifiers,
      '',
      k.native,
    );
    p.key = k.key;
    p.code = k.code;
    p.handled = handled;
    this.send('key', p);
  }

  leave(section: number, native: SeqNativeEvent | null) {
    const p = this.base(
      'leave',
      'pointer',
      this.partOf(null, section),
      null,
      Number.NaN,
      Number.NaN,
      0,
      -1,
      0,
      '',
      native,
    );
    this.send('leave', p);
  }

  /** seqspec:interaction#streams: a UI's chrome, through the same streams. */
  readonly chrome: SeqInteraction['chrome'] = (stream, part, native, place, opts) => {
    if (!this.subjects[stream].observed) return;
    const d = native as DomLike;
    const p = this.base(
      stream,
      'chrome',
      part,
      null,
      place.x,
      place.y,
      opts?.buttons ?? num(d.buttons, 0),
      opts?.button ?? -1,
      opts?.modifiers ?? 0,
      '',
      native,
    );
    p.place = place;
    if (stream === 'click') p.count = opts?.count ?? 1;
    this.send(stream, p);
  };

  dispose() {
    for (const n of NAMES) this.subjects[n].complete();
  }
}
