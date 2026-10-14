/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Input and base events (seqspec:interaction#input-events). The engine (or a test) fills reused, DOM-free input records
// and hands them to `ctx.input`; the dispatcher hit-tests them against the headless layout, writes `ctx.pointer` and
// emits base events, as Mol*'s interaction behaviors. It decides nothing else: who owns a gesture is the subscribers'
// claims, in priority order, and what a gesture does is behaviors'. Every record is reused: a move allocates nothing,
// and the element under the pointer is compared by `HitKey`, straight-line (IN2's hot-path design).

import { type ConfigManager, SeqViewerConfig } from "./config.js";
import { Hit, HitKey, HitRegion } from "./hit.js";
import type { HoverManager, PointerState } from "./hover.js";
import type { LayoutManager } from "./layout.js";
import type { Scheduler } from "./scheduler.js";
import { event, type SeqEvent, type SeqEventEmitter, type WritableLive } from "./state.js";
import type { DragState } from "./ui.js";

export const PointerType = { down: 0, move: 1, up: 2, cancel: 3 } as const;
/** Bits of `modifiers`. */
export const Modifier = { shift: 1, ctrl: 2, alt: 4, meta: 8 } as const;
/** `PointerInput.button`, as the DOM numbers it; -1 on a move. */
export const Button = { primary: 0, middle: 1, secondary: 2 } as const;
export const DragPhase = { start: 0, move: 1, end: 2, cancel: 3 } as const;

/** One record per engine, reused. x and y are CSS px in the section's surface (outside it during a captured drag). */
export interface PointerInput {
  /** `PointerType`. */
  type: number;
  section: number;
  x: number;
  y: number;
  buttons: number;
  button: number;
  /** `Modifier` bits. */
  modifiers: number;
  /** A press's or a release's time (ms); a move's is left as it was (the engine reads no time per move). */
  time: number;
  pointerId: number;
  /**
   * During a press, a point outside the pressing surface is clamped to its bounds (x and y above); these are how far
   * past them it was (seqspec:interaction#input-events), 0 inside: what a scroll past an edge reads.
   */
  overX: number;
  overY: number;
}
/** Deltas in CSS px (the engine converts line and page modes). */
export interface WheelInput {
  section: number;
  x: number;
  y: number;
  deltaX: number;
  deltaY: number;
  modifiers: number;
  time: number;
}
export interface KeyInput {
  /** The focused surface's section. */
  section: number;
  /** `KeyboardEvent.key`. */
  key: string;
  /** `KeyboardEvent.code`: the physical key, for the keys a modifier changes (Alt+= is not `=` on macOS). */
  code: string;
  modifiers: number;
  /** The event's target is an editable element: the key is the page's, never the viewer's. */
  typing: boolean;
  time: number;
}
/** The records' doubles are created as doubles (0.5): a Smi field turning Double deprecates the record's map. */
export const Inputs = {
  pointer: (): PointerInput => ({
    type: 1,
    section: -1,
    x: 0.5,
    y: 0.5,
    buttons: 0,
    button: -1,
    modifiers: 0,
    time: 0.5,
    pointerId: 0,
    overX: 0.5,
    overY: 0.5,
  }),
  wheel: (): WheelInput => ({ section: -1, x: 0.5, y: 0.5, deltaX: 0.5, deltaY: 0.5, modifiers: 0, time: 0.5 }),
  key: (): KeyInput => ({ section: -1, key: "", code: "", modifiers: 0, typing: false, time: 0.5 }),
};

export interface InputSink {
  pointer(e: PointerInput): void;
  /** True: handled; the engine prevents the page's default. */
  wheel(e: WheelInput): boolean;
  key(e: KeyInput): boolean;
  /** A surface gained or lost keyboard focus. */
  focus(section: number, focused: boolean): void;
  leave(section: number): void;
  /** A press or drag in progress ends as a cancel; `section`: only if it started on that surface. */
  cancel(section?: number): void;
  /** The engine is about to draw: a re-hit waiting for the frame (`InputDispatcher.moved`) runs now, drawn with it. */
  flush?(): void;
}

/** The first subscriber to `claim` owns the gesture; later ones see `owner` set. */
export class PointerEvt {
  input: PointerInput = Inputs.pointer();
  hit: Hit = Hit.create();
  owner = "";
  claim(name: string) {
    if (this.owner === "") this.owner = name;
  }
}
export interface ClickEvt {
  readonly input: PointerInput;
  readonly hit: Hit;
  /** 2: a double click, 3: a triple (then the count starts again). */
  count: number;
}
export class DragEvt extends PointerEvt {
  /** `DragPhase`. */
  phase = 0;
  /** The hit at the press. */
  start: Hit = Hit.create();
  /** The press's button (moves carry none). */
  button = 0;
  /** Since the last drag event, and since the press (doubles from birth). */
  dx = 0.5;
  dy = 0.5;
  totalX = 0.5;
  totalY = 0.5;
}
export interface WheelEvt {
  readonly input: WheelInput;
  readonly hit: Hit;
  handled: boolean;
}
export interface KeyEvt {
  readonly input: KeyInput;
  handled: boolean;
}
export interface InteractionEvents {
  /** Every pointer input, with its hit. */
  readonly pointer: SeqEvent<PointerEvt>;
  /** The element under the pointer changed. */
  readonly hover: SeqEvent<PointerEvt>;
  /** A press and release within the drag slop. */
  readonly click: SeqEvent<ClickEvt>;
  /** start, move, end, cancel; move and end only for a claimed drag. */
  readonly drag: SeqEvent<DragEvt>;
  /**
   * A query, asked by a click's subscriber once per click (`InputDispatcher.askDoubleClick`), never per move: a
   * subscriber whose double click means something at the hit (a zoom, an activation) claims it.
   */
  readonly doubleClick: SeqEvent<PointerEvt>;
  readonly wheel: SeqEvent<WheelEvt>;
  readonly key: SeqEvent<KeyEvt>;
  readonly focus: SeqEvent<{ section: number; focused: boolean }>;
  readonly leave: SeqEvent<{ section: number }>;
}

export interface Trigger {
  readonly button?: "primary" | "secondary" | "middle";
  readonly shift?: boolean;
  readonly ctrl?: boolean;
  readonly alt?: boolean;
  readonly meta?: boolean;
}
export type Binding = readonly Trigger[];
const bits = (t: Trigger) =>
  (t.shift ? Modifier.shift : 0) |
  (t.ctrl ? Modifier.ctrl : 0) |
  (t.alt ? Modifier.alt : 0) |
  (t.meta ? Modifier.meta : 0);
export const Binding = {
  /** Exact modifiers; a trigger with no button matches any button. Allocation-free. */
  matches(t: Trigger, button: number, modifiers: number): boolean {
    return bits(t) === modifiers && (t.button === undefined || Button[t.button] === button);
  },
  match(b: Binding, button: number, modifiers: number): boolean {
    for (let i = 0; i < b.length; i++) if (Binding.matches(b[i], button, modifiers)) return true;
    return false;
  },
};

export interface KeyTrigger {
  readonly key?: string;
  readonly code?: string;
  readonly shift?: boolean;
  readonly ctrl?: boolean;
  readonly alt?: boolean;
  readonly meta?: boolean;
}
export type KeyBinding = readonly KeyTrigger[];
export const KeyBinding = {
  /**
   * Exact Ctrl, Meta and Alt (Ctrl/Cmd with `+` stays the browser's zoom); Shift only for a named key
   * (`key.length > 1`: Enter, F10), since a printable key already has it in its value. `code` matches the physical
   * key, `key` the value. Allocation-free.
   */
  matches(t: KeyTrigger, k: KeyInput): boolean {
    const mask = k.key.length > 1 ? -1 : ~Modifier.shift;
    if ((bits(t) & mask) !== (k.modifiers & mask)) return false;
    return t.code !== undefined ? t.code === k.code && (t.key === undefined || t.key === k.key) : t.key === k.key;
  },
  match(b: KeyBinding, k: KeyInput): boolean {
    for (let i = 0; i < b.length; i++) if (KeyBinding.matches(b[i], k)) return true;
    return false;
  },
};

export interface InputHost {
  readonly layout: LayoutManager;
  readonly pointer: WritableLive<PointerState>;
  readonly hover: HoverManager;
  readonly config: ConfigManager;
  readonly scheduler: Scheduler;
  /** False while the viewer is paused or disposed: no events. */
  active(): boolean;
  setDrag(drag: DragState | null): void;
  onError(e: unknown): void;
}

export class InputDispatcher implements InputSink {
  readonly events: { readonly [K in keyof InteractionEvents]: SeqEventEmitter<EventOf<K>> };
  private readonly host: InputHost;
  private readonly pe = new PointerEvt();
  private readonly de = new DragEvt();
  private readonly ce: ClickEvt = { input: Inputs.pointer(), hit: Hit.create(), count: 1 };
  /** The double-click query's record, reused. */
  private readonly dq = new PointerEvt();
  /** In a release the engine never delivered (`unseenRelease`), and in a keyboard's click (`synthClick`). */
  private unseen = false;
  private synthetic = false;
  private readonly we: WheelEvt = { input: Inputs.wheel(), hit: Hit.create(), handled: false };
  private readonly ke: KeyEvt = { input: Inputs.key(), handled: false };
  private readonly fe = { section: -1, focused: false };
  private readonly le = { section: -1 };
  private readonly hit = Hit.create();
  /** The last hovered element's `HitKey`; -2: none since a leave or a cancel. */
  private readonly last = new Int32Array(HitKey.SIZE).fill(-2);
  /**
   * The press (x, y), the last drag event's, the last click's time, x and y, and the press's time: Float64Array so
   * nothing boxes. The last click's section and button beside it.
   */
  private readonly at = new Float64Array(8);
  /** The last click's section, button and count. */
  private readonly clicked = new Int32Array(3).fill(-1);
  private pressed = false;
  private dragging = false;
  private owner = "";
  /** The pointer that pressed, and the surface it pressed on: the drag follows both. */
  private pressId = 0;
  private pressSection = -1;
  /**
   * The press was claimed on down (its drag starts on any move) and has stayed within the slop: released soon
   * enough, it is also a click (a double click wobbles).
   */
  private still = false;
  /** Dispatches in progress, and a section rebuilt during one (-1: none): `rehit` waits for their end. */
  private depth = 0;
  private stale = -1;
  /** The section the dispatcher is hit testing now (-1: none): a rebuild inside its own hit is not re-hit. */
  private hitting = -1;
  /** A moved section's re-hit, waiting for the frame (`moved`; -1: none), and that frame's id. */
  private moving = -1;
  private frame = -1;
  /** In a moved section's re-hit: a move its hover's subscribers make is not hit again (no loop). */
  private following = false;

  constructor(host: InputHost) {
    this.host = host;
    const e = <T>() => event<T>(host.onError);
    this.events = {
      pointer: e(),
      hover: e(),
      click: e(),
      drag: e(),
      doubleClick: e(),
      wheel: e(),
      key: e(),
      focus: e(),
      leave: e(),
    };
    this.at[4] = Number.NEGATIVE_INFINITY;
  }

  pointer(e: PointerInput) {
    if (!this.host.active()) return;
    // One pointer at a time; a move with no button down while pressed, or a down again (from any pointer: a touch whose
    // release the page lost leaves the next touch, with a new id, ending it), follows a release we never saw: the
    // press ends first, released at the last point it saw, so its gesture completes. Another pointer's moves and ups
    // are not the press's.
    if (this.pressed) {
      const own = e.pointerId === this.pressId;
      if (!own && e.type !== PointerType.down) return;
      if (e.type === PointerType.down || (e.type === PointerType.move && e.buttons === 0)) this.unseenRelease(e);
    }
    this.depth++;
    try {
      this.dispatch(e);
    } finally {
      this.depth--;
    }
    if (this.depth === 0 && this.stale >= 0) this.rehit(this.stale);
  }

  /**
   * A release the engine never delivered: dispatched as a release at the last point the press saw (its last drag
   * event's, else the press's), on the record itself, which is then restored; it ends a drag and is never a click.
   * Edge-rate.
   */
  private unseenRelease(e: PointerInput) {
    const type = e.type;
    const x = e.x;
    const y = e.y;
    const button = e.button;
    const buttons = e.buttons;
    e.type = PointerType.up;
    e.x = this.at[2];
    e.y = this.at[3];
    e.button = this.de.button;
    e.buttons = 0;
    this.depth++;
    this.unseen = true;
    try {
      this.dispatch(e);
    } finally {
      this.unseen = false;
      this.depth--;
      e.type = type;
      e.x = x;
      e.y = y;
      e.button = button;
      e.buttons = buttons;
    }
  }

  /**
   * The section's layout was rebuilt (`LayoutHost.rebuilt`): what is under a still pointer is hit again, as a move, so
   * the hover and its label follow without the pointer moving. Inside the dispatcher's own hit it is ignored (that hit
   * reads the new layout); during a dispatch (a click's handler rebuilt it) it waits for the end. Edge-rate. The same
   * element again is refreshed (`HoverManager.refresh`), unless `refresh` is false (a view that only moved).
   */
  rehit(section: number, refresh = true) {
    if (section === this.hitting) return;
    if (this.depth > 0) {
      // Several rebuilt during one dispatch: the pointer's section is the one to hit again.
      if (this.stale < 0 || section === this.host.pointer.current.section) this.stale = section;
      return;
    }
    this.stale = -1;
    const p = this.host.pointer.current;
    // Not during a press, and not over the keyboard's hover, which stays as it is.
    if (this.pressed || !p.inside || p.section !== section || !this.host.active()) return;
    if (this.host.hover.live.current.owner === "keyboard") return;
    const h = this.hit;
    h.x = p.x;
    h.y = p.y;
    this.hitting = section;
    p.inside = this.host.layout.hitTestAt(section, h);
    this.hitting = -1;
    this.pe.hit = h;
    this.pe.owner = "";
    if (!HitKey.same(h, this.last)) {
      HitKey.write(h, this.last);
      this.events.hover.emit(this.pe);
    } else if (refresh) this.host.hover.refresh();
  }

  /**
   * Section `section`'s view moved (a wheel, a zoom, a key, a programmatic or `follow` move, an animation's frame) or
   * its content changed (regions arrived); -1: the pointer's section. Under a still pointer, the content moved: it is
   * hit again, once, at the next frame (or the engine's draw before it, `flush`), and the hover follows it. No frame
   * while the pointer is elsewhere, pressed, or the keyboard holds the hover. Allocation-free; frame-rate.
   */
  moved(section: number) {
    const p = this.host.pointer.current;
    if (section < 0) section = p.section;
    if (this.following || this.pressed || !p.inside || p.section !== section || !this.host.active()) return;
    if (this.host.hover.live.current.owner === "keyboard") return;
    this.moving = section;
    if (this.frame < 0) this.frame = this.host.scheduler.requestFrame(this.onFrame);
  }
  /**
   * The moved section's re-hit now, if one waits (the engine draws next: the hover is drawn with the view). Its frame
   * stays requested and then finds nothing: a scroll's frames request no more of them, and cancel none.
   */
  flush() {
    const s = this.moving;
    if (s < 0) return;
    this.moving = -1;
    this.following = true;
    try {
      this.rehit(s, false);
    } finally {
      this.following = false;
    }
  }
  private readonly onFrame = () => {
    this.frame = -1;
    this.flush();
  };

  private dispatch(e: PointerInput) {
    const t = e.type;
    const h = this.hit;
    // A press keeps its surface wherever the pointer goes, its point clamped to the surface's bounds (how far past
    // them kept in `overX`, `overY`).
    const section = this.pressed ? this.pressSection : e.section;
    this.clamp(e, section, this.pressed && t !== PointerType.down);
    h.x = e.x;
    h.y = e.y;
    this.hitting = section;
    const inside = this.host.layout.hitTestAt(section, h);
    this.hitting = -1;
    const p = this.host.pointer.current;
    p.section = section;
    p.x = e.x;
    p.y = e.y;
    p.inside = inside;
    this.host.pointer.bump();
    const pe = this.pe;
    pe.input = e;
    pe.hit = h;
    pe.owner = t === PointerType.down ? "" : this.owner;
    if (t === PointerType.down) {
      const s = this.de.start;
      s.x = e.x;
      s.y = e.y;
      this.host.layout.hitTestAt(section, s);
      this.pressed = true;
      this.dragging = false;
      this.pressId = e.pointerId;
      this.pressSection = section;
    }
    if (this.pressed) this.de.input = e;
    this.events.pointer.emit(pe);
    if (!HitKey.same(h, this.last)) {
      HitKey.write(h, this.last);
      this.events.hover.emit(pe);
    }
    const at = this.at;
    if (t === PointerType.down) {
      this.owner = pe.owner;
      this.still = pe.owner !== "";
      this.de.button = e.button;
      at[0] = at[2] = e.x;
      at[1] = at[3] = e.y;
      at[7] = this.host.scheduler.now();
    } else if (t === PointerType.move && this.pressed) {
      if (this.dragging && !this.still) this.drag(e, DragPhase.move);
      else {
        const slop = this.host.config.get(SeqViewerConfig.Input.DragSlopPx);
        const dx = e.x - at[0];
        const dy = e.y - at[1];
        const out = dx * dx + dy * dy > slop * slop;
        if (out) this.still = false;
        if (this.dragging) this.drag(e, DragPhase.move);
        else if (this.owner !== "" || out) this.drag(e, DragPhase.start);
      }
    } else if (t === PointerType.up && this.pressed) {
      const slop = this.host.config.get(SeqViewerConfig.Input.DragSlopPx);
      const dx = e.x - at[0];
      const dy = e.y - at[1];
      // A claimed press's drag ends first; within the slop and the click time, its release is a click too. A release
      // never seen is no click (its time and place are unknown).
      const click =
        !this.unseen &&
        (!this.dragging ||
          (this.still &&
            dx * dx + dy * dy <= slop * slop &&
            this.host.scheduler.now() - at[7] <= this.host.config.get(SeqViewerConfig.Input.DoubleClickMs)));
      if (this.dragging) this.drag(e, DragPhase.end);
      if (click) this.click(e, section);
      this.pressed = this.dragging = this.still = false;
      this.owner = "";
    } else if (t === PointerType.cancel) this.cancel();
  }

  /** A pressed point into section `section`'s bounds (`clamp`), else its overshoot 0. Allocation-free. */
  private clamp(e: PointerInput, section: number, clamp: boolean) {
    e.overX = 0;
    e.overY = 0;
    const ref = clamp ? this.host.layout.sectionRef(section) : null;
    const l = ref === null ? null : this.host.layout.get(ref)?.current;
    if (!l) return;
    const x = e.x < 0 ? 0 : e.x > l.width ? l.width : e.x;
    const y = e.y < 0 ? 0 : e.y > l.height ? l.height : e.y;
    e.overX = e.x - x;
    e.overY = e.y - y;
    e.x = x;
    e.y = y;
  }

  /**
   * The press ends without a release: the viewer paused (`setActive(false)`), the surface it started on went away
   * (`section`), the engine cancelled it, or a move came with no button down. An owned drag gets phase cancel and
   * `ui.drag` clears; the next move is a hover again.
   */
  cancel(section?: number) {
    if (!this.pressed || (section !== undefined && section !== this.pressSection)) return;
    const d = this.de;
    if (this.dragging && this.owner !== "") {
      d.phase = DragPhase.cancel;
      d.hit = this.hit;
      d.owner = this.owner;
      d.dx = d.dy = d.totalX = d.totalY = 0.5;
      this.events.drag.emit(d);
      this.host.setDrag(null);
    }
    this.pressed = this.dragging = this.still = false;
    this.owner = "";
    this.last.fill(-2);
  }

  private drag(e: PointerInput, phase: number) {
    const d = this.de;
    const at = this.at;
    d.phase = phase;
    d.input = e;
    d.hit = this.hit;
    d.dx = e.x - at[2];
    d.dy = e.y - at[3];
    d.totalX = e.x - at[0];
    d.totalY = e.y - at[1];
    at[2] = e.x;
    at[3] = e.y;
    d.owner = this.owner;
    if (phase === DragPhase.start) this.dragging = true;
    // An unclaimed drag does nothing: its moves and end are not emitted.
    if (phase === DragPhase.start || d.owner !== "") this.events.drag.emit(d);
    if (phase === DragPhase.start && d.owner !== "") {
      this.owner = d.owner;
      this.host.setDrag({ section: d.start.section, owner: d.owner, region: d.start.region });
    } else if (phase >= DragPhase.end && d.owner !== "") this.host.setDrag(null);
  }

  private click(e: PointerInput, section: number) {
    const at = this.at;
    const now = this.host.scheduler.now();
    const slop = this.host.config.get(SeqViewerConfig.Input.DragSlopPx);
    const button = this.de.button;
    const near = Math.abs(e.x - at[5]) <= slop && Math.abs(e.y - at[6]) <= slop;
    // The same button on the same surface, each click within the double-click time of the one before.
    const again =
      near &&
      this.clicked[0] === section &&
      this.clicked[1] === button &&
      now - at[4] <= this.host.config.get(SeqViewerConfig.Input.DoubleClickMs);
    const c = this.ce as { input: PointerInput; hit: Hit; count: number };
    c.input = e;
    c.hit = this.hit;
    c.count = again ? this.clicked[2] + 1 : 1;
    this.clicked[2] = c.count;
    // A triple click starts the count again.
    at[4] = c.count >= 3 ? Number.NEGATIVE_INFINITY : now;
    at[5] = e.x;
    at[6] = e.y;
    this.clicked[0] = section;
    this.clicked[1] = button;
    this.events.click.emit(c);
  }

  wheel(e: WheelInput): boolean {
    if (!this.host.active()) return false;
    const w = this.we as { input: WheelInput; hit: Hit; handled: boolean };
    this.depth++;
    try {
      w.hit.x = e.x;
      w.hit.y = e.y;
      this.host.layout.hitTestAt(e.section, w.hit);
      w.input = e;
      w.handled = false;
      this.events.wheel.emit(w);
    } finally {
      this.depth--;
    }
    if (this.depth === 0 && this.stale >= 0) this.rehit(this.stale);
    return w.handled;
  }

  key(e: KeyInput): boolean {
    if (!this.host.active() || e.typing) return false;
    const k = this.ke as { input: KeyInput; handled: boolean };
    k.input = e;
    k.handled = false;
    this.events.key.emit(k);
    return k.handled;
  }

  /**
   * Whether a double click at click `c`'s hit means something to a behavior (`doubleClick`): its name, or "" for
   * none, and "" for a keyboard's click (it has no second). Edge-rate: once per click.
   */
  askDoubleClick(c: ClickEvt): string {
    if (this.synthetic || !this.host.active()) return "";
    const q = this.dq;
    q.input = c.input;
    q.hit = c.hit;
    q.owner = "";
    this.events.doubleClick.emit(q);
    return q.owner;
  }

  /**
   * A click with no press: the keyboard's Enter on its cursor's cell (seqspec:viewer#behaviors); never half of a double
   * click.
   */
  synthClick(c: ClickEvt) {
    if (!this.host.active()) return;
    this.synthetic = true;
    try {
      this.events.click.emit(c);
    } finally {
      this.synthetic = false;
    }
  }

  focus(section: number, focused: boolean) {
    if (!this.host.active()) return;
    this.fe.section = section;
    this.fe.focused = focused;
    this.events.focus.emit(this.fe);
  }

  leave(section: number) {
    if (!this.host.active()) return;
    const p = this.host.pointer;
    p.current.section = -1;
    p.current.inside = false;
    p.bump();
    this.last.fill(-2);
    this.hit.region = HitRegion.none;
    this.le.section = section;
    this.events.leave.emit(this.le);
  }

  dispose() {
    for (const e of Object.values(this.events)) e.complete();
  }
}

type EventOf<K extends keyof InteractionEvents> = InteractionEvents[K] extends SeqEvent<infer T> ? T : never;
