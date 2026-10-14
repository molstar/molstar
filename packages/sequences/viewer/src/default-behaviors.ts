/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The default behaviors that follow state and commands (seqspec:viewer#behaviors): hover, its label and the cursor, the
// focus and scroll-group followers and the player's stop on input (`commands` is CM1's). The gestures are in
// gestures.ts. Each reads its params through `forSection` when it acts, so a `setParams` needs no re-creation, and
// each undoes what it set when removed (`add`'s teardowns).

import type { Closure } from "@molstar/seq-core/format";
import type { TeardownLogic } from "rxjs";
import { type BehaviorHandle, type SeqBehaviorProvider, SeqBehaviors } from "./behaviors.js";
import type { SeqViewerContext } from "./context.js";
import { DefaultHeaderHoverParams, type HeaderHoverParams } from "./header-hover.js";
import { type Hit, HitRegion, HoverKind } from "./hit.js";
import { DragPhase, type InputDispatcher, PointerType } from "./input.js";
import { RING_WINDOW, ringWindowAt } from "./ring-window.js";
import type { Origin, SeqEvent, Unsubscribable } from "./state.js";
import type { RowHoverParams } from "./ui.js";
import type { Viewport } from "./viewport.js";
import { WINDOW_FRAME, WindowKind, WindowPart, windowFrame, windowIdle, windowKind, windowPart } from "./windows.js";

export interface BehaviorKit<P> {
  readonly ctx: SeqViewerContext;
  readonly h: BehaviorHandle<P>;
  /** Subscribes at the behavior's priority; released on dispose. */
  on<T>(ev: SeqEvent<T>, fn: (e: T) => void): void;
  /** A subscription, or a teardown run on dispose (what the behavior set: its hover, its cursors, its marks). */
  add(t: TeardownLogic): void;
}

/** A behavior from a setup run at `register`; everything it adds is released on dispose. */
export function defineBehavior<P extends object>(
  name: string,
  category: SeqBehaviorProvider["category"],
  priority: number,
  params: P,
  description: string,
  setup: (k: BehaviorKit<P>) => void,
): SeqBehaviorProvider<P> {
  const display = {
    name: name[0].toUpperCase() + name.slice(1).replace(/[A-Z]/g, (c) => ` ${c.toLowerCase()}`),
    description,
  };
  return {
    name,
    display,
    category,
    priority,
    params,
    create(ctx, h) {
      const subs = SeqBehaviors.handler();
      const kit: BehaviorKit<P> = {
        ctx,
        h,
        on: (ev, fn) => void subs.add(ev.subscribe(fn, { priority })),
        add: (t) => void subs.add(t),
      };
      return { register: () => setup(kit), update: () => true, dispose: () => subs.unsubscribe() };
    },
  };
}

/** `HitRegion` names by value: params name regions. */
export const RegionName = Object.keys(HitRegion) as readonly (keyof typeof HitRegion)[];
export type RegionName = keyof typeof HitRegion;
export const viewportOf = (ctx: SeqViewerContext, section: number) =>
  ctx.viewports.get(ctx.layout.sectionRef(section) ?? "") as Viewport | null;
/** The context's diagnostics and its dispatcher, not on its public interface. */
export const internal = (ctx: SeqViewerContext) =>
  ctx as unknown as {
    log(level: "info" | "warning" | "error", message: string, data?: unknown): void;
    report(code: string, message: string): void;
    activate(section: string, track: string, origin: Origin): void;
    openMenu(section: string, track: string, origin: Origin): void;
    readonly input: InputDispatcher;
  };

const ALL_BUT_SCROLLBAR: RegionName[] = ["data", "column", "header", "gutter", "dataColumn", "dataColumnHeader"];

export const HoverBehavior = defineBehavior(
  "hover",
  "interaction",
  0,
  {
    regions: [...ALL_BUT_SCROLLBAR, "navigator", "vnav"] as readonly RegionName[],
    whileDragging: "clear" as "clear" | "keep",
    /** `"element"`: a nucleotide lights alone (seqspec:interaction#hover). */
    closure: "containing" as Closure,
  },
  "The element under the pointer is the hover; cleared on leave and during a drag another behavior owns.",
  ({ ctx, h, on, add }) => {
    let dragging = false;
    const set = (hit: Hit) => {
      const p = h.forSection(hit.section);
      if (p?.regions.includes(RegionName[hit.region])) ctx.hover.setFromHit(hit, "pointer", p.closure);
      else ctx.hover.clear("pointer");
    };
    on(ctx.interaction.hover, (e) => {
      if (!dragging) set(e.hit);
    });
    on(ctx.interaction.leave, () => ctx.hover.clear("pointer"));
    on(ctx.interaction.drag, (e) => {
      if (e.phase === DragPhase.start && e.owner !== "" && h.forSection(e.start.section)?.whileDragging === "clear") {
        dragging = true;
        ctx.hover.clear("pointer");
      } else if (e.phase >= DragPhase.end && dragging) {
        dragging = false;
        set(e.hit);
      }
    });
    add(() => ctx.hover.clear("pointer"));
  },
);

/**
 * What a header's hover stands for (seqspec:interaction#hover): a feature layer's shown items (a split row's, a rows
 * track's row's, a group's tracks'), a continuous layer's row band alone, the row band alone past `maxItems`; and
 * whether linked views get the set. Removed or replaced in the viewer spec, a header's hover is the header alone.
 */
export const HeaderHoverBehavior = defineBehavior(
  "headerHover",
  "interaction",
  0,
  { ...DefaultHeaderHoverParams } as HeaderHoverParams,
  "A header's hover lights what it labels: a feature track's shown items (a split row's, a group's tracks'), a continuous track's row band; the row band alone past maxItems; forward gives linked views the set.",
  ({ ctx, h, add }) => {
    ctx.hover.setHeaderHover((section) => h.forSection(section));
    add(() => ctx.hover.setHeaderHover(null));
  },
);

/**
 * A data hover's row band (seqspec:interaction#hover): where the pointer is, drawn on the overlay (`ui.rowHover`). It
 * never changes the hover (its kind, its loci, what linked views get), so it reads as no header's hover; and it is
 * drawn at half a header hover's strength, so a header's still reads as the whole track (David, 2026-10-13).
 */
export const HoverRowBehavior = defineBehavior(
  "hoverRow",
  "interaction",
  0,
  { band: "entry", data: 0.25, header: 0.5 } as RowHoverParams,
  "A hover on a track's data (linear) lights its row band: the track's box, a split row's, a rows track's row, or its " +
    "renderer's own (a matrix's hovered band); none on the ruler, a dot plot and navigators. A drawing only: the " +
    "hover, its loci and what linked views get stay the element's.",
  ({ ctx, h, add }) => {
    const set = () => ctx.ui.setRowHover((section) => h.forSection(section));
    set();
    // A new function per params change: the overlay redraws.
    add(h.params.subscribe(set));
    add(() => ctx.ui.setRowHover(null));
  },
);

export const HoverLabelBehavior = defineBehavior(
  "hoverLabel",
  "interaction",
  0,
  {
    regions: ["data", "column", "header", "dataColumn"] as readonly RegionName[],
    whileDragging: "hide" as "hide" | "show",
  },
  "Computes the hovered element's label (labels.current); none over chrome and while dragging.",
  ({ ctx, h, add }) => {
    const update = () => {
      const s = ctx.hover.live.current;
      const p = h.forSection(s.section);
      const show =
        p !== null &&
        !(ctx.ui.drag.value !== null && p.whileDragging === "hide") &&
        s.kind !== HoverKind.none &&
        p.regions.includes(RegionName[s.region]);
      ctx.labels.current.set(show ? ctx.labels.compute(s) : null);
    };
    // The hovered element, a drag, and what a label is made of: the providers and the view (a numbering
    // switch is a new view). Never a move: `compute` caches on all three.
    add(ctx.hover.live.version.subscribe(update));
    add(ctx.ui.drag.subscribe(update));
    add(ctx.labels.providers.subscribe(update));
    add(ctx.graph.assembled.subscribe(update));
    // A region arriving (`compute` is made again only for a label that waits on one: a streamed row's tile).
    add(ctx.streams.version.subscribe(update));
    add(() => ctx.labels.current.set(null));
  },
);

// Read once: a member read per pointer event allocated under the test runner's module wrapper.
const RING_KIND = WindowKind.ring;
const NAVIGATOR_KIND = WindowKind.navigator;
const PART_NONE = WindowPart.none;
/** The cursor for a window part (windows.ts): the gesture a press there starts. */
const WINDOW_CURSORS = ["", "grab", "ew-resize", "ew-resize", "ns-resize", "ns-resize", "pointer"];

export const CursorBehavior = defineBehavior(
  "cursor",
  "interaction",
  0,
  {
    byRegion: {
      data: "grab",
      column: "crosshair",
      header: "pointer",
      navigator: "pointer",
      vnav: "pointer",
    } as Readonly<Partial<Record<RegionName, string>>>,
    /** A rows layer's header off its tree's nodes (a press there toggles nothing); over a node, `byRegion.header`. */
    rowHeader: "default",
    /** A dot plot's data, where a drag selects a rectangle (`dotPlotSelect`). */
    dotPlot: "crosshair",
    dragging: "grabbing",
  },
  "The surface's cursor by region (a window's body, edges and outside on the navigator, overview and scrollbar), and " +
    "for a whole drag the press's (`grabbing` where it was `grab`); `grab` on data only where dragPan is on; on a " +
    "rows layer's header, `rowHeader` off its tree's nodes; on a dot plot's data, `dotPlot`. With the cursor, the " +
    "window part lit (`ui.window`: a grip, the body), the press's for a whole window drag.",
  ({ ctx, h, on, add }) => {
    let section = -1;
    let region = 0;
    /** The window's kind and part under the point (`partAt`), and the cursor and window part at the press, kept for
     * the drag. */
    let kind = -1;
    let part: number = WindowPart.none;
    /** On a rows layer's header row, not over a tree node. */
    let offNode = false;
    /** On a dot plot's data (`dotPlot`). */
    let plot = false;
    let pressed = "default";
    let heldKind = -1;
    let heldPart: number = WindowPart.none;
    /** The window frame of a section's lookup, reused (windows.ts: lookups, no allocation). */
    const g = new Float64Array(WINDOW_FRAME);
    /** The window kind `partAt` found: its region's, or on rings the ring's window (ring-window.ts); -1: none. */
    let found = -1;
    /** The navigator under the point has nowhere to go (windows.ts `windowIdle`): inert, so no part, no `pointer`. */
    let idle = false;
    const rg = new Float64Array(RING_WINDOW);
    /** Over the navigator, the vertical overview, the scrollbar or a ring's window: its part; `none` elsewhere. */
    const partAt = (sec: number, reg: number, at: Hit): number => {
      const k = windowKind(reg);
      found = k;
      idle = false;
      if (k < 0) {
        const part = ringWindowAt(ctx, sec, at, rg);
        found = part === PART_NONE ? -1 : RING_KIND;
        return part;
      }
      const ref = ctx.layout.sectionRef(sec);
      const v = ref === null ? null : ctx.layout.axis(ref);
      const l = ref === null ? null : ctx.layout.get(ref)?.current;
      if (ref === null || !v || !l) return WindowPart.none;
      if (k === NAVIGATOR_KIND) {
        const vp = viewportOf(ctx, sec);
        idle = vp !== null && windowIdle(vp.live.current, vp.axis);
        if (idle) return WindowPart.none;
      }
      windowFrame(v, l, k, g);
      return windowPart(k, at, g, ctx.layout.props(ref).verticalZoom);
    };
    const apply = () => {
      const ref = ctx.layout.sectionRef(section);
      const p = h.forSection(section);
      const d = ctx.ui.drag.value;
      // A window's drag keeps its press's part lit (a press outside the navigator's or overview's window jumps there and
      // pans it), wherever the pointer goes; else the part under the pointer.
      const held =
        d !== null &&
        d.section === section &&
        heldKind >= 0 &&
        (heldKind === RING_KIND ? d.owner === "ringWindow" : windowKind(d.region) === heldKind);
      ctx.ui.setWindow(section, held ? heldKind : kind, held ? heldPart : part);
      if (ref === null || !p) return;
      if (d && d.section === section) {
        // The whole drag keeps the press's cursor, wherever the pointer goes (a move never resets it).
        const pan = d.owner === "dragPan" || pressed === "grab";
        return ctx.ui.setCursor(ref, pan ? p.dragging : pressed);
      }
      let c =
        (idle ? "default" : null) ??
        (part !== WindowPart.none ? WINDOW_CURSORS[part] : null) ??
        (offNode && region === HitRegion.header
          ? p.rowHeader
          : plot && region === HitRegion.data
            ? p.dotPlot
            : p.byRegion[RegionName[region]]) ??
        "default";
      if (
        c === "grab" &&
        part === WindowPart.none &&
        !ctx.behaviors
          .get<{ dragRegions: readonly string[] }>("dragPan")
          ?.forSection(section)
          ?.dragRegions.includes(RegionName[region])
      )
        c = "default";
      ctx.ui.setCursor(ref, c);
    };
    // Followed during a drag too (apply keeps the drag's cursor), so its end shows what is under the pointer.
    on(ctx.interaction.hover, (e) => {
      section = e.hit.section;
      region = e.hit.region;
      part = partAt(section, region, e.hit);
      kind = found;
      offNode = e.hit.displayRow >= 0 && e.hit.treeNode < 0;
      const view = ctx.graph.assembled.value?.view;
      const type = view && e.hit.layer >= 0 ? view.resolved.layers.type[e.hit.layer] : -1;
      plot = type >= 0 && view?.set.strings.get(type) === "dot_plot";
      apply();
    });
    // A window's body and edges are not hover elements: its part is followed per move (setCursor emits on change).
    on(ctx.interaction.pointer, (e) => {
      if (e.input.type === PointerType.down) {
        const ref = ctx.layout.sectionRef(e.hit.section);
        pressed = ref === null ? "default" : ctx.ui.cursor(ref).value;
        // As before the press (the cursor shown), which a window's own press handler (first) may have moved.
        const seen = e.hit.section === section && region === e.hit.region;
        heldPart = seen ? part : partAt(e.hit.section, e.hit.region, e.hit);
        heldKind = seen ? kind : found;
        if (heldPart === WindowPart.page && heldKind !== 2) heldPart = WindowPart.pan;
        return;
      }
      const next = partAt(e.hit.section, e.hit.region, e.hit);
      if (next === part && e.hit.section === section) return;
      section = e.hit.section;
      region = e.hit.region;
      kind = found;
      part = next;
      apply();
    });
    on(ctx.interaction.leave, (e) => {
      region = HitRegion.none;
      kind = -1;
      part = WindowPart.none;
      offNode = plot = false;
      if (e.section === section) apply();
    });
    add(ctx.ui.drag.subscribe(apply));
    add(() => {
      for (const ref of ctx.viewports.sections.value) ctx.ui.setCursor(ref, "default");
      ctx.ui.setWindow(-1, -1, WindowPart.none);
    });
  },
);

export const FocusViewportBehavior = defineBehavior(
  "focusViewport",
  "view",
  0,
  { animate: true, pad: 0.1 },
  "Moves the viewports of the sections that show a changed focus, over its transition.",
  ({ ctx, h, on }) => {
    on(ctx.focus.changed, (c) => {
      const p = h.forSection(ctx.layout.sectionIndex(c.section));
      const vp = ctx.viewports.get(c.section);
      if (!p || !vp || !c.entry) return;
      const opts = {
        animate: p.animate && c.transitionMs !== 0,
        durationMs: c.transitionMs > 0 ? c.transitionMs : undefined,
        origin: c.origin,
        pad: p.pad,
      };
      if (c.entry.blocks) vp.zoomToBlocks(c.entry.blocks, opts);
      else vp.zoomTo(c.entry.block, c.entry.start, c.entry.stop, opts);
    });
  },
);

export const ScrollGroupsBehavior = defineBehavior(
  "scrollGroups",
  "view",
  0,
  {},
  "Copies the zoom and horizontal scroll among a scroll group's members; each keeps its vertical scroll.",
  ({ ctx, add }) => {
    let subs: Unsubscribable[] = [];
    const release = () => {
      for (const s of subs) s.unsubscribe();
      subs = [];
    };
    add(
      ctx.scrollGroups.groups.subscribe((groups) => {
        release();
        for (const members of groups.values()) {
          // The member moved last by anyone but a reassembly (seqspec:interaction#replays-and-carries): a `graph`
          // change leads nothing, and the member it moved is set again from this one.
          let leader: (typeof members)[number] | null = null;
          for (const m of members) {
            let replay = true;
            subs.push(
              m.live.version.subscribe(() => {
                // A change made by this behavior is never copied, so nothing echoes.
                const origin = m.live.current.origin;
                if (replay || origin === "behavior") return;
                if (origin !== "graph") leader = m;
                else if (leader && leader !== m) return void m.follow(leader, "behavior");
                for (let i = 0; i < members.length; i++) if (members[i] !== m) members[i].follow(m, "behavior");
              }),
            );
            replay = false;
          }
        }
      }),
    );
    add(release);
  },
);

export const SnapshotControlsBehavior = defineBehavior(
  "snapshotControls",
  "view",
  0,
  { stopOnInput: true },
  "A press or a key on a surface stops the player.",
  ({ ctx, h, on }) => {
    const stop = (section: number) => {
      if (ctx.snapshots.playing.value && h.forSection(section)?.stopOnInput) ctx.snapshots.stop();
    };
    on(ctx.interaction.pointer, (e) => {
      if (e.input.type === PointerType.down) stop(e.input.section);
    });
    on(ctx.interaction.key, (e) => stop(e.input.section));
  },
);
