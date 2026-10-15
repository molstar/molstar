/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Fixtures for input, behaviors and labels (IN2): LY1's opened contexts with the input records a host or `/canvas`
// would fill, and the subscriptions a removed behavior must leave.

import type { Loaders, SeqViewSpecDoc } from '@molstar/seq-core/format';
import { fixtureLoaders, MolstarFixture } from '../../core/_test/format/load-fixtures.js';
import {
  createSeqViewer,
  DefaultRegistry,
  DefaultSeqViewerSpec,
  Inputs,
  type LayerRenderer,
  Modifier,
  manualScheduler,
  PointerType,
  type SeqViewerContext,
  type SeqViewerRegistryEntry,
  type SeqViewerSpec,
} from '../src/index.js';
import type { Layouts } from '../src/layout.js';
import { fixtureMeasurer } from '../src/testing/index.js';
import { WINDOW_FRAME, windowFrame } from '../src/windows.js';
import { helpers, PRESETS, renderers } from './layout-fixtures.js';
export interface DrivenOptions {
  size?: Record<string, readonly [number, number]>;
  registry?: SeqViewerRegistryEntry[];
  renderers?: LayerRenderer[];
  loaders?: Loaders;
  /** False: without the built-in label providers (their strings are one allocation per hovered element by design). */
  labels?: false;
  /** The viewer spec's behaviors as a host composes them (default: `DefaultSeqViewerSpec`'s). */
  behaviors?: (list: SeqViewerSpec['behaviors']) => SeqViewerSpec['behaviors'];
}
export interface PointerOpts {
  shift?: boolean;
  alt?: boolean;
  ctrl?: boolean;
  button?: number;
  pointerId?: number;
  /** Outside the surface: the engine listens on the document during a drag. */
  onDocument?: boolean;
}
const mods = (o: { shift?: boolean; alt?: boolean; ctrl?: boolean; meta?: boolean }) =>
  (o.shift ? Modifier.shift : 0) |
  (o.alt ? Modifier.alt : 0) |
  (o.ctrl ? Modifier.ctrl : 0) |
  (o.meta ? Modifier.meta : 0);

/** A context showing `doc` with every default behavior, every viewport sized (1000 × 400 unless given). */
export async function driven(doc: SeqViewSpecDoc, o: DrivenOptions = {}) {
  const scheduler = manualScheduler();
  const defaults = DefaultSeqViewerSpec();
  const spec = o.behaviors ? { ...defaults, behaviors: o.behaviors(defaults.behaviors) } : defaults;
  const ctx = createSeqViewer(
    {
      ...spec,
      registry: [
        {
          ...DefaultRegistry,
          presets: undefined,
          renderers: undefined,
          labels: o.labels === false ? undefined : DefaultRegistry.labels,
        },
        { svs: MolstarFixture },
        { renderers: o.renderers ?? renderers() },
        { presets: PRESETS },
        ...(o.registry ?? []),
      ],
    },
    { scheduler, loaders: o.loaders ?? fixtureLoaders(), text: fixtureMeasurer() },
  );
  const done = ctx.setDocument(doc);
  for (const ref of ctx.viewports.sections.value) {
    const [w, h] = o.size?.[ref] ?? [1000, 400];
    ctx.viewports.get(ref)?.setSize(w, h, 1);
  }
  const r = await done;
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  const base = helpers(ctx);
  const p = Inputs.pointer();
  const w = Inputs.wheel();
  const k = Inputs.key();
  const pointer = (type: keyof typeof PointerType, section: string, x: number, y: number, opts: PointerOpts = {}) => {
    p.type = PointerType[type];
    p.section = ctx.layout.sectionIndex(section);
    p.x = x;
    p.y = y;
    p.button = type === 'move' ? -1 : (opts.button ?? 0);
    p.buttons = type === 'up' || type === 'cancel' ? 0 : 1 << (opts.button ?? 0);
    p.modifiers = mods(opts);
    p.pointerId = opts.pointerId ?? 0;
    p.time = scheduler.now();
    ctx.input.pointer(p);
  };
  /** A press, moves in `steps` and a release. */
  const drag = (section: string, x0: number, y0: number, x1: number, y1: number, opts: PointerOpts = {}, steps = 4) => {
    pointer('down', section, x0, y0, opts);
    for (let i = 1; i <= steps; i++)
      pointer('move', section, x0 + ((x1 - x0) * i) / steps, y0 + ((y1 - y0) * i) / steps, opts);
    pointer('up', section, x1, y1, opts);
  };
  const click = (section: string, x: number, y: number, opts: PointerOpts = {}) => {
    pointer('down', section, x, y, opts);
    pointer('up', section, x, y, opts);
  };
  const wheel = (section: string, x: number, y: number, deltaX: number, deltaY: number, opts: PointerOpts = {}) => {
    Object.assign(w, { section: ctx.layout.sectionIndex(section), x, y, deltaX, deltaY, modifiers: mods(opts) });
    return ctx.input.wheel(w);
  };
  /** A key press; `code` is the physical key (`Equal`), which a modifier leaves while `key` changes (Alt+= on macOS). */
  const key = (
    section: string,
    name: string,
    opts: { shift?: boolean; alt?: boolean; ctrl?: boolean; meta?: boolean; typing?: boolean; code?: string } = {},
  ) => {
    Object.assign(k, {
      section: ctx.layout.sectionIndex(section),
      key: name,
      code: opts.code ?? '',
      modifiers: mods(opts),
      typing: !!opts.typing,
    });
    return ctx.input.key(k);
  };
  /** The surface gained (Tab into it) or lost keyboard focus. */
  const focus = (section: string, focused = true) => ctx.input.focus(ctx.layout.sectionIndex(section), focused);
  /** Moves the clock: animations and settle timers. */
  const advance = (ms: number) => scheduler.advance(ms);
  /** A section's window frame (windows.ts `windowFrame`) of `kind` (`WindowKind`), a new array. */
  const frame = (section: string, kind: number) => {
    const g = new Float64Array(WINDOW_FRAME);
    windowFrame(
      ctx.layout.axis(section) as NonNullable<ReturnType<typeof ctx.layout.axis>>,
      base.layout(section),
      kind,
      g,
    );
    return g;
  };
  return {
    ctx,
    scheduler,
    layouts: ctx.layout as Layouts,
    ...base,
    pointer,
    drag,
    click,
    wheel,
    key,
    focus,
    advance,
    frame,
  };
}
export type Driven = Awaited<ReturnType<typeof driven>>;

/** Every subscription behaviors can hold: base events, the states they follow. */
export function subscriptions(ctx: SeqViewerContext): number {
  const count = (s: object) => (s as { subscriberCount: number }).subscriberCount;
  let n = 0;
  for (const e of Object.values(ctx.input.events)) n += count(e);
  for (const s of [ctx.hover.live.version, ctx.ui.drag, ctx.focus.changed, ctx.scrollGroups.groups]) n += count(s);
  for (const s of [ctx.labels.providers, ctx.graph.assembled, ctx.ui.window.version]) n += count(s);
  for (const ref of ctx.viewports.sections.value) n += count(ctx.viewports.get(ref)?.live.version as object);
  return n;
}
