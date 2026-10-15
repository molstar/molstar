/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Commands (seqspec:viewer#commands): the core's only request mechanism, as Mol*'s PluginCommand with stable namespaced
// ids. A command is an id, a params table (the format's tables) and its subscribers; it is registered like any registry
// entry (`commands`, seqspec:viewer#registration), so `run` takes a request as plain data, an invocation
// (seqspec:format#descriptions-and-invocations). The built-ins
// and the `commands` behavior that handles them (seqspec:viewer#behaviors; from default-behaviors.ts) are here too.
// Edge-rate:
// a dispatch allocates its promise.

import {
  bool,
  type CheckContext,
  type Closure,
  caseOf,
  checkParams,
  type FieldType,
  type Invocation,
  int,
  isPlainObject,
  json,
  list,
  literal,
  nullable,
  OptionalField,
  object,
  type ParamsSchema,
  RequiredField,
  ref,
  report,
  SimpleParamsSchema,
  str,
  union,
} from '@molstar/seq-core/format';
import { Extents, Loci } from '@molstar/seq-core/model';
import type { SeqViewerContext } from './context.js';
import { defineBehavior, internal } from './default-behaviors.js';
import type { HighlightStyle } from './highlight.js';
import type { SeqLociAny } from './lineage.js';
import type { Coordinates } from './projection.js';
import type { RangeArgs } from './refs.js';
import type { SeqViewerRegistry } from './registry.js';
import type { Origin } from './state.js';
import type { Viewports } from './viewport.js';

/** A command is callable: `cmd(ctx, params)` dispatches it on that context. */
export interface SeqCommand<T> {
  readonly id: string;
  /** What `run` checks an invocation's params against; the reference lists it. */
  readonly params: ParamsSchema;
  /** What it does, for the reference; may be empty. */
  readonly description: string;
  (ctx: { readonly commands: CommandManager }, params: T): Promise<void>;
}
/** Ids are namespaced: the viewer's `seqv:*`, an extension's `<ns>:<name>`. */
export function SeqCommand<T>(id: string, params: ParamsSchema, description = ''): SeqCommand<T> {
  const cmd: SeqCommand<T> = Object.assign(
    (ctx: { readonly commands: CommandManager }, p: T) => ctx.commands.dispatch(cmd, p),
    { id, params, description },
  );
  return cmd;
}

export interface CommandManager {
  /** Typed, from code. Resolves after every subscriber's action; with none, at once. */
  dispatch<T>(cmd: SeqCommand<T>, params: T): Promise<void>;
  /**
   * Plain data: the registered command by id, its params checked by its table, then dispatched with `origin` in the
   * params of a command whose table takes one. False: unknown id or invalid params, nothing run, nothing reported; or
   * a target's range args outside its ref's domain (`seqv.target.range`, reported). Rejects as `dispatch` does (a
   * subscriber's rejection, a disposed context).
   */
  run(invocation: Invocation, origin: Origin): Promise<boolean>;
  subscribe<T>(cmd: SeqCommand<T>, action: (params: T) => unknown): { unsubscribe(): void };
}

/** What a command acts on; in an invocation `refs` or `coordinates` (their `runs` a list of numbers), `loci` from code. */
export type SeqTarget =
  | { readonly refs: readonly string[]; readonly args?: RangeArgs }
  | { readonly loci: SeqLociAny }
  | { readonly coordinates: Coordinates };

/** `Coordinates` as data: systems, and (system, start, stop) triples into them. */
const shape = object({ systems: RequiredField(list(str), ''), runs: RequiredField(list(int), '') });
const coordinates: FieldType<Coordinates> = {
  name: shape.name,
  check(v, path, ctx) {
    const n = ctx.issues.length;
    shape.check(v, path, ctx);
    if (ctx.issues.length > n || !isPlainObject(v)) return;
    const [systems, runs] = [v.systems as string[], v.runs as number[]];
    let ok = runs.length % 3 === 0;
    for (let i = 0; ok && i < runs.length; i += 3)
      ok = runs[i] >= 0 && runs[i] < systems.length && runs[i + 1] <= runs[i + 2];
    if (!ok) report(ctx, 'svs.params.type', path, 'runs are (system index, start, stop) triples, start ≤ stop');
  },
};
const target = union(
  object({
    refs: RequiredField(
      list(ref('unit', 'sequence', 'ligand', 'complex', 'alignment', 'row', 'annotation', 'track', 'layer', 'block')),
      'document refs, derived refs included',
    ),
    args: OptionalField(
      nullable(
        object({
          start: OptionalField(nullable(int), null, 'an index'),
          stop: OptionalField(nullable(int), null, 'an index'),
          beg: OptionalField(nullable(str), null, 'a label in `numbering`'),
          end: OptionalField(nullable(str), null, 'a label in `numbering`'),
          numbering: OptionalField(nullable(str), null, "the labels' numbering"),
          items: OptionalField(nullable(list(str)), null, 'item ids of an annotation ref'),
        }),
      ),
      null,
      'narrows each ref',
    ),
  }),
  object({ coordinates: RequiredField(coordinates, 'runs in coordinate systems by key (seqspec:interaction#focus)') }),
);
const closure = OptionalField(literal<Closure>('element', 'containing'), 'element', 'closed before it lights or acts');
const origin = OptionalField(str, 'api', "the caller's: `run` sets it");
const animate = OptionalField(bool, true, 'false: at once');

/** The built-ins, handled by the `commands` behavior (below). */
export const SeqCommands = {
  Highlight: SeqCommand<{ channel: string; target: SeqTarget | null; style?: HighlightStyle; closure?: Closure }>(
    'seqv:highlight',
    SimpleParamsSchema({
      channel: RequiredField(str, "the owner's channel; others are left as they are"),
      target: RequiredField(nullable(target), 'null clears the channel'),
      style: OptionalField(literal<HighlightStyle>('highlight', 'accent', 'band'), 'highlight', ''),
      closure,
    }),
    'Lights a target on one highlight channel, projected onto every section that shows it; null clears the channel.',
  ),
  Select: SeqCommand<{
    target: SeqTarget | null;
    mode?: 'set' | 'add' | 'remove' | 'toggle';
    closure?: Closure;
    origin?: Origin;
  }>(
    'seqv:select',
    SimpleParamsSchema({
      target: RequiredField(nullable(target), 'null clears the selection'),
      mode: OptionalField(literal('set', 'add', 'remove', 'toggle'), 'set', ''),
      closure,
      origin,
    }),
    'Sets, adds to, removes from or toggles the selection; null clears it.',
  ),
  Focus: SeqCommand<{ target: SeqTarget; section?: string | null; animate?: boolean; origin?: Origin }>(
    'seqv:focus',
    SimpleParamsSchema({
      target: RequiredField(target, ''),
      section: OptionalField(nullable(str), null, 'null: every section that shows it'),
      animate,
      origin,
    }),
    'Brings a target into view in one section or in every section that shows it.',
  ),
  Reveal: SeqCommand<{ target: SeqTarget; section?: string | null; animate?: boolean; origin?: Origin }>(
    'seqv:reveal',
    SimpleParamsSchema({
      target: RequiredField(target, ''),
      section: OptionalField(nullable(str), null, 'null: every section that shows it'),
      animate,
      origin,
    }),
    'Brings a target into view without selecting, focusing or zooming: text scrolls its lines in, tracks pan at their cell; nothing where it is shown already.',
  ),
  SetSnapshot: SeqCommand<{ ref: string | null; animate?: boolean }>(
    'seqv:setSnapshot',
    SimpleParamsSchema({
      ref: RequiredField(nullable(str), 'a snapshot ref, "next" or "prev"; null: the base'),
      animate,
    }),
    'Shows a snapshot of the document, the next or the previous one, or the base.',
  ),
  SetView: SeqCommand<{ ref: string }>(
    'seqv:setView',
    SimpleParamsSchema({ ref: RequiredField(str, 'a view ref') }),
    'Switches to another view of the document; nothing is reloaded.',
  ),
  /**
   * The `layer.colorScheme` action (seqspec:interaction#color-schemes-at-runtime): refused (nothing done) for a scheme
   * that does not fit.
   */
  ColorScheme: SeqCommand<{ layer: string; scheme: string | null; write?: boolean; origin?: Origin }>(
    'seqv:colorScheme',
    SimpleParamsSchema({
      layer: RequiredField(ref('layer'), 'the layer whose color it sets'),
      scheme: RequiredField(nullable(str), 'a registered color scheme that fits the layer; null: its own color'),
      write: OptionalField(bool, false, "true: into the document (a patch on the layer's color); false: UI state"),
      origin,
    }),
    'Colors a layer by a registered color scheme, or back to its own color; refused (nothing done) for a scheme that does not fit.',
  ),
  /** The `layer.colorRange` action: the layer's color scale over [lo, hi], wider or narrower than its data; null: its own. */
  ColorRange: SeqCommand<{ layer: string; range: [number, number] | null; origin?: Origin }>(
    'seqv:colorRange',
    SimpleParamsSchema({
      layer: RequiredField(ref('layer'), 'the layer whose color range it sets'),
      range: RequiredField(json, "[lo, hi] with lo < hi; null: the layer's own range"),
      origin,
    }),
    "Sets a layer's color scale over [lo, hi], wider or narrower than its data; null: its own range.",
  ),
  /**
   * The `view.legend` action (seqspec:interaction#the-legend): the legend open, closed, or toggled; anything else is
   * refused.
   */
  Legend: SeqCommand<{ open: boolean | 'toggle'; origin?: Origin }>(
    'seqv:legend',
    SimpleParamsSchema({ open: RequiredField(json, 'true, false or "toggle"'), origin }),
    'Opens, closes or toggles the legend; any other value is refused.',
  ),
  /**
   * The `stream.retry` action (seqspec:interaction#the-default-ui): a track's failed regions requested again now, or
   * every track's (ST3b).
   */
  Retry: SeqCommand<{ track?: string | null; origin?: Origin }>(
    'seqv:retry',
    SimpleParamsSchema({
      track: OptionalField(nullable(ref('track')), null, 'null: every track with a failed region'),
      origin,
    }),
    "Requests a track's failed stream regions again now, or every track's.",
  ),
  RunAction: SeqCommand<{ name: string; target: unknown; origin?: Origin }>(
    'seqv:runAction',
    SimpleParamsSchema({
      name: RequiredField(str, 'a registered action'),
      target: OptionalField(json, null, ''),
      origin,
    }),
    'Runs a registered action on a target, as a menu or a hover card would.',
  ),
} as const;

type Action = (params: unknown) => unknown;

export class Commands implements CommandManager {
  private readonly subs = new Map<string, readonly Action[]>();
  private readonly pending = new Set<(reason: Error) => void>();
  private disposed = false;
  private readonly registry: SeqViewerRegistry;
  private readonly onError: (id: string, error: unknown) => void;

  constructor(registry: SeqViewerRegistry, onError: (id: string, error: unknown) => void) {
    this.registry = registry;
    this.onError = onError;
  }

  subscriberCount(cmd?: SeqCommand<never>): number {
    if (cmd) return this.subs.get(cmd.id)?.length ?? 0;
    let n = 0;
    for (const l of this.subs.values()) n += l.length;
    return n;
  }

  subscribe<T>(cmd: SeqCommand<T>, action: (params: T) => unknown) {
    if (this.disposed) throw new Error('commands: the context is disposed');
    const a = action as Action;
    this.subs.set(cmd.id, [...(this.subs.get(cmd.id) ?? []), a]);
    let done = false;
    return {
      unsubscribe: () => {
        if (done) return;
        done = true;
        const rest = (this.subs.get(cmd.id) ?? []).filter((x) => x !== a);
        if (rest.length > 0) this.subs.set(cmd.id, rest);
        else this.subs.delete(cmd.id);
      },
    };
  }

  run(invocation: Invocation, origin: Origin): Promise<boolean> {
    if (this.disposed) return Promise.reject(new Error('commands: run on a disposed context'));
    const id = invocation?.command;
    const cmd = typeof id === 'string' ? this.registry.get('commands', id) : undefined;
    if (!cmd) return Promise.resolve(false);
    const params = invocation.params ?? {};
    const ctx: CheckContext = {
      issues: [],
      registries: this.registry.svs,
      refs: [],
      lenient: false,
      from: '',
      depth: 0,
    };
    checkParams(cmd.params, params, '', ctx);
    if (ctx.issues.some((i) => i.severity === 'error')) return Promise.resolve(false);
    const s = cmd.params;
    const fields =
      s.type === 'simple' ? s.fields : caseOf(s, (params as Record<string, unknown>)[s.discriminator])?.fields;
    const p = fields && Object.hasOwn(fields, 'origin') ? { ...params, origin } : params;
    return this.dispatch(cmd, p).then(
      () => true,
      (e: unknown) => {
        if (e instanceof TargetRangeError) return false;
        throw e;
      },
    );
  }

  dispatch<T>(cmd: SeqCommand<T>, params: T): Promise<void> {
    if (this.disposed) return Promise.reject(new Error(`commands: "${cmd.id}" dispatched on a disposed context`));
    const actions = this.subs.get(cmd.id);
    if (!actions) return Promise.resolve();
    return new Promise<void>((resolve, reject) => {
      this.pending.add(reject);
      const run = async () => {
        // In subscription order, each awaited; the first rejection rejects the dispatch.
        for (const action of actions) await action(params);
      };
      run().then(
        () => {
          this.pending.delete(reject);
          resolve();
        },
        (error: unknown) => {
          this.pending.delete(reject);
          if (!this.disposed) this.onError(cmd.id, error);
          reject(error);
        },
      );
    });
  }

  /** Rejects every pending dispatch; later dispatches reject. */
  dispose() {
    this.disposed = true;
    for (const reject of this.pending) reject(new Error('commands: the context was disposed'));
    this.pending.clear();
    this.subs.clear();
  }
}

/** A command's range args outside its ref's domain: reported as `seqv.target.range`, and `run` resolves false. */
export class TargetRangeError extends Error {}

/**
 * A target as the managers take it: loci (refs with range args resolved and joined by kind: one loci a kind, an
 * alignment's apart, in the order first named) or `Coordinates` (runs as data made typed); null clears; undefined: it
 * resolves to nothing (logged, not reported to the reader). Range args outside a ref's domain throw
 * `TargetRangeError` (never clipped, as a document's focus is refused).
 */
function targetOf(ctx: SeqViewerContext, t: SeqTarget | null): readonly SeqLociAny[] | Coordinates | null | undefined {
  if (!t) return null;
  if ('loci' in t) return [t.loci];
  if ('coordinates' in t) {
    const { systems, runs } = t.coordinates;
    return { systems, runs: runs instanceof Int32Array ? runs : Int32Array.from(runs) };
  }
  const out: SeqLociAny[] = [];
  const outside = (message: string) => {
    throw new TargetRangeError(`${message}: give start and stop within it`);
  };
  for (const ref of t.refs) {
    const l = ctx.refs.toLoci(ref, t.args ?? undefined, outside);
    if (!l) continue;
    const k = out.findIndex(
      (o) => o.kind === l.kind && (!('alignment' in o) || o.alignment === (l as typeof o).alignment),
    );
    if (k < 0) out.push(l);
    else out[k] = Loci.union(out[k], l as never) as SeqLociAny;
  }
  if (out.length === 0) internal(ctx).log('info', `a command's target resolves to nothing: ${t.refs.join(', ')}`);
  return out.length ? out : undefined;
}
const isList = (t: readonly SeqLociAny[] | Coordinates): t is readonly SeqLociAny[] => Array.isArray(t);
/** A highlight channel's parts: `channel` for the first kind, `channel/1`, … for the others. */
const partChannel = (channel: string, k: number) => (k === 0 ? channel : `${channel}/${k}`);

export const CommandsBehavior = defineBehavior(
  'commands',
  'commands',
  0,
  {},
  'Handles the built-in commands (SeqCommands) with the managers.',
  ({ ctx, add }) => {
    const C = SeqCommands;
    // The parts each command channel set: only those are cleared (a host's own `<channel>/1` is not a part).
    const partsOf = new Map<string, number>();
    add(
      ctx.commands.subscribe(C.Highlight, (p) => {
        const t = targetOf(ctx, p.target);
        if (t === undefined) return;
        // Several kinds: a channel each (`channel/1`, …), cleared with it; a later target's fewer parts drop the rest.
        const parts: readonly (SeqLociAny | Coordinates)[] = t === null ? [] : isList(t) ? t : [t];
        for (let k = 0; k < parts.length; k++)
          ctx.highlight.set(partChannel(p.channel, k), parts[k], { style: p.style, closure: p.closure });
        const before = Math.max(1, partsOf.get(p.channel) ?? 0);
        for (let k = parts.length; k < before; k++) ctx.highlight.clear(partChannel(p.channel, k));
        if (parts.length > 0) partsOf.set(p.channel, parts.length);
        else partsOf.delete(p.channel);
      }),
    );
    add(
      ctx.commands.subscribe(C.Select, (p) => {
        const t = targetOf(ctx, p.target);
        const [origin, mode] = [p.origin ?? 'api', p.mode ?? 'set'];
        if (t === null) return ctx.selection.clear(origin);
        if (!t) return;
        // Closed before the verb; `Coordinates` first as loci in the shown set, one per kind (a set: the first sets).
        const parts = isList(t) ? t : (ctx.viewports as Viewports).projector.lociOf(t);
        if (mode === 'set' && parts.length === 0) ctx.selection.clear(origin);
        for (let i = 0; i < parts.length; i++) {
          const loci = ctx.selection.close(parts[i], p.closure ?? 'element');
          ctx.selection[mode === 'set' && i > 0 ? 'add' : mode](loci, origin);
        }
      }),
    );
    add(
      ctx.commands.subscribe(C.Focus, (p) => {
        const t = targetOf(ctx, p.target);
        const opts = {
          section: p.section ?? undefined,
          transitionMs: p.animate === false ? 0 : undefined,
          origin: p.origin,
        };
        if (!t) return;
        let shown: boolean;
        if (!isList(t) || t.length === 1) shown = ctx.focus.fromLoci(isList(t) ? t[0] : t, opts);
        else {
          // Several kinds framed together: their projections' hull.
          const ext = Extents.create(4);
          for (const part of t) (ctx.viewports as Viewports).projector.fromLoci(part, ext);
          shown = ctx.focus.fromExtents(ext, opts);
        }
        if (!shown) internal(ctx).log('info', 'a focus target is shown nowhere');
      }),
    );
    add(
      ctx.commands.subscribe(C.Reveal, (p) => {
        const t = targetOf(ctx, p.target);
        if (!t) return;
        const opts = { section: p.section, animate: p.animate !== false, origin: p.origin };
        if (!ctx.reveal(isList(t) ? t : [t], opts)) internal(ctx).log('info', 'a reveal target is shown nowhere');
      }),
    );
    add(
      ctx.commands.subscribe(C.SetSnapshot, (p) =>
        p.ref === 'next' || p.ref === 'prev'
          ? ctx.snapshots.next(p.ref === 'next' ? 1 : -1, { animate: p.animate })
          : ctx.snapshots.apply(p.ref, { animate: p.animate }),
      ),
    );
    add(ctx.commands.subscribe(C.SetView, (p) => ctx.views.setCurrent(p.ref)));
    add(ctx.commands.subscribe(C.RunAction, (p) => ctx.actions.run(p.name, p.target, p.origin)));
    add(
      ctx.commands.subscribe(C.Retry, (p) => {
        if (p.track) return void ctx.actions.run('stream.retry', p.track, p.origin);
        const m = ctx.graph.assembled.value?.view;
        for (let k = 0; m && k < m.resolved.layers.count; k++)
          ctx.streams.retry(m.set.strings.get(m.resolved.layers.ref[k]));
      }),
    );
    add(
      ctx.commands.subscribe(C.Legend, (p) => {
        if (!ctx.actions.run('view.legend', { open: p.open }, p.origin))
          internal(ctx).log('info', `the legend's "open" is true, false or "toggle", not ${JSON.stringify(p.open)}`);
      }),
    );
    add(
      ctx.commands.subscribe(C.ColorRange, (p) => {
        if (!ctx.actions.run('layer.colorRange', { layer: p.layer, range: p.range }, p.origin))
          internal(ctx).log('info', `the color range of "${p.layer}" is [lo, hi] with lo < hi, or null`);
      }),
    );
    add(
      ctx.commands.subscribe(C.ColorScheme, (p) => {
        const t = { layer: p.layer, scheme: p.scheme, write: p.write ?? false };
        if (!ctx.actions.run('layer.colorScheme', t, p.origin))
          internal(ctx).log('info', `color scheme "${p.scheme}" does not fit layer "${p.layer}"`);
      }),
    );
  },
);
