/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The behavior manager (seqspec:viewer#behaviors): policy over base events, states and commands, with params, added and
// removed at runtime. Lifecycle as Mol*'s PluginBehavior: `create` builds it, `register` subscribes, `dispose`
// releases.

import { Subscription } from 'rxjs';
import type { SeqViewerContext } from './context.js';
import type { PresetEntry } from './presets.js';
import { type SeqState, state, type WritableSeqState } from './state.js';

export interface SeqBehaviorProvider<P extends object = object> {
  /** Built-ins bare ("select"); an extension's `<ns>:<name>`. */
  readonly name: string;
  readonly display: { readonly name: string; readonly description?: string };
  readonly category: 'interaction' | 'view' | 'link' | 'commands';
  /** Its order in base events. */
  readonly priority: number;
  /** Defaults; data only. */
  readonly params: P;
  create(ctx: SeqViewerContext, handle: BehaviorHandle<P>): SeqBehavior<P>;
}
export interface SeqBehavior<P> {
  register(): void;
  /** Applies params in place; false or absent: the behavior is disposed and created again. */
  update?(params: P): boolean;
  dispose(): void;
}
export interface BehaviorHandle<P> {
  readonly name: string;
  readonly params: SeqState<P>;
  setParams(p: Partial<P>): void;
  /** With the section's preset overrides (`PresetEntry.behaviors`); null: off there. Cached per section. */
  forSection(section: number): P | null;
  /** Idempotent. */
  remove(): void;
}
export interface BehaviorInfo {
  readonly name: string;
  readonly category: string;
}
export interface BehaviorManager {
  readonly list: SeqState<readonly BehaviorInfo[]>;
  /** A taken name throws. A behavior that throws is disposed and reported; the returned remove still works. */
  add<P extends object>(b: SeqBehaviorProvider<P>, params?: Partial<P>): () => void;
  get<P extends object>(name: string): BehaviorHandle<P> | null;
}

export const SeqBehaviors = {
  /** Collects a behavior's subscriptions (events, states, commands) so `dispose` releases them in one call. */
  handler: () => new Subscription(),
};

class Handle<P extends object> implements BehaviorHandle<P> {
  readonly name: string;
  readonly provider: SeqBehaviorProvider<P>;
  readonly paramsState: WritableSeqState<P>;
  behavior: SeqBehavior<P> | null = null;
  removed = false;
  /** Per section number: the preset and params the merged params were built from. */
  private readonly presets: (PresetEntry | null)[] = [];
  private readonly bases: P[] = [];
  private readonly merged: (P | null)[] = [];
  private readonly manager: Behaviors;

  constructor(manager: Behaviors, provider: SeqBehaviorProvider<P>, params: P) {
    this.manager = manager;
    this.provider = provider;
    this.name = provider.name;
    this.paramsState = state(params);
  }
  get params(): SeqState<P> {
    return this.paramsState;
  }
  setParams(p: Partial<P>) {
    if (this.removed) return;
    const prev = this.paramsState.value;
    this.paramsState.set({ ...prev, ...p });
    if (this.paramsState.value !== prev) this.manager.apply(this, this.paramsState.value);
  }
  forSection(section: number): P | null {
    const params = this.paramsState.value;
    if (section < 0) return params;
    const preset = this.manager.presetOf(section);
    if (this.bases[section] === params && this.presets[section] === preset) return this.merged[section];
    const over = preset?.behaviors?.[this.name];
    const out = over === false ? null : over ? ({ ...params, ...over } as P) : params;
    this.presets[section] = preset;
    this.bases[section] = params;
    this.merged[section] = out;
    return out;
  }
  remove() {
    this.manager.remove(this);
  }
}

export class Behaviors implements BehaviorManager {
  private readonly handles = new Map<string, Handle<object>>();
  private readonly listState = state<readonly BehaviorInfo[]>([]);
  private readonly ctx: SeqViewerContext;
  private readonly onError: (name: string, error: unknown) => void;

  constructor(ctx: SeqViewerContext, onError: (name: string, error: unknown) => void) {
    this.ctx = ctx;
    this.onError = onError;
  }

  presetOf(section: number): PresetEntry | null {
    return this.ctx.layout.presetOf(section);
  }

  get list(): SeqState<readonly BehaviorInfo[]> {
    return this.listState;
  }

  add<P extends object>(b: SeqBehaviorProvider<P>, params?: Partial<P>): () => void {
    if (this.ctx.disposed) throw new Error('SeqViewerContext: used after dispose');
    if (this.handles.has(b.name)) throw new Error(`behaviors: "${b.name}" is already added`);
    const h = new Handle<P>(this, b, { ...b.params, ...params });
    this.handles.set(h.name, h as unknown as Handle<object>);
    this.publish();
    this.start(h);
    return () => h.remove();
  }

  get<P extends object>(name: string): BehaviorHandle<P> | null {
    return (this.handles.get(name) as unknown as Handle<P> | undefined) ?? null;
  }

  /** Creates and registers; a throw disposes it and reports. */
  private start<P extends object>(h: Handle<P>) {
    try {
      h.behavior = h.provider.create(this.ctx, h);
      h.behavior.register();
    } catch (e) {
      this.fail(h, e);
    }
  }

  apply<P extends object>(h: Handle<P>, params: P) {
    const b = h.behavior;
    if (!b) return;
    try {
      if (b.update?.(params)) return;
    } catch (e) {
      this.fail(h, e);
      return;
    }
    this.stop(h);
    this.start(h);
  }

  private stop(h: Handle<object>) {
    const b = h.behavior;
    h.behavior = null;
    try {
      b?.dispose();
    } catch (e) {
      this.onError(h.name, e);
    }
  }

  private fail(h: Handle<object>, error: unknown) {
    this.onError(h.name, error);
    this.remove(h);
  }

  remove(h: Handle<object>) {
    if (h.removed) return;
    h.removed = true;
    this.stop(h);
    h.paramsState.complete();
    if (this.handles.get(h.name) === h) this.handles.delete(h.name);
    this.publish();
  }

  private publish() {
    this.listState.set([...this.handles.values()].map((h) => ({ name: h.name, category: h.provider.category })));
  }

  /** In reverse creation order; completes the list. */
  dispose() {
    for (const h of [...this.handles.values()].reverse()) this.remove(h);
    this.listState.complete();
  }
}
