/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Snapshots with their player (seqspec:viewer#snapshots-and-playback): the snapshot layer of the document
// (seqspec:viewer#the-document), entered with its view and its `transition_ms`, and stepped on the scheduler's clock.
// `current` mirrors the layer, so a snapshot the document leaves (it no longer applies) is left here too, once.

import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { Issue } from "@molstar/seq-core/model";
import type { Documents, Views } from "./document.js";
import type { Scheduler } from "./scheduler.js";
import { type Origin, type SeqState, state } from "./state.js";

export interface SnapshotInfo {
  readonly ref: string;
  readonly title: string | null;
  readonly description: string | null;
  readonly durationMs: number;
  readonly transitionMs: number;
  /** False after a base patch it no longer applies over (`svs.snapshot.invalid`): the player skips it. */
  readonly valid: boolean;
}
export interface SnapshotManager {
  readonly list: SeqState<readonly SnapshotInfo[]>;
  /** Null: the base. */
  readonly current: SeqState<string | null>;
  readonly playing: SeqState<boolean>;
  apply(ref: string | null, opts?: { animate?: boolean }): Promise<void>;
  next(dir: 1 | -1, opts?: { animate?: boolean }): Promise<void>;
  play(opts?: { loop?: boolean }): void;
  stop(): void;
}

export interface SnapshotHost {
  readonly documents: Documents;
  readonly views: Views;
  readonly scheduler: Scheduler;
  /** Resolves when the graph settles. */
  update(): Promise<void>;
  /** While a snapshot is applied: the transition its focus changes carry (focus.ts); null after. */
  enter(transition: { ms: number; origin: Origin } | null): void;
  report(issue: Issue): void;
}

/**
 * The view a snapshot enters: when its ops set `root.default_view`, that value (null: the first view), detected by the
 * key's presence, not by comparing with the base's; undefined: it sets none (`apply` then returns to the base's view
 * only from the view the snapshot it leaves set).
 */
function viewOf(base: SeqViewSpecDoc | null, ref: string): string | null | undefined {
  const root = base?.root.ref ?? "root";
  let out: string | null | undefined;
  for (const op of base?.snapshots?.find((s) => s.ref === ref)?.ops ?? []) {
    if (op.op === "set" && op.ref === root && op.params && Object.hasOwn(op.params, "default_view"))
      out = (op.params.default_view as string | null) ?? null;
  }
  return out;
}

export class Snapshots implements SnapshotManager {
  readonly list = state<readonly SnapshotInfo[]>([]);
  readonly current = state<string | null>(null);
  readonly playing = state(false);
  private readonly host: SnapshotHost;
  private timer = -1;
  private loop = false;
  private paused = false;
  /** No step or timer after dispose. */
  private disposed = false;

  constructor(host: SnapshotHost) {
    this.host = host;
    const d = host.documents;
    d.base.subscribe((doc) => this.list.set(infos(doc, d.invalidSnapshots())));
    d.snapshot.subscribe((ref) => this.current.set(ref, d.snapshot.origin ?? "graph"));
  }

  /** Over the base (never over another snapshot); its view and focus changes enter with its `transition_ms`. */
  async apply(ref: string | null, opts: { animate?: boolean } = {}) {
    if (this.disposed) return;
    const h = this.host;
    const info = this.list.value.find((s) => s.ref === ref);
    if (ref !== null && !info?.valid) {
      const k = this.list.value.indexOf(info as SnapshotInfo);
      const issue = info
        ? Issue.warning("svs.snapshot.invalid", `/snapshots/${k}`, `snapshot "${ref}" no longer applies`)
        : Issue.error("svs.ref.unresolved", "/snapshots", `no snapshot "${ref}"`);
      return h.report(issue);
    }
    h.enter({ ms: opts.animate === false ? 0 : (info?.transitionMs ?? 0), origin: "snapshot" });
    const base = h.documents.base.value;
    const left = this.current.value === null ? undefined : viewOf(base, this.current.value);
    const first = () => h.views.list.value[0]?.ref ?? "";
    // Still on the view the snapshot it leaves set (the reader has not picked another).
    const onLeft = left !== undefined && h.views.current.value === (left ?? first());
    h.documents.setSnapshot(ref, "snapshot");
    // Leaving a snapshot that set a view for one (or the base) that sets none returns to the base's view: the change
    // between the two is a `default_view` change too (seqspec:format#snapshots).
    let view = ref === null ? undefined : viewOf(base, ref);
    if (view === undefined && onLeft) view = (base?.root.params?.default_view as string | null) ?? null;
    if (view !== undefined) h.views.setCurrent(view ?? first(), "snapshot");
    try {
      await h.update();
    } finally {
      h.enter(null);
    }
  }

  async next(dir: 1 | -1, opts: { animate?: boolean } = {}) {
    const ref = this.neighbor(dir, false);
    if (ref !== undefined) await this.apply(ref, opts);
  }

  /** Stays `duration_ms` on each, timed from the end of its transition and graph update; skips invalid ones. */
  play(opts: { loop?: boolean } = {}) {
    if (this.playing.value || this.disposed) return;
    this.loop = opts.loop === true;
    this.playing.set(true);
    if (this.current.value === null) void this.advance();
    else this.dwell(this.info()?.durationMs ?? 0);
  }
  stop() {
    this.clear();
    if (!this.disposed) this.playing.set(false);
  }
  /** `setActive`: a hidden viewer's player pauses, and resumes its current snapshot's dwell. */
  setActive(on: boolean) {
    this.paused = !on;
    if (!on) this.clear();
    else if (this.playing.value && this.timer < 0) this.dwell(this.info()?.durationMs ?? 0);
  }

  private readonly advance = async () => {
    this.timer = -1;
    if (this.disposed) return;
    const ref = this.neighbor(1, this.loop);
    if (ref === undefined || ref === null) return this.stop();
    await this.apply(ref);
    const info = this.info();
    if (this.playing.value && info) this.dwell(info.transitionMs + info.durationMs);
  };
  private dwell(ms: number) {
    this.clear();
    if (!this.paused && !this.disposed) this.timer = this.host.scheduler.setTimer(this.advance, ms);
  }
  private clear() {
    if (this.timer >= 0) this.host.scheduler.clearTimer(this.timer);
    this.timer = -1;
  }
  private info() {
    return this.list.value.find((s) => s.ref === this.current.value);
  }
  /** The next valid snapshot in a direction; undefined at an end (wrapping with `loop`). */
  private neighbor(dir: 1 | -1, loop: boolean): string | undefined {
    const list = this.list.value;
    const at = list.findIndex((s) => s.ref === this.current.value);
    for (let k = 1; k <= list.length; k++) {
      let i = (at < 0 && dir < 0 ? list.length : at) + dir * k;
      if (i < 0 || i >= list.length) {
        if (!loop) return undefined;
        i = (i + list.length * 2) % list.length;
      }
      if (list[i].valid) return list[i].ref;
    }
    return undefined;
  }

  dispose() {
    this.disposed = true;
    this.clear();
    this.playing.set(false);
    this.list.complete();
    this.current.complete();
    this.playing.complete();
  }
}

/** The base's snapshots, valid unless its validation or its last patch said otherwise (no check of their own). */
function infos(doc: SeqViewSpecDoc | null, invalid: ReadonlySet<string>): SnapshotInfo[] {
  return (doc?.snapshots ?? []).map((s) => ({
    ref: s.ref,
    title: s.title ?? null,
    description: s.description ?? null,
    durationMs: s.duration_ms ?? 1000,
    transitionMs: s.transition_ms ?? 0,
    valid: !invalid.has(s.ref),
  }));
}
