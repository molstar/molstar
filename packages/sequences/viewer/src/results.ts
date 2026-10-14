/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Results (seqspec:interaction#results, LS2): coordinates an application found, shown and stepped through with no
// search logic. Every named set lights on `results:<name>` (its loci as one union per kind, its `Coordinates` merged),
// the current result is focused and accented on `results:<name>:current`; several sets coexist and `step` walks the
// current one. Results are held as given and carried to every new set as the selection's loci are. Bookmarks step
// through a set of their own (`bookmarks:<section>`, not lit: their band lights them). Edge rate: a set, a step.

import type { SeqSet } from "@molstar/seq-core/model";
import { SeqSelection } from "@molstar/seq-core/model/query";
import type { FocusManager } from "./focus.ts";
import type { HighlightManager } from "./highlight.ts";
import { carry, type SeqLociAny } from "./lineage.ts";
import { type Coordinates, isCoordinates } from "./projection.ts";
import { type Origin, type SeqState, state } from "./state.ts";

export type Result = SeqLociAny | Coordinates;
export interface ResultsState {
  readonly name: string;
  readonly count: number;
  /** The current result; -1: none yet. */
  readonly index: number;
  readonly labels: readonly string[];
  /** Lit and focused in this section only; null: every section that shows them. */
  readonly section: string | null;
}
export interface ResultsOptions {
  /** How every result lights (default `highlight`); the current one is always accented. */
  readonly style?: "highlight" | "accent";
  readonly labels?: readonly string[];
  readonly section?: string | null;
  /** False: not lit as a set (another channel lights them); the current one still is. Default true. */
  readonly light?: boolean;
}
export interface ResultsManager {
  /** Replaces the set `name` (null clears it); its index starts at -1. */
  set(name: string, results: readonly Result[] | null, opts?: ResultsOptions): void;
  /** The last set or stepped; null: none. */
  readonly current: SeqState<ResultsState | null>;
  /** Every set by name. */
  readonly sets: SeqState<ReadonlyMap<string, ResultsState>>;
  /** The next result of `name` (default the current set's), round: focused and accented. */
  step(dir: 1 | -1, opts?: { name?: string; origin?: Origin }): void;
  /**
   * Result `index` of `name`: focused and accented, and the current one for `step`. Going to the current one again
   * clears it (`clearCurrent`), so a click on the current result's pill unselects it.
   */
  go(name: string, index: number, opts?: { origin?: Origin }): void;
  /** No current result in `name`: its accent goes, the focus stays; `step` starts again from the first or last. */
  clearCurrent(name: string): void;
  /** The set's results as held (carried to the shown set). */
  results(name: string): readonly Result[];
  clear(name: string): void;
}

interface Held {
  results: readonly Result[];
  state: ResultsState;
  readonly style: "highlight" | "accent";
  readonly light: boolean;
}

/** Results as few targets as their kinds allow: loci of the set as one per kind, `Coordinates` merged into one. */
function lit(results: readonly Result[], set: SeqSet | null): Result[] {
  const loci = results.filter((r): r is SeqLociAny => !isCoordinates(r) && r.set === set);
  const coords = results.filter(isCoordinates);
  const out: Result[] = set && loci.length ? [...SeqSelection.toLoci(SeqSelection.fromLoci(set, loci))] : [];
  if (coords.length) {
    const systems: string[] = [];
    const runs: number[] = [];
    for (const c of coords)
      for (let r = 0; r + 2 < c.runs.length; r += 3) {
        const key = c.systems[c.runs[r]];
        const at = systems.indexOf(key);
        runs.push(at < 0 ? systems.push(key) - 1 : at, c.runs[r + 1], c.runs[r + 2]);
      }
    out.push({ systems, runs: Int32Array.from(runs) });
  }
  return out;
}

export class Results implements ResultsManager {
  readonly current = state<ResultsState | null>(null);
  readonly sets = state<ReadonlyMap<string, ResultsState>>(new Map());
  private readonly held = new Map<string, Held>();
  /** The channels each set lights (`results:<name>`, then `results:<name>/<k>` for a second kind on). */
  private readonly lit = new Map<string, number>();
  private readonly highlight: HighlightManager;
  private readonly focus: FocusManager;
  /** The shown set, as results are lit on it. */
  private readonly shown: () => SeqSet | null;
  private carriedTo: SeqSet | null = null;

  constructor(highlight: HighlightManager, focus: FocusManager, shown: () => SeqSet | null) {
    this.highlight = highlight;
    this.focus = focus;
    this.shown = shown;
  }

  set(name: string, results: readonly Result[] | null, opts: ResultsOptions = {}) {
    if (!results) return this.clear(name);
    const st = { name, count: results.length, index: -1, labels: opts.labels ?? [], section: opts.section ?? null };
    this.held.set(name, { results, state: st, style: opts.style ?? "highlight", light: opts.light ?? true });
    this.highlight.clear(`results:${name}:current`);
    this.light(name);
    this.publish(name);
  }

  step(dir: 1 | -1, opts: { name?: string; origin?: Origin } = {}) {
    const h = this.held.get(opts.name ?? this.current.value?.name ?? "");
    const n = h?.state.count ?? 0;
    if (!h || n === 0) return;
    const i = h.state.index;
    this.show(h.state.name, i < 0 ? (dir > 0 ? 0 : n - 1) : (i + dir + n) % n, opts);
  }

  go(name: string, index: number, opts: { origin?: Origin } = {}) {
    if (this.held.get(name)?.state.index === index) this.clearCurrent(name);
    else this.show(name, index, opts);
  }

  clearCurrent(name: string) {
    const h = this.held.get(name);
    if (!h || h.state.index < 0) return;
    h.state = { ...h.state, index: -1 };
    this.highlight.clear(`results:${name}:current`);
    this.publish(name);
  }

  /** Result `index` current: focused and accented (a step of one result goes round to itself, not clearing it). */
  private show(name: string, index: number, opts: { origin?: Origin }) {
    const h = this.held.get(name);
    const r = h?.results[index];
    if (!h || !r) return;
    h.state = { ...h.state, index };
    const section = h.state.section;
    this.highlight.set(`results:${name}:current`, r, { style: "accent", section });
    this.focus.fromLoci(r, { ...(section ? { section } : {}), origin: opts.origin ?? "api" });
    this.publish(name);
  }

  results(name: string): readonly Result[] {
    return this.held.get(name)?.results ?? [];
  }

  clear(name: string) {
    if (!this.held.delete(name)) return;
    this.highlight.clear(`results:${name}:current`);
    this.unlight(name, 0);
    const next = new Map(this.sets.value);
    next.delete(name);
    this.sets.set(next);
    if (this.current.value?.name === name) this.current.set(null);
  }

  /**
   * Every result carried to the shown set (dropped when nothing of it is left), and the sets lit again on it. A carry
   * is no one's move (seqspec:interaction#replays-and-carries): the sets are published once, `current` stays the set
   * last set or stepped, and only a set whose count changed loses its current result and accent.
   */
  update(set: SeqSet | null) {
    if (set === this.carriedTo) return;
    this.carriedTo = set;
    if (this.held.size === 0) return;
    const sets = new Map(this.sets.value);
    for (const [name, h] of this.held) {
      if (set)
        h.results = h.results.flatMap((r): Result[] => {
          const c = isCoordinates(r) ? r : carry(r, set);
          return c ? [c] : [];
        });
      if (h.results.length !== h.state.count) {
        h.state = { ...h.state, count: h.results.length, index: -1 };
        this.highlight.clear(`results:${name}:current`);
      }
      this.light(name);
      sets.set(name, h.state);
    }
    this.sets.set(sets, "graph");
    const current = this.current.value && this.held.get(this.current.value.name);
    if (current && current.state !== this.current.value) this.current.set(current.state, "graph");
  }

  private light(name: string) {
    const h = this.held.get(name) as Held;
    const targets = h.light ? lit(h.results, this.shown()) : [];
    const section = h.state.section;
    for (let k = 0; k < targets.length; k++)
      this.highlight.set(k === 0 ? `results:${name}` : `results:${name}/${k}`, targets[k], { style: h.style, section });
    this.unlight(name, targets.length);
  }
  private unlight(name: string, from: number) {
    for (let k = from; k < (this.lit.get(name) ?? 0); k++)
      this.highlight.clear(k === 0 ? `results:${name}` : `results:${name}/${k}`);
    if (from > 0) this.lit.set(name, from);
    else this.lit.delete(name);
  }
  private publish(name: string) {
    const st = (this.held.get(name) as Held).state;
    this.sets.set(new Map(this.sets.value).set(name, st));
    this.current.set(st);
  }

  dispose() {
    this.current.complete();
    this.sets.complete();
  }
}
