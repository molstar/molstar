/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `Search` part's logic, one for both UIs (seqspec:interaction#the-default-ui,
// seqspec:interaction#slots-and-chrome; LS2): the UI spec's `search(text, signal)` is called once the text settles, the
// previous call aborted, and its answer goes to `ctx.results.set("search", …)`; the part steps through it (Enter,
// Shift+Enter, its buttons) and keeps it as a location set (`results.save`), which it can also remove again ("Clear
// kept", where it was made, so removing never needs the legend open). The viewer ships no search of its own:
// what a hit is stays the application's.

import { sectionHighlightEntries } from "./chrome.js";
import type { SeqViewerContext } from "./context.js";
import type { Result } from "./results.js";
import { state } from "./state.js";

/** What the part's "Keep" button does, as its hint says. */
export const KEEP_HINT = "Keep these results highlighted";
/** What its "Clear kept" button does. */
export const CLEAR_KEPT_HINT = "Remove every kept search highlight";
/** The results set the part fills. */
export const SEARCH_RESULTS = "search";
/** How long the text rests before the search runs. */
export const SEARCH_SETTLE_MS = 200;

/** The UI spec's search: an application's, over `@molstar/seq-core/algorithms/search` or anything else. */
export type SearchFn = (
  text: string,
  signal: AbortSignal,
) => Promise<{
  readonly results: readonly Result[];
  readonly labels?: readonly string[];
  /** The search stopped at its cap: there are more than it returned (the status says `n+`). */
  readonly more?: boolean;
}>;

/** What the part shows. */
export interface SearchView {
  readonly text: string;
  /** A search is in flight (or waiting for the text to settle). */
  readonly pending: boolean;
  readonly count: number;
  /** The current result; -1: none yet. */
  readonly index: number;
  /** The search's last failure; null: none. */
  readonly error: string | null;
  /** The answer was capped (`SearchFn`'s `more`). */
  readonly more: boolean;
  /** Sets kept by a search in the shown view's highlights ("Clear kept" removes them). */
  readonly kept: number;
}

export class SearchBox {
  readonly view = state<SearchView>({
    text: "",
    pending: false,
    count: 0,
    index: -1,
    error: null,
    more: false,
    kept: 0,
  });
  private readonly ctx: SeqViewerContext;
  private readonly search: SearchFn;
  private timer = -1;
  private abort: AbortController | null = null;
  private readonly subs: { unsubscribe(): void }[];

  constructor(ctx: SeqViewerContext, search: SearchFn) {
    this.ctx = ctx;
    this.search = search;
    this.subs = [
      ctx.results.sets.subscribe((sets) => {
        const st = sets.get(SEARCH_RESULTS);
        this.view.update({ count: st?.count ?? 0, index: st?.index ?? -1 });
      }),
      // Edge-rate: a new view (a set kept or removed, here or in the legend) recounts them.
      ctx.graph.assembled.subscribe(() => {
        const kept = this.keptSets().length;
        if (kept !== this.view.value.kept) this.view.update({ kept });
      }),
    ];
  }

  /** The box's text changed: searched once it settles; empty clears the results. */
  input(text: string): void {
    const s = this.ctx.scheduler;
    if (this.timer >= 0) s.clearTimer(this.timer);
    this.abort?.abort();
    this.abort = null;
    this.view.update({ text, pending: text.length > 0, error: null, more: false });
    if (text.length === 0) {
      this.timer = -1;
      this.ctx.results.clear(SEARCH_RESULTS);
      return;
    }
    this.timer = s.setTimer(() => this.run(text), SEARCH_SETTLE_MS);
  }

  private async run(text: string) {
    this.timer = -1;
    const abort = new AbortController();
    this.abort = abort;
    try {
      const r = await this.search(text, abort.signal);
      if (abort.signal.aborted || this.ctx.disposed) return;
      this.ctx.results.set(SEARCH_RESULTS, r.results, { labels: r.labels });
      this.view.update({ pending: false, error: null, more: r.more === true });
    } catch (e) {
      if (abort.signal.aborted) return;
      this.view.update({ pending: false, error: e instanceof Error ? e.message : String(e) });
    } finally {
      if (this.abort === abort) this.abort = null;
    }
  }

  step(dir: 1 | -1): void {
    this.ctx.results.step(dir, { name: SEARCH_RESULTS, origin: "keyboard" });
  }

  /**
   * "Keep": the results saved into the current view's first section as a set named after the query ("A", "zinc"); the
   * same query kept again replaces its set.
   */
  private saveTarget() {
    const section = this.ctx.viewports.sections.value[0];
    const text = this.view.value.text.trim();
    if (!section) return null;
    const slug = text.toLowerCase().replace(/[^a-z0-9_-]+/g, "-");
    return { section, name: SEARCH_RESULTS, ref: `${section}.keep.${slug || "search"}`, label: text || "Search" };
  }
  canSave(): boolean {
    const t = this.saveTarget();
    return !!t && this.ctx.actions.can("results.save", t);
  }
  save(): boolean {
    const t = this.saveTarget();
    return !!t && this.ctx.actions.run("results.save", t, "pointer");
  }

  /** The highlights entries a search kept (their data `<section>.keep.<query>`), in the shown view. */
  private keptSets() {
    return sectionHighlightEntries(this.ctx).filter((e) => e.data.startsWith(`${e.section}.keep.`));
  }
  /** Removes every kept search set: one `highlight.remove` per section, one patch, so one Undo brings them back. */
  clearKept(): void {
    const by = new Map<string, string[]>();
    for (const e of this.keptSets()) by.set(e.section, [...(by.get(e.section) ?? []), e.data]);
    for (const [section, data] of by) this.ctx.actions.run("highlight.remove", { section, data }, "pointer");
  }

  /** "2 / 14" ("2 / 2000+" when capped), "no results", "searching…", or the failure. */
  status(v: SearchView = this.view.value): string {
    if (v.error) return v.error;
    if (v.pending) return "searching…";
    if (v.text.length === 0) return "";
    if (v.count === 0) return "no results";
    return `${v.index < 0 ? "–" : v.index + 1} / ${v.count}${v.more ? "+" : ""}`;
  }

  dispose(): void {
    if (this.timer >= 0) this.ctx.scheduler.clearTimer(this.timer);
    this.abort?.abort();
    for (const s of this.subs) s.unsubscribe();
    this.view.complete();
  }
}
