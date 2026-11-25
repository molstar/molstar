/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The search model (_spec/README.md "Search"): DOM-free, one search a viewer. The query, the status, the hits (held
// columnar, their targets' addresses), the count past the held ones, the current hit, the issues and whether the result
// is stale; two highlight channels of its own (every hit, the current one); stepping revealed through the address model;
// a saved set a ref, one patch. What an application builds its box and its results list on, with the commands. The
// matchers are the registry's (`search-matcher`), the built-ins a lazy chunk loaded by the first run.

import { ReactiveModel, type ReadonlyState, State } from '@molstar/core/reactive';
import type { Issue } from '@molstar/core/repr/tree';
import type { Location } from '@molstar/seqx-core/query';
import type { SearchHits, SearchMatcherEntry, SearchStrands } from '@molstar/seqx-core/search';
import type { ViewerConfig } from '@molstar/seqx-viewer/config';
import type {
  AddressModel,
  Emphasized,
  HighlightModel,
  SelectionModel,
  SpecModel,
  StoreModel,
  ViewModel,
} from '@molstar/seqx-viewer/models';

/** A search's query as run. */
export interface SearchQuery {
  readonly pattern: string;
  /** A registered matcher's name (`degenerate`, `exact`, `regex`, an application's `<id>:<name>`). */
  readonly matcher: string;
  readonly strands: SearchStrands;
  /** An address, or `selection`: what the search is limited to (the whole of what the view draws when absent). */
  readonly within?: string;
}

export type SearchStatus = 'idle' | 'running' | 'done' | 'cancelled' | 'failed';

/** The held hits and their targets: a hit's `target` indexes `systems`. */
export interface SearchResult {
  readonly hits: SearchHits;
  /** The searched systems' addresses, in the search's order. */
  readonly systems: readonly string[];
  /** Their lengths (a hit across a circle's origin ends past it). */
  readonly lengths: readonly number[];
  /** Their words (a label, else an id). */
  readonly labels: readonly string[];
}

/** What a run, a step or a save says: its words and its issues. */
export interface SearchOutcome {
  readonly message: string;
  readonly issues: readonly Issue[];
}

/** What the model reads: the viewer's models, its configuration, its registry's matchers, its notices. */
export interface SearchDeps {
  readonly spec: SpecModel;
  readonly view: ViewModel;
  readonly store: StoreModel;
  readonly selection: SelectionModel;
  readonly highlight: HighlightModel;
  readonly addresses: AddressModel;
  readonly config: ViewerConfig;
  /** The registry's matchers by name (read at each run: a late extension's too). */
  matchers(): Readonly<Record<string, SearchMatcherEntry>>;
  /** Says a run's end and progress, and a step (the announcer's notice). */
  notify(text: string): void;
}

const NoResult: SearchResult = {
  hits: {
    length: 0,
    target: new Int32Array(0),
    start: new Int32Array(0),
    end: new Int32Array(0),
    strand: new Int8Array(0),
  },
  systems: [],
  lengths: [],
  labels: [],
};

const NoHits: SearchOutcome = { message: 'No hits', issues: [] };

/** The model's work, a lazy chunk the first run loads (every step comes after one). */
let ops: typeof import('./search-ops.js') | undefined;
const loadOps = async () => {
  ops ??= await import('./search-ops.js');
  return ops;
};

export class SearchModel extends ReactiveModel {
  /** @internal */
  readonly _query = new State<SearchQuery | undefined>(undefined);
  /** @internal */
  readonly _status = new State<SearchStatus>('idle');
  /** @internal */
  readonly _hits = new State<SearchResult>(NoResult);
  /** @internal */
  readonly _count = new State(0);
  /** @internal */
  readonly _complete = new State(true);
  /** @internal */
  readonly _current = new State(-1);
  /** @internal */
  readonly _issues = new State<readonly Issue[]>([]);
  /** @internal */
  readonly _stale = new State(false);

  /** The query of the last run (undefined: none, or cleared). */
  readonly query: ReadonlyState<SearchQuery | undefined> = this._query.asReadonly;
  readonly status: ReadonlyState<SearchStatus> = this._status.asReadonly;
  /** The held hits (at most `search.maxHits`), in order: what the channel draws and stepping steps through. */
  readonly hits: ReadonlyState<SearchResult> = this._hits.asReadonly;
  /** The hits found: past the held ones when the bound was reached (exact within the step budget). */
  readonly count: ReadonlyState<number> = this._count.asReadonly;
  /** Whether the result means what it says: false while running, past a bound, a system skipped, or cancelled. */
  readonly complete: ReadonlyState<boolean> = this._complete.asReadonly;
  /** The current hit's index in `hits`, −1 for none. */
  readonly current: ReadonlyState<number> = this._current.asReadonly;
  readonly issues: ReadonlyState<readonly Issue[]> = this._issues.asReadonly;
  /** New letters for a searched system since the run: still drawn and said, until the next run. */
  readonly stale: ReadonlyState<boolean> = this._stale.asReadonly;

  /** @internal */
  _run: AbortController | undefined;
  /** @internal */
  _runs = 0;
  /** The channels as last set by the model: a change it did not make (the legend's removal) clears the search. */
  /** @internal */
  _drawn: Emphasized | undefined;
  /** The searched systems' entries at the run: another entry later is new letters (stale). */
  /** @internal */
  _entries: unknown[] = [];

  /** @internal */
  readonly deps: SearchDeps;

  constructor(deps: SearchDeps) {
    super();
    this.deps = deps;
    const { highlight, spec } = deps;
    // the legend's channels row removing `search` runs `search.clear`: a channel and the search never disagree
    this.listenPermanently(highlight.channel(this.channel), (v) => {
      if (this._drawn && !v) this.clear();
    });
    // results follow the spec: a removed system's hits dropped, new letters stale, never a run of its own
    this.listenPermanently(spec.model, () => ops?.followModel(this));
  }

  /** The channel every hit is drawn on (`search.channel`). */
  get channel(): string {
    return this.deps.config.search.channel;
  }

  /** The channel the current hit is drawn on (`search.currentChannel`). */
  get currentChannel(): string {
    return this.deps.config.search.currentChannel;
  }

  /**
   * Runs a search: the pattern compiled first (a compile error is the outcome's message and issue, and nothing changes),
   * a run in progress cancelled, the result replaced as batches arrive. Resolves when done, cancelled or failed.
   */
  async run(q: {
    readonly pattern: string;
    readonly matcher?: string;
    readonly strands?: SearchStrands;
    readonly within?: string;
  }): Promise<SearchOutcome> {
    // a later run or a clear, even before this one's chunk lands, makes it moot
    const id = ++this._runs;
    const o = await loadOps();
    return id === this._runs ? o.runSearch(this, q, id) : { message: '', issues: [] };
  }

  /** Cancels a run in progress, clears the query, the hits and both channels. */
  clear(): void {
    this._run?.abort();
    this._run = undefined;
    this._runs++;
    this._query.set(undefined);
    this._status.set('idle');
    this._hits.set(NoResult);
    this._count.set(0);
    this._complete.set(true);
    this._current.set(-1);
    this._issues.set([]);
    this._stale.set(false);
    this._entries = [];
    this._drawn = undefined;
    this.deps.highlight.clear(this.channel);
    this.deps.highlight.clear(this.currentChannel);
  }

  /** The next hit becomes current and is revealed (`focus`: zoomed to); from none, the first in or after the view. */
  next(options: { readonly focus?: boolean } = {}): SearchOutcome {
    return ops?.step(this, 1, options.focus ?? false) ?? NoHits;
  }

  /** The previous hit becomes current and is revealed; from none, the last before the view's end. */
  previous(options: { readonly focus?: boolean } = {}): SearchOutcome {
    return ops?.step(this, -1, options.focus ?? false) ?? NoHits;
  }

  /**
   * That hit becomes current and is revealed; the current one again clears it (EC-161); an index outside the hits is
   * refused, never clipped (EC-160).
   */
  go(index: number, options: { readonly focus?: boolean } = {}): SearchOutcome {
    return ops?.go(this, index, options) ?? NoHits;
  }

  /** A hit's Location: its range on its system, two ranges across a circle's origin (one hit). */
  locationOf(i: number): Location | undefined {
    return ops?.locationOf(this, i);
  }

  /**
   * A hit's address: its system's, then its range in ordinal labels (copy link's numbering); a hit across a circle's
   * origin one range with `from > to` (EC-10). Made on demand, never stored a hit.
   */
  addressOf(i: number): string | undefined {
    return ops?.addressOf(this, i);
  }

  /** A hit in words: its system's label, its ordinal positions, its strand ("pUC19 396 to 401, minus strand"). */
  describe(i: number): string {
    return ops?.describe(this, i) ?? '';
  }

  /**
   * The hits as a ref (one patch, one undo step: EC-162): `to` their addresses, `label` `Search: GAATTC (14)`. Refused,
   * never cut: an incomplete result, past `search.maxSaved` hits.
   */
  async save(options: { readonly id?: string; readonly label?: string } = {}): Promise<SearchOutcome> {
    return (await loadOps()).save(this, options);
  }

  /** @internal Whether the model was retired (its work's check). */
  get _gone(): boolean {
    return this.isRetired;
  }

  /** Ends a run in progress with the model. */
  protected override onRetire(): void {
    this._run?.abort();
  }
}
