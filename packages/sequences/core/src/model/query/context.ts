/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The query context (seqspec:seqql#the-context, seqspec:seqql#performance): the set, the inputs the host passes, the
// limits and the abort signal, the issues a run reports, the slots a test reads, and the builders a run reuses.

import { type QueryIssue, seqqlSeverity } from "@molstar/seq-core/query";
import type { SeqSet } from "../set.ts";
import { type BuildGuard, SelectionBuilder } from "./kernels.ts";
import type { SeqSelection } from "./selection.ts";

export interface QueryLimits {
  /** Positions visited by scans and tests. */
  readonly positions: number;
  /** Runs in any selection a run builds. */
  readonly runs: number;
  readonly groups: number;
  /** 0: none. Checked every 1,024 runs and per system. */
  readonly timeMs: number;
}
export const DEFAULT_LIMITS: QueryLimits = { positions: 50_000_000, runs: 1_000_000, groups: 1_000_000, timeMs: 0 };

export interface QueryOptions {
  readonly set: SeqSet;
  /** Named selections (`(input selection)`). */
  readonly inputs?: Readonly<Record<string, SeqSelection>>;
  readonly limits?: Partial<QueryLimits>;
  readonly signal?: AbortSignal;
  /** A host's own names for refs (a Context entity's id), tried before the set's keys; null: not one of its names. */
  readonly resolveRef?: (ref: string, set: SeqSet) => SeqSelection | null;
}

/**
 * The element a test reads (seqspec:seqql#elements-and-groups), one reused object: a system (kind, id, row), a
 * position, an item, a group.
 */
export interface QuerySlot {
  kind: number;
  id: number;
  row: number;
  index: number;
  annotation: number;
  item: number;
  group: number;
  selection: SeqSelection | null;
}

/** An error that ends a run; its issue is the run's. */
export class QueryAbort extends Error {
  readonly issue: QueryIssue;
  constructor(issue: QueryIssue) {
    super(issue.message);
    this.issue = issue;
  }
}

/** An issue at its code's registered severity (seqspec:seqql#errors). */
const issue = (code: string, path: string, message: string): QueryIssue => ({
  code,
  severity: seqqlSeverity(code),
  path,
  message,
  at: null,
  open: null,
  signature: null,
});

/**
 * Where a number property writes its value for a relation (a double returned from a call is boxed). One for every run
 * (runs are synchronous); read through a context field, it got the optimized tests thrown away at each collection
 * (QR2b's heap guards).
 */
export const F64 = new Float64Array(1);

/** Builders kept between runs, so their arrays are grown once (seqspec:seqql#the-selection). */
const POOL: SelectionBuilder[] = [];

export class QueryContext implements BuildGuard {
  readonly set: SeqSet;
  readonly inputs: Readonly<Record<string, SeqSelection>>;
  readonly limits: QueryLimits;
  readonly options: QueryOptions;
  readonly issues: QueryIssue[] = [];
  readonly slot: QuerySlot = {
    kind: -1,
    id: -1,
    row: -1,
    index: -1,
    annotation: -1,
    item: -1,
    group: -1,
    selection: null,
  };
  /** Positions visited so far by scans and tests. */
  positions = 0;
  /** Inside `each`, the group being run on (`(current)`); null outside. */
  current: SeqSelection | null = null;
  /** Set by `and` for a scan that is its direct operand with no `:in`, and taken by it. */
  scope: SeqSelection | null = null;
  /** seqspec:seqql#empty-results: the missing values tests met, the last one's path; the first node that gave nothing from something, the
   * selection it had last (its input), what the runtimes said of theirs, and the printer of a node by its path. */
  missed = 0;
  missedAt = "";
  nonEmpty = 0;
  emptyAt: string | null = null;
  emptyFrom: SeqSelection | null = null;
  lastNonEmpty: SeqSelection | null = null;
  why: Map<string, string> | null = null;
  say: (path: string) => string = (path) => path;
  private readonly started = Date.now();

  constructor(o: QueryOptions) {
    this.options = o;
    this.set = o.set;
    this.inputs = o.inputs ?? {};
    this.limits = { ...DEFAULT_LIMITS, ...o.limits };
  }
  get maxRuns() {
    return this.limits.runs;
  }
  get maxGroups() {
    return this.limits.groups;
  }

  warn(code: string, path: string, message: string): void {
    this.issues.push(issue(code, path, message));
  }
  fail(code: string, path: string, message: string): never {
    throw new QueryAbort(issue(code, path, message));
  }
  over(what: "runs" | "groups", n: number): never {
    const limit = what === "runs" ? this.limits.runs : this.limits.groups;
    const say = `this builds ${n.toLocaleString("en-US")} ${what} or more (limit ${limit.toLocaleString("en-US")})`;
    this.fail(`seqql.limit.${what}`, "", `${say}; narrow first: intersect with a sequence or a range`);
  }
  /** A test met a missing value. */
  miss(path: string): void {
    this.missed++;
    this.missedAt = path;
  }
  /** What a node that gave nothing knew, said if it is the one named. */
  because(path: string, detail: string): void {
    this.why ??= new Map();
    this.why.set(path, detail);
  }
  /** A binder's `:where` passed nothing of `n` elements: how many, and the value missing since `missed0`. */
  tested(path: string, n: number, what: string, missed0: number): void {
    const m = this.missed - missed0;
    const miss =
      m > 0
        ? `; ${this.say(this.missedAt)} was missing ${m === 1 ? "once" : `${m.toLocaleString("en-US")} times`}`
        : "";
    const none = n === 1 ? `the one ${what} did not pass` : `none of ${n.toLocaleString("en-US")} ${what}s passed`;
    if (n > 0) this.because(path, `${none} :where${miss}`);
  }
  /** Counts `n` positions a scan or test is about to visit, against the limit. */
  visit(n: number, path: string, what: string): void {
    this.positions += n;
    const max = this.limits.positions;
    if (this.positions > max) {
      const say = `this visits ${this.positions.toLocaleString("en-US")} positions in ${what}`;
      this.fail("seqql.limit.positions", path, `${say} (limit ${max.toLocaleString("en-US")}); narrow first`);
    }
  }
  /** The clock and the abort signal: every 1,024 runs and per system. */
  tick(): void {
    if (this.options.signal?.aborted) this.fail("seqql.aborted", "", "the query was stopped");
    const ms = this.limits.timeMs;
    if (ms > 0 && Date.now() - this.started > ms) this.fail("seqql.limit.time", "", `the query took over ${ms} ms`);
  }

  /** A builder from the pool every run shares (runs are synchronous); `finish` builds its selection and returns it. */
  builder(): SelectionBuilder {
    return (POOL.pop() ?? new SelectionBuilder()).reset(this);
  }
  finish(b: SelectionBuilder): SeqSelection {
    const s = b.build(this.set);
    if (POOL.length < 16) POOL.push(b.reset(null));
    return s;
  }
}
