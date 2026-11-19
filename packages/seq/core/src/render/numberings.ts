/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, type Column, type Issue, issueOf } from '@molstar/core/repr/tree';
import type { Coverage, DataEntry, RangesList, Want } from '@molstar/seqx-core/data';
import { labelsNumbering, type Numbering, type QueryModel, type QuerySystem } from '@molstar/seqx-core/query';
import type { NumberingEntry, StaticModel } from '@molstar/seqx-core/resolve';
import { DefaultRenderConfig } from './config.js';
import { RenderIssues } from './issues.js';
import type { StoreReader } from './layer-type.js';

// A source's numbering (_spec/README.md "Render", Numberings from a source): its labels are its rows' (a row a
// position: its `position` or its `start`; the label its `field`), read whole on the system at level 1, held once
// every tile of the system is in (labels need no order, so a half-read numbering would snap and print wrong), and kept
// while the source's coverage is the same (an evicted tile never makes it ask again). Bounded: a system past
// `maxNumberingPositions` is neither asked for nor held. Whoever reads it asks for its tiles: a ruler at layout, the
// viewer's addresses when a form names it.

/** What a source's numbering is now: its labels held, pending, or blank with why. */
export type SourceLabels =
  | { readonly status: 'held'; readonly numbering: Numbering; readonly issues: readonly Issue[] }
  | { readonly status: 'pending' }
  | { readonly status: 'blank'; readonly issue: Issue };

export interface SourceLabelsOptions {
  /** Positions it may cover at most (`RenderConfig.maxNumberingPositions`). */
  readonly maxPositions?: number;
  /** Asks for its tiles (a layer's `ctx.want`, an address model's wanter); without it only what is held is read. */
  readonly want?: (w: Want) => void;
}

/** What a source's numbering read last: the coverage and entries it was read from (compared in place), its result. */
interface Read {
  readonly coverage: Coverage | undefined;
  readonly entries: readonly (DataEntry | undefined)[];
  readonly result: SourceLabels;
}

const reads = new WeakMap<NumberingEntry, Read>();
const Pending: SourceLabels = { status: 'pending' };
const NoEntries: readonly (DataEntry | undefined)[] = [];

/** A cell's label: a category's text, a number's spelling; null for none. */
function cellLabel(c: Column, i: number): string | null {
  if (c instanceof CategoryColumn) return c.get(i);
  const v = (c as ArrayLike<number>)[i];
  return v === undefined || Number.isNaN(v) ? null : String(v);
}

function blank(
  entry: NumberingEntry,
  code: 'seq.render.numbering.limit' | 'seq.render.numbering.source',
  message: string,
): SourceLabels {
  return {
    status: 'blank',
    issue: issueOf(RenderIssues, code, message, { address: entry.address, hint: RenderIssues.codes[code].hint }),
  };
}

function rowsIssue(entry: NumberingEntry, message: string): Issue {
  return issueOf(RenderIssues, 'seq.render.numbering.rows', message, {
    address: entry.address,
    hint: RenderIssues.codes['seq.render.numbering.rows'].hint,
  });
}

/** Whether `read` was made from this coverage and these entries. */
function same(read: Read | undefined, coverage: Coverage | undefined, entries: readonly (DataEntry | undefined)[]) {
  return (
    read !== undefined &&
    read.coverage === coverage &&
    read.entries.length === entries.length &&
    read.entries.every((e, i) => e === entries[i])
  );
}

/**
 * A source's numbering from the store. Held labels are kept while the source's coverage is the same, whatever
 * becomes of its tiles; else pending until every tile of its system is in, blank when a tile failed, the source has
 * no system of the sequence's name or the system is past the bound. A held numbering's issues: rows past the end
 * (left out), rows for another count of positions than the sequence has.
 */
export function readSourceNumbering(
  model: StaticModel,
  entry: NumberingEntry,
  data: StoreReader,
  options: SourceLabelsOptions = {},
): SourceLabels {
  if (!entry.source || !entry.access) return Pending;
  const max = options.maxPositions ?? DefaultRenderConfig.maxNumberingPositions;
  const coverage = data.coverage(entry.source);
  const prev = reads.get(entry);
  if (prev && prev.coverage === coverage && prev.result.status === 'held') return prev.result;
  const keep = (entries: readonly (DataEntry | undefined)[], result: SourceLabels) => {
    reads.set(entry, { coverage, entries, result });
    return result;
  };
  const name = entry.access.system;
  const own = model.systems.get(entry.system)?.length;
  const want = (end: number) => options.want?.({ source: entry.source!, system: name, start: 0, end, level: 1 });
  if (own !== undefined && own > max) {
    if (same(prev, coverage, NoEntries)) return prev!.result;
    return keep(NoEntries, blank(entry, 'seq.render.numbering.limit', `The numbering covers ${own} positions`));
  }
  if (!coverage) {
    // not attached yet: asked for where its length is known, so it comes once the source opens
    if (own !== undefined) want(own);
    return Pending;
  }
  const sys = coverage.systems.find((s) => s.name === name || s.aliases?.includes(name));
  const length = own ?? sys?.length;
  if (!sys || length === undefined) {
    if (same(prev, coverage, NoEntries)) return prev!.result;
    return keep(
      NoEntries,
      blank(
        entry,
        'seq.render.numbering.source',
        sys ? `The source does not say how long ${name} is` : `The source has no system ${name}`,
      ),
    );
  }
  if (length > max) {
    if (same(prev, coverage, NoEntries)) return prev!.result;
    return keep(NoEntries, blank(entry, 'seq.render.numbering.limit', `The numbering covers ${length} positions`));
  }
  const w: Want = { source: entry.source, system: name, start: 0, end: length, level: 1 };
  want(length);
  const entries = data.keysOf(w).map((k) => data.get(k));
  if (same(prev, coverage, entries)) return prev!.result;
  const failed = entries.find((e) => e?.status === 'failed');
  if (failed?.status === 'failed')
    return keep(
      entries,
      blank(entry, 'seq.render.numbering.source', `A tile of the source failed: ${failed.issue.message}`),
    );
  if (!entries.every((e) => e?.status === 'loaded' || e?.status === 'empty')) return keep(entries, Pending);
  const labels: (string | null)[] = new Array(length).fill(null);
  // the positions a row names (a row with no label too: `null` is no label there, as inline)
  const rowed = new Uint8Array(length);
  const field = entry.field ?? 'label';
  let past = 0;
  let labelled = 0;
  for (const e of entries) {
    if (e?.status !== 'loaded') continue;
    const t = e.tables.ranges as RangesList | undefined;
    if (t?.form !== 'list') continue;
    const column = field === 'id' ? t.id : t.fields?.find((f) => f.name === field)?.column;
    for (let k = 0; k < t.start.length; k++) {
      const p = t.start[k];
      if (p < 0 || p >= length) {
        past++;
        continue;
      }
      // a row a position; the first row of a position keeps it
      if (rowed[p]) continue;
      rowed[p] = 1;
      labelled++;
      if (column) labels[p] = cellLabel(column, k);
    }
  }
  const issues: Issue[] = [];
  if (past) issues.push(rowsIssue(entry, `${past} rows of the source are past the sequence’s ${length} positions`));
  if (labelled !== length) issues.push(rowsIssue(entry, `The source has rows for ${labelled} of ${length} positions`));
  const numbering = labelsNumbering(entry.id, { length, get: (i) => labels[i] ?? null });
  return keep(entries, { status: 'held', numbering, issues });
}

/**
 * A query model whose systems also answer their sources' numberings (_spec/README.md "Render", Numberings from a
 * source): held labels as a numbering, labels not in yet as pending (a form naming them waits; `want` asks for their
 * tiles), and a blank one left out of the ids. A unit's molecule is wrapped too, so a molecule's source numbering
 * numbers its units. Everything else is `base`'s.
 */
export function withSourceNumberings(
  base: QueryModel,
  model: StaticModel,
  data: StoreReader,
  options: SourceLabelsOptions = {},
): QueryModel {
  if (!model.hasSourceNumberings) return base;
  const systems = new Map<string, QuerySystem | undefined>();
  const read = (address: string, id: string): SourceLabels | undefined => {
    const entry = model.numberingOf(address, id);
    return entry?.source ? readSourceNumbering(model, entry, data, options) : undefined;
  };
  const wrap = (s: QuerySystem): QuerySystem => {
    let out = systems.get(s.address);
    if (out) return out;
    const address = s.address;
    out = {
      address,
      length: s.length,
      ...(s.cyclic !== undefined ? { cyclic: s.cyclic } : {}),
      ...(s.parent ? { parent: { start: s.parent.start, system: wrap(s.parent.system) } } : {}),
      numberingIds: () =>
        s.numberingIds().filter((id) => {
          const r = read(address, id);
          return r === undefined || r.status !== 'blank';
        }),
      numbering: (id) => {
        const own = s.numbering(id);
        if (own) return own;
        const r = read(address, id);
        return r?.status === 'held' ? r.numbering : undefined;
      },
      numberingPending: (id) => {
        const r = read(address, id);
        return r ? r.status === 'pending' : (s.numberingPending?.(id) ?? false);
      },
    };
    systems.set(address, out);
    return out;
  };
  return {
    system(address) {
      if (systems.has(address)) return systems.get(address);
      const s = base.system(address);
      if (!s) systems.set(address, undefined);
      return s && wrap(s);
    },
    annotation: (address) => base.annotation(address),
    ...(base.node ? { node: (owner, key) => base.node!(owner, key) } : {}),
    ...(base.chain ? { chain: (owner, id) => base.chain!(owner, id) } : {}),
    ...(base.row ? { row: (rows, key) => base.row!(rows, key) } : {}),
    ...(base.rowsPending ? { rowsPending: (rows) => base.rowsPending!(rows) } : {}),
  };
}
