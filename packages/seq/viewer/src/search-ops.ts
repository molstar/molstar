/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The search model's work (_spec/README.md "Search"): a run, the steps, a hit's address and words, a save, the channels
// and following the spec. A lazy chunk the first run loads: a viewer that never searches loads none of it, and every
// step comes after a run.

import { type Issue, issueOf } from '@molstar/core/repr/tree';
import { DomHost } from '@molstar/ui/dom/host';
import { addRefPatch, RefsId } from '@molstar/seqx-core/lang';
import { AddressIssues, addressFor, EmptyKeys, type Location, location } from '@molstar/seqx-core/query';
import {
  compileFailed,
  HitList,
  IncompleteSearchCodes,
  loadMatcher,
  type SearchHits,
  SearchIssues,
  type SearchStrands,
  searchBounds,
} from '@molstar/seqx-core/search';
import { positionAt } from '@molstar/seqx-core/view';
import type { Emphasized } from '@molstar/seqx-viewer/models';
import { HiddenMessage, HiddenTrackMessage } from './models/address.js';
import type { SearchModel, SearchOutcome, SearchQuery, SearchResult } from './search.js';
import { searchTargets } from './search-targets.js';

const fmt = (n: number) => n.toLocaleString('en-US');
const plural = (n: number, one: string) => `${fmt(n)} ${one}${n === 1 ? '' : 's'}`;
const StrandWords: Readonly<Record<number, string>> = { 1: 'plus strand', [-1]: 'minus strand', 0: 'both strands' };

/**
 * Runs a search: the pattern compiled first (a compile error is the outcome's message and issue, and nothing changes),
 * a run in progress cancelled, the result replaced as batches arrive. Resolves when done, cancelled or failed.
 */
export async function runSearch(
  m: SearchModel,
  q: {
    readonly pattern: string;
    readonly matcher?: string;
    readonly strands?: SearchStrands;
    readonly within?: string;
  },
  id: number,
): Promise<SearchOutcome> {
  const { config } = m.deps;
  const name = q.matcher ?? config.search.matcher;
  const query: SearchQuery = {
    pattern: q.pattern,
    matcher: name,
    strands: q.strands ?? config.search.strands,
    ...(q.within !== undefined ? { within: q.within } : {}),
  };
  const entry = m.deps.matchers()[name];
  if (!entry) {
    const issue = issueOf(SearchIssues, 'seq.search.matcher.missing', `No search matcher ${name}`, {
      hint: `Registered: ${Object.keys(m.deps.matchers()).join(', ') || 'none'}`,
    });
    return { message: issue.message, issues: [issue] };
  }
  const matcher = await loadMatcher(entry);
  if (id !== m._runs || m._gone) return { message: '', issues: [] };
  const bounds = searchBounds(boundsOf(m));
  const compiled = matcher.compile(query.pattern, { strands: query.strands }, bounds);
  if (compileFailed(compiled)) return { message: compiled[0]?.message ?? 'Not searched', issues: compiled };

  // the targets: what the view draws, or what `within` names
  const within = withinOf(m, query.within);
  if (within && 'message' in within) return within;
  const model = m.deps.spec.model.value;
  const { targets, unknown } = searchTargets(model, m.deps.view.resolved.value, m.deps.store.store, within?.location);
  if (id !== m._runs || m._gone) return { message: '', issues: [] };
  m._run?.abort();
  const ac = new AbortController();
  m._run = ac;
  const held = new HitList(64);
  const result = (): SearchResult => ({
    hits: held,
    systems: targets.map((t) => t.system),
    lengths: targets.map((t) => t.length),
    labels: targets.map((t) => t.label),
  });
  const issues: Issue[] = unknown.map((a) =>
    issueOf(SearchIssues, 'seq.search.letters.failed', `${a} is not searched: its length is not known yet`, {
      address: a,
      hint: 'Search again once its source has opened',
    }),
  );
  m._query.set(query);
  m._status.set('running');
  m._current.set(-1);
  m._count.set(0);
  m._complete.set(false);
  m._stale.set(false);
  m._issues.set(issues);
  m._entries = targets.map((t) => model.systems.get(t.system));
  m._hits.set(result());
  paint(m);

  // progress said once, after a second, never per batch (a live region queues)
  let scanned = 0;
  const total = targets.reduce((n, t) => n + t.length, 0);
  const progress = DomHost.setTimeout(() => {
    if (m._run === ac && m._status.value === 'running')
      m.deps.notify(`Searching: ${Math.round((100 * scanned) / Math.max(1, total))} % read`);
  }, 1000);

  let count = 0;
  let failed: unknown;
  await new Promise<void>((resolve) => {
    const sub = matcher
      .run(compiled, targets, {
        signal: ac.signal,
        bounds,
        // a run never holds a frame: a turn of the event loop between batches
        yield: () => new Promise<void>((r) => setTimeout(r, 0)),
      })
      .subscribe({
        next: (b) => {
          if (m._run !== ac) return;
          scanned = b.scanned;
          const room = bounds.maxHits - held.length;
          held.append(b.hits, 0, Math.min(room, b.hits.length));
          count += b.hits.length + b.counted;
          if (b.issues.length) issues.push(...b.issues);
          m._count.set(count);
          if (b.issues.length) m._issues.set([...issues]);
          if (b.hits.length) {
            m._hits.set(result());
            paint(m);
          }
        },
        error: (e) => {
          failed = e;
          resolve();
        },
        complete: () => resolve(),
      });
    ac.signal.addEventListener('abort', () => {
      sub.unsubscribe();
      resolve();
    });
  });
  DomHost.clearTimeout(progress);
  if (m._run !== ac) return { message: '', issues: [] };
  m._run = undefined;
  if (ac.signal.aborted) {
    m._status.set('cancelled');
    return { message: '', issues };
  }
  if (failed !== undefined) {
    m._status.set('failed');
    const message = failed instanceof Error ? failed.message : 'The search failed';
    m.deps.notify(message);
    return { message, issues };
  }
  const incomplete = issues.some((i) => IncompleteSearchCodes.has(i.code));
  m._complete.set(!incomplete);
  m._status.set('done');
  const message = doneWords(query, count, held.length, targets.length);
  m.deps.notify(message);
  return { message, issues };
}

/** The words a finished run says. */
function doneWords(q: SearchQuery, count: number, held: number, systems: number): string {
  if (!count) return `No hits for ${q.pattern}`;
  if (held < count) return `${fmt(held)} of ${fmt(count)} hits held for ${q.pattern}`;
  return `${plural(count, 'hit')} for ${q.pattern} on ${plural(systems, 'sequence')}`;
}

/** The search's bounds from the configuration. */
function boundsOf(m: SearchModel) {
  const s = m.deps.config.search;
  return {
    maxPattern: s.maxPattern,
    maxRegexChars: s.maxRegexChars,
    maxRegexPositions: s.maxRegexPositions,
    maxRegexDepth: s.maxRegexDepth,
    maxHitLength: s.maxHitLength,
    maxLetters: s.maxLetters,
    maxSteps: s.maxSteps,
    maxHits: s.maxHits,
    chunkLetters: s.chunkLetters,
  };
}

/** What `within` names as a Location, or why not; undefined for none. */
function withinOf(m: SearchModel, within: string | undefined): { location: Location } | SearchOutcome | undefined {
  if (within === undefined) return undefined;
  if (within === 'selection') {
    const s = m.deps.selection.value.value.shown;
    return s.parts.length ? { location: s } : { message: 'Nothing selected', issues: [] };
  }
  const located = m.deps.addresses.locate(within);
  if (located.resolution.status !== 'resolved' || !located.location.parts.length)
    return { message: `Nothing to search in ${within}`, issues: located.resolution.issues };
  return { location: located.location };
}

export function step(m: SearchModel, dir: 1 | -1, focus: boolean): SearchOutcome {
  const n = m._hits.value.hits.length;
  if (!n) return { message: 'No hits', issues: [] };
  const cur = m._current.value;
  let i: number;
  let wrapped = false;
  if (cur < 0) i = fromView(m, dir);
  else {
    i = cur + dir;
    if (i >= n || i < 0) {
      i = dir === 1 ? 0 : n - 1;
      wrapped = true;
    }
  }
  return goTo(m, i, focus, wrapped ? (dir === 1 ? 'wrapped to the first hit' : 'wrapped to the last hit') : '');
}

/**
 * The hit a step starts at with none current: the first at or after the start of the window of the first leading
 * section showing hits (`next`), the last before its end (`previous`). Why from the view: on a genome "next" means
 * next from where the reader is.
 */
function fromView(m: SearchModel, dir: 1 | -1): number {
  const { hits, systems } = m._hits.value;
  const n = hits.length;
  const resolved = m.deps.view.resolved.value;
  for (const s of resolved?.sections ?? []) {
    if (m.deps.view.isFollowing(s.address)) continue;
    const t = systems.indexOf(s.axis.system);
    if (t < 0) continue;
    const v = m.deps.view.get(s.address);
    if (!v || v.pending) continue;
    const a = positionAt(s.axis, v.start, 'start');
    const b = positionAt(s.axis, v.end, 'end');
    if (a === undefined || b === undefined) continue;
    const lo = Math.min(a, b);
    const hi = Math.max(a, b);
    if (dir === 1) {
      for (let i = 0; i < n; i++)
        if (hits.target[i] > t || (hits.target[i] === t && hits.start[i] >= Math.floor(lo))) return i;
      return 0;
    }
    for (let i = n - 1; i >= 0; i--)
      if (hits.target[i] < t || (hits.target[i] === t && hits.end[i] <= Math.ceil(hi))) return i;
    return n - 1;
  }
  return dir === 1 ? 0 : n - 1;
}

/**
 * That hit becomes current and is revealed; the current one again clears it (EC-161); an index outside the hits is
 * refused, never clipped (EC-160).
 */
export function go(m: SearchModel, index: number, options: { readonly focus?: boolean } = {}): SearchOutcome {
  const n = m._hits.value.hits.length;
  if (!Number.isInteger(index) || index < 0 || index >= n)
    return { message: n ? `No hit ${index + 1}: there ${n === 1 ? 'is 1' : `are ${fmt(n)}`}` : 'No hits', issues: [] };
  if (index === m._current.value) {
    m._current.set(-1);
    paint(m);
    return { message: 'No current hit', issues: [] };
  }
  return goTo(m, index, options.focus ?? false, '');
}

function goTo(m: SearchModel, i: number, focus: boolean, wrapped: string): SearchOutcome {
  m._current.set(i);
  paint(m);
  const address = addressOf(m, i);
  // reveal through the address model: every leading section showing it pans (a circle by the shortest arc, EC-12), a
  // hit in a collapsed stretch lands at its break; one step of the view's history (EC-159: pans only)
  const o = address ? (focus ? m.deps.addresses.focus(address) : m.deps.addresses.reveal(address)) : undefined;
  const extra = o && (o.message === HiddenMessage || o.message === HiddenTrackMessage) ? `; ${o.message}` : '';
  const message = `Hit ${fmt(i + 1)} of ${fmt(m._count.value)}: ${describe(m, i)}${wrapped ? `, ${wrapped}` : ''}${extra}`;
  m.deps.notify(message);
  return { message, issues: o?.issues ?? [] };
}

/** A hit's Location: its range on its system, two ranges across a circle's origin (one hit). */
export function locationOf(m: SearchModel, i: number): Location | undefined {
  const { hits, systems, lengths } = m._hits.value;
  if (i < 0 || i >= hits.length) return undefined;
  const L = lengths[hits.target[i]];
  const s = hits.start[i];
  const e = hits.end[i];
  const ranges = e > L ? [0, e - L, s, L] : [s, e];
  return location({ parts: [{ system: systems[hits.target[i]], ranges }] });
}

/**
 * A hit's address: its system's, then its range in ordinal labels (copy link's numbering); a hit across a circle's
 * origin one range with `from > to` (EC-10). Made on demand, never stored a hit.
 */
export function addressOf(m: SearchModel, i: number): string | undefined {
  const loc = locationOf(m, i);
  if (!loc) return undefined;
  return addressFor(
    m.deps.spec.shown.value,
    m.deps.addresses.queryModel(),
    { type: 'location', location: loc },
    { ...m.deps.config.address, prefer: ['ids'], numbering: 'ordinal' },
  );
}

/** A hit in words: its system's label, its ordinal positions, its strand ("pUC19 396 to 401, minus strand"). */
export function describe(m: SearchModel, i: number): string {
  const { hits, labels, lengths } = m._hits.value;
  if (i < 0 || i >= hits.length) return '';
  const L = lengths[hits.target[i]];
  const from = hits.start[i] + 1;
  const to = hits.end[i] > L ? hits.end[i] - L : hits.end[i];
  const strand = hits.strand[i];
  const alphabet = m.deps.spec.model.value.systems.get(m._hits.value.systems[hits.target[i]])?.alphabet;
  const words = alphabet === 'protein' ? '' : `, ${StrandWords[strand]}`;
  return `${labels[hits.target[i]]} ${fmt(from)} to ${fmt(to)}${words}`;
}

/**
 * The hits as a ref (one patch, one undo step: EC-162): `to` their addresses, `label` `Search: GAATTC (14)`. Refused,
 * never cut: an incomplete result, past `search.maxSaved` hits.
 */
export async function save(
  m: SearchModel,
  options: { readonly id?: string; readonly label?: string } = {},
): Promise<SearchOutcome> {
  const q = m._query.value;
  const { hits } = m._hits.value;
  if (!q || !hits.length) return { message: 'No hits to save', issues: [] };
  if (!m._complete.value || m._status.value !== 'done') {
    const issue = issueOf(SearchIssues, 'seq.search.save.incomplete', 'The result is incomplete: not saved');
    return { message: issue.message, issues: [issue] };
  }
  // new letters since the run: the hits may not be where they were (design Decision 5, "results may be out of date")
  if (m._stale.value) {
    const issue = issueOf(SearchIssues, 'seq.search.save.stale', 'The result may be out of date: not saved', {
      hint: SearchIssues.codes['seq.search.save.stale'].hint!,
    });
    return { message: issue.message, issues: [issue] };
  }
  const max = m.deps.config.search.maxSaved;
  if (m._count.value > max) {
    const issue = issueOf(
      SearchIssues,
      'seq.search.save.limit',
      `${fmt(m._count.value)} hits: a saved set holds at most ${fmt(max)}`,
      { hint: SearchIssues.codes['seq.search.save.limit'].hint! },
    );
    return { message: issue.message, issues: [issue] };
  }
  const to: string[] = [];
  for (let i = 0; i < hits.length; i++) {
    const a = addressOf(m, i);
    if (a !== undefined && !to.includes(a)) to.push(a);
  }
  // the strand is not in a Location: said in the label when every hit is on one
  let strand = hits.strand[0];
  for (let i = 1; i < hits.length && strand !== 2; i++) if (hits.strand[i] !== strand) strand = 2;
  const { systems } = m._hits.value;
  const alphabetOf = (i: number) => m.deps.spec.model.value.systems.get(systems[hits.target[i]])?.alphabet;
  let nucleotide = false;
  for (let i = 0; i < hits.length && !nucleotide; i++) nucleotide = alphabetOf(i) !== 'protein';
  const on = nucleotide && strand !== 2 ? `, ${StrandWords[strand]}` : '';
  const label = options.label ?? `Search: ${q.pattern} (${fmt(m._count.value)}${on})`;
  const outcome = m.deps.spec.patch('Save the search', (root) => {
    const taken = new Set(root.children?.find((c) => c.id === RefsId)?.children?.map((c) => c.id));
    let n = 1;
    while (options.id === undefined && taken.has(`search-${n}`)) n++;
    return addRefPatch(root, { id: options.id ?? `search-${n}`, to, label });
  });
  const done = await outcome;
  return done.ok
    ? { message: `Saved ${plural(to.length, 'hit')} as ${label}`, issues: [] }
    : { message: done.issues[0]?.message ?? 'Not saved', issues: done.issues };
}

/** Both channels from the hits and the current one (one value each; the model's own). */
export function paint(m: SearchModel): void {
  const { highlight } = m.deps;
  const { hits, systems, lengths } = m._hits.value;
  if (!hits.length) {
    m._drawn = undefined;
    highlight.clear(m.channel);
    highlight.clear(m.currentChannel);
    return;
  }
  const by = new Map<string, number[]>();
  for (let i = 0; i < hits.length; i++) {
    const system = systems[hits.target[i]];
    let r = by.get(system);
    if (!r) by.set(system, (r = []));
    const L = lengths[hits.target[i]];
    if (hits.end[i] > L) r.push(hits.start[i], L, 0, hits.end[i] - L);
    else r.push(hits.start[i], hits.end[i]);
  }
  const loc = location({ parts: [...by].map(([system, r]) => ({ system, ranges: sortedRanges(r) })) });
  const all: Emphasized = { elements: EmptyKeys, location: loc, shown: loc, style: 'tint' };
  m._drawn = all;
  highlight.set(m.channel, all);
  const c = locationOf(m, m._current.value);
  if (c) highlight.set(m.currentChannel, { elements: EmptyKeys, location: c, shown: c, style: 'outline' });
  else highlight.clear(m.currentChannel);
}

/** Results follow the spec: a removed system's hits dropped (said), new letters for a searched one stale. */
export function followModel(m: SearchModel): void {
  const r = m._hits.value;
  if (!r.systems.length) return;
  const model = m.deps.spec.model.value;
  const gone = r.systems.filter((s) => !model.systems.has(s));
  const changed = r.systems.some((s, k) => model.systems.has(s) && model.systems.get(s) !== m._entries[k]);
  if (changed && !m._stale.value) {
    m._stale.set(true);
    m.deps.notify('Search results may be out of date');
  }
  if (!gone.length) return;
  const keep = r.systems.map((s) => !gone.includes(s));
  const kept = { ...r, hits: filterHits(r.hits, (t) => keep[t]) };
  const cur = m._current.value;
  m._current.set(cur >= 0 && keep[r.hits.target[cur]] ? countBefore(r.hits, cur, (t) => keep[t]) : -1);
  m._count.set(m._count.value - (r.hits.length - kept.hits.length));
  m._hits.set(kept);
  m._issues.set([
    ...m._issues.value,
    ...gone.map((a) =>
      issueOf(AddressIssues, 'seq.address.missing', `${a} is gone: its hits are dropped`, { address: a }),
    ),
  ]);
  paint(m);
}

/** Ranges as sorted `[start, end)` pairs. */
function sortedRanges(r: number[]): number[] {
  const pairs: [number, number][] = [];
  for (let i = 0; i < r.length; i += 2) pairs.push([r[i], r[i + 1]]);
  pairs.sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  return pairs.flat();
}

function filterHits(h: SearchHits, keep: (target: number) => boolean): SearchHits {
  const idx: number[] = [];
  for (let i = 0; i < h.length; i++) if (keep(h.target[i])) idx.push(i);
  return {
    length: idx.length,
    target: Int32Array.from(idx, (i) => h.target[i]),
    start: Int32Array.from(idx, (i) => h.start[i]),
    end: Int32Array.from(idx, (i) => h.end[i]),
    strand: Int8Array.from(idx, (i) => h.strand[i]),
  };
}

function countBefore(h: SearchHits, cur: number, keep: (target: number) => boolean): number {
  let n = 0;
  for (let i = 0; i < cur; i++) if (keep(h.target[i])) n++;
  return n;
}
