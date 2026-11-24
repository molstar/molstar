/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A matcher of kernels (_spec/README.md "Search", Matchers): each target's letters read in chunks through its reader,
// overlapping by the longest hit less one, folded (upper case, U as T on nucleotides) with the gaps left out, the kernel
// run on each strand (the minus strand the reverse complement of the letters, its hits mapped back), and one batch a
// chunk, sorted, so the batches are a prefix of the final order. Every bound checked before the work it sizes.

import { Observable } from '@molstar/core/reactive';
import type { Issue } from '@molstar/core/repr/tree';
import { issueOf } from '@molstar/core/repr/tree';
import { sort } from '@molstar/core/data/util/sort';
import { complementCode } from '@molstar/seqx-core/resolve/codon';
import { HitList, NoHits } from './hits.js';
import { SearchIssues } from './issues.js';
import {
  EdgeEnd,
  EdgeStart,
  type SearchBatch,
  type SearchBounds,
  type SearchCompiled,
  type SearchKernel,
  type SearchMatcher,
  type SearchOptions,
  type SearchRunContext,
  type SearchTarget,
  type StepBudget,
} from './types.js';

/** What a matcher of kernels is made of: its compile, and a kernel for each target (or why the target is skipped). */
export interface KernelMatcherSpec<C extends SearchCompiled> {
  readonly name: string;
  readonly title: string;
  readonly alphabets?: SearchMatcher['alphabets'];
  compile(pattern: string, options: SearchOptions, bounds: SearchBounds): C | readonly Issue[];
  /** The kernel for a target: it sees the target's letters folded, gaps left out. An issue skips the target. */
  kernelOf(compiled: C, target: SearchTarget, bounds: SearchBounds): SearchKernel | Issue;
}

/** Folds of the stored letters by code: upper case, U as T on nucleotides; 0 for a gap (`-`, `.`), left out. */
function foldTable(nucleotide: boolean): Uint8Array {
  const t = new Uint8Array(256);
  for (let c = 1; c < 256; c++) t[c] = c >= 97 && c <= 122 ? c - 32 : c;
  t[45] = 0;
  t[46] = 0;
  if (nucleotide) t[85] = 84;
  return t;
}
const NucleotideFold = foldTable(true);
const ProteinFold = foldTable(false);

/** A chunk's buffers, grown and reused: the folded letters, their positions, the reverse complement, the hits. */
class Buffers {
  letters = new Uint8Array(0);
  at = new Int32Array(0);
  minus = new Uint8Array(0);
  readonly hits = new HitList(64);
  reserve(n: number): void {
    if (this.letters.length >= n) return;
    const c = Math.max(n, this.letters.length * 2);
    const letters = new Uint8Array(c);
    letters.set(this.letters);
    const at = new Int32Array(c);
    at.set(this.at);
    this.letters = letters;
    this.at = at;
    this.minus = new Uint8Array(c);
  }
}

/** The strand's order in a sort: `1`, `0`, `-1`. */
const strandRank = (s: number) => (s === 1 ? 0 : s === 0 ? 1 : 2);

function compareHits(h: HitList, i: number, j: number): number {
  return h.start[i] - h.start[j] || h.end[i] - h.end[j] || strandRank(h.strand[i]) - strandRank(h.strand[j]) || 0;
}

function swapHits(h: HitList, i: number, j: number): void {
  let t = h.start[i];
  h.start[i] = h.start[j];
  h.start[j] = t;
  t = h.end[i];
  h.end[i] = h.end[j];
  h.end[j] = t;
  t = h.strand[i];
  h.strand[i] = h.strand[j];
  h.strand[j] = t;
}

/** A kernel's answer or a skip. */
const isKernel = (k: SearchKernel | Issue): k is SearchKernel => typeof (k as SearchKernel).find === 'function';

/**
 * A matcher of kernels: an application's in-process matcher writes a kernel and nothing else (one in a worker or a
 * server writes `run`). The built-in matchers are made so.
 */
export function kernelMatcher<C extends SearchCompiled>(spec: KernelMatcherSpec<C>): SearchMatcher<C> {
  return {
    name: spec.name,
    title: spec.title,
    ...(spec.alphabets ? { alphabets: spec.alphabets } : {}),
    compile: (pattern, options, bounds) => spec.compile(pattern, options, bounds),
    run: (compiled, targets, ctx) =>
      new Observable<SearchBatch>((sub) => {
        const ended = new AbortController();
        const abort = () => ended.abort();
        if (ctx.signal.aborted) ended.abort();
        else ctx.signal.addEventListener('abort', abort);
        runKernels(spec, compiled, targets, ctx, ended.signal, (b) => sub.next(b)).then(
          () => sub.complete(),
          // a cancelled run ends with no error
          (e) => (ended.signal.aborted ? sub.complete() : sub.error(e)),
        );
        return () => {
          ctx.signal.removeEventListener('abort', abort);
          ended.abort();
        };
      }),
  };
}

/** The letters of `[from, to)` of a target, positions past its length read from its start (a circle). */
async function readUnrolled(target: SearchTarget, from: number, to: number, signal: AbortSignal): Promise<Uint8Array> {
  const L = target.length;
  if (to <= L) return target.letters.read(from, to, signal);
  const parts: Uint8Array[] = [];
  for (let p = from; p < to; ) {
    const lap = Math.floor(p / L) * L;
    const end = Math.min(to, lap + L);
    parts.push(await target.letters.read(p - lap, end - lap, signal));
    p = end;
  }
  const out = new Uint8Array(to - from);
  let k = 0;
  for (const part of parts) {
    out.set(part.subarray(0, Math.min(part.length, out.length - k)), k);
    k += part.length;
  }
  return out;
}

async function runKernels<C extends SearchCompiled>(
  spec: KernelMatcherSpec<C>,
  compiled: C,
  targets: readonly SearchTarget[],
  ctx: SearchRunContext,
  signal: AbortSignal,
  emit: (b: SearchBatch) => void,
): Promise<void> {
  const bounds = ctx.bounds;
  const budget: StepBudget = { steps: 0, limit: bounds.maxSteps };
  const buf = new Buffers();
  let letters = 0;
  let scanned = 0;
  let held = 0;
  let saidHits = false;
  const skip = (issue: Issue) => emit({ hits: NoHits, counted: 0, scanned, issues: [issue] });

  for (let t = 0; t < targets.length; t++) {
    if (signal.aborted) return;
    const target = targets[t];
    const L = target.length;
    if (L <= 0) continue;
    const kernel = spec.kernelOf(compiled, target, bounds);
    if (!isKernel(kernel)) {
      skip(kernel);
      continue;
    }
    const regions = target.within ?? Int32Array.of(0, L);
    let size = 0;
    for (let r = 0; r < regions.length; r += 2)
      size += Math.max(0, Math.min(L, regions[r + 1]) - Math.max(0, regions[r]));
    // the letters a run reads, checked before the target is read: the one that would pass the bound is skipped
    if (letters + size > bounds.maxLetters) {
      skip(
        issueOf(
          SearchIssues,
          'seq.search.letters.limit',
          `${target.system} is not searched: past the letters a search reads`,
          {
            address: target.system,
            hint: SearchIssues.codes['seq.search.letters.limit'].hint!,
          },
        ),
      );
      continue;
    }
    letters += size;
    const cyclic = target.cyclic && !target.within;
    const nucleotide = target.alphabet !== 'protein';
    const minus = compiled.strands === 'both' && nucleotide;
    const fold = nucleotide ? NucleotideFold : ProteinFold;
    const over = Math.max(0, (cyclic ? Math.min(kernel.maxLength, L) : kernel.maxLength) - 1);
    const strand = nucleotide ? 1 : 0;

    regions: for (let r = 0; r < regions.length; r += 2) {
      const rs = Math.max(0, regions[r]);
      const re = Math.min(L, regions[r + 1]);
      for (let a = rs; a < re; ) {
        const b = Math.min(re, a + bounds.chunkLetters);
        if (signal.aborted) return;
        // the chunk and its overlap, read on until the overlap holds `over` letters past its gaps (or the end)
        const hardEnd = cyclic ? b + over : re;
        let read = a;
        let n = 0;
        let starts = 0;
        let after = 0;
        let want = Math.min(hardEnd, b + over);
        try {
          while (read < want) {
            const raw = await readUnrolled(target, read, want, signal);
            if (signal.aborted) return;
            buf.reserve(n + raw.length);
            const lettersBuf = buf.letters;
            const at = buf.at;
            for (let k = 0; k < raw.length; k++) {
              const f = fold[raw[k]];
              if (f === 0) continue;
              const p = read + k;
              lettersBuf[n] = f;
              at[n++] = p;
              if (p < b) starts++;
              else after++;
            }
            read = want;
            if (after < over && read < hardEnd) want = Math.min(hardEnd, read + Math.max(over - after, 1024));
          }
        } catch (e) {
          if (signal.aborted) return;
          skip(
            issueOf(SearchIssues, 'seq.search.letters.failed', `The letters of ${target.system} could not be read`, {
              address: target.system,
              ...(e instanceof Error ? { hint: e.message } : {}),
            }),
          );
          break regions;
        }
        // the word steps, charged before the work: past the budget the run ends, incomplete
        const cost = n * kernel.words * (minus ? 2 : 1);
        if (budget.steps + cost > budget.limit) {
          skip(issueOf(SearchIssues, 'seq.search.steps.limit', 'The search ended at the work it may take'));
          return;
        }
        budget.steps += cost;

        const hits = buf.hits;
        hits.clear();
        const map = buf.at;
        const edges = (a === 0 ? EdgeStart : 0) | (!cyclic && read >= L ? EdgeEnd : 0);
        let within = kernel.find(buf.letters, n, 0, starts, edges, budget, (s, e) => {
          const start = map[s];
          const end = map[e - 1] + 1;
          if (end - start <= L) hits.push(t, start, end, strand);
        });
        if (within && minus) {
          const m = buf.minus;
          const l = buf.letters;
          for (let k = 0; k < n; k++) m[k] = complementCode(l[n - 1 - k]);
          const back = ((edges & EdgeEnd) !== 0 ? EdgeStart : 0) | ((edges & EdgeStart) !== 0 ? EdgeEnd : 0);
          within = kernel.find(m, n, 0, n, back, budget, (s, e) => {
            const start = map[n - e];
            const end = map[n - s - 1] + 1;
            if (start < b && end - start <= L) hits.push(t, start, end, -1);
          });
        }
        const batch = new HitList(Math.min(hits.length, 1024) || 1);
        let counted = 0;
        const issues: Issue[] = [];
        if (hits.length) {
          sort(hits, 0, hits.length, compareHits, swapHits);
          for (let i = 0; i < hits.length; i++) {
            let s = hits.strand[i];
            // a site the same on both strands is one hit
            if (s === 1 && i + 1 < hits.length && hits.strand[i + 1] === -1) {
              if (hits.start[i + 1] === hits.start[i] && hits.end[i + 1] === hits.end[i]) {
                s = 0;
                i++;
              }
            }
            if (held < bounds.maxHits) {
              batch.push(t, hits.start[i], hits.end[i], s);
              held++;
            } else counted++;
          }
        }
        if (counted && !saidHits) {
          saidHits = true;
          issues.push(
            issueOf(
              SearchIssues,
              'seq.search.hits.limit',
              `More hits than the ${bounds.maxHits.toLocaleString('en-US')} held`,
            ),
          );
        }
        scanned += b - a;
        emit({ hits: batch.length ? batch.trimmed() : NoHits, counted, scanned, issues });
        if (!within) {
          skip(issueOf(SearchIssues, 'seq.search.steps.limit', 'The search ended at the work it may take'));
          return;
        }
        await ctx.yield?.();
        a = b;
      }
    }
  }
}
