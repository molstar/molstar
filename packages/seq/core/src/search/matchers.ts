/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The built-in matchers (_spec/README.md "Search", Built-in matchers): `exact`, `degenerate` (IUPAC, the default) and
// `regex`, each a matcher of kernels over a system's stored letters, case-insensitive, U as T on nucleotides. Degenerate
// matching is subset matching: a sequence's ambiguity letter is a hit only where every base it may be is allowed. The
// system's alphabet decides what a pattern may hold; a target that cannot is skipped, said.

import type { Issue } from '@molstar/core/repr/tree';
import { issueOf } from '@molstar/core/repr/tree';
import type { SequenceAlphabet } from '@molstar/core/util/sequence-constants';
import { Iupac } from '@molstar/seqx-core/resolve/codon';
import { SearchIssues } from './issues.js';
import { kernelMatcher } from './kernel-matcher.js';
import { classKernel, classMasks, exactKernel } from './kernels.js';
import { glushkov, parseRegex, type RegexParsed, type RegexTables, regexKernel, regexTables } from './regex.js';
import type { SearchBounds, SearchCompiled, SearchMatcher, SearchOptions, SearchTarget } from './types.js';

const A = 65;
const Z = 90;
const Star = 42;
const isNucleotide = (a: SequenceAlphabet) => a !== 'protein';

/** The letters each alphabet holds in a pattern: IUPAC's on nucleotides, A to Z and `*` on proteins. */
const NucleotideLetters = new Set(Object.keys(Iupac).map((k) => k.charCodeAt(0)));
const ProteinLetters = new Set([...Array.from({ length: 26 }, (_, i) => A + i), Star]);
/** Under `degenerate` a protein pattern's U is refused: U is T on nucleotides, never a residue's set. */
const ProteinDegenerate = new Set([...ProteinLetters].filter((c) => c !== 85));

/** The bases or residues a letter may be, as bits: a sequence letter's set must lie inside its pattern letter's. */
const Other = 1 << 30;
function nucleotideSets(): Int32Array {
  const t = new Int32Array(256).fill(Other);
  const bit: Record<string, number> = { A: 1, C: 2, G: 4, T: 8 };
  for (const [k, bases] of Object.entries(Iupac)) {
    let s = 0;
    for (const b of bases) s |= bit[b];
    t[k.charCodeAt(0)] = s;
  }
  return t;
}
function proteinSets(): Int32Array {
  const t = new Int32Array(256).fill(Other);
  for (let c = A; c <= Z; c++) t[c] = 1 << (c - A);
  const of = (s: string) => [...s].reduce((m, x) => m | (1 << (x.charCodeAt(0) - A)), 0);
  t[66] = of('DN'); // B
  t[90] = of('EQ'); // Z
  t[74] = of('IL'); // J
  t[88] = (1 << 26) - 1; // X: any residue
  t[Star] = 1 << 26;
  return t;
}
const NucleotideSets = nucleotideSets();
const ProteinSets = proteinSets();

/** A pattern's codes, upper case, whitespace left out; or its issue. */
function patternCodes(pattern: string, bounds: SearchBounds): Uint8Array | readonly Issue[] {
  const text = pattern.replace(/\s+/g, '').toUpperCase();
  if (!text.length) return [issueOf(SearchIssues, 'seq.search.pattern.empty', 'The pattern is empty')];
  if (text.length > bounds.maxPattern)
    return [
      issueOf(
        SearchIssues,
        'seq.search.pattern.size',
        `The pattern has ${text.length} letters; a pattern has at most ${bounds.maxPattern}`,
      ),
    ];
  const out = new Uint8Array(text.length);
  for (let i = 0; i < text.length; i++) {
    const c = text.charCodeAt(i);
    if ((c < A || c > Z) && c !== Star)
      return [
        issueOf(SearchIssues, 'seq.search.pattern.syntax', 'The pattern holds a character that is no letter', {
          offset: i,
        }),
      ];
    out[i] = c;
  }
  return out;
}

/** The alphabet check: the first letter a target's alphabet does not hold, as its issue; undefined when it holds all. */
function alphabetIssue(codes: ArrayLike<number>, target: SearchTarget, holds: ReadonlySet<number>): Issue | undefined {
  for (let i = 0; i < codes.length; i++)
    if (!holds.has(codes[i]))
      return issueOf(
        SearchIssues,
        'seq.search.alphabet',
        `${target.system} is not searched: its alphabet (${target.alphabet}) has no ${String.fromCharCode(codes[i])}`,
        { address: target.system },
      );
  return undefined;
}

interface CodesCompiled extends SearchCompiled {
  readonly codes: Uint8Array;
}

const strandsOf = (o: SearchOptions) => o.strands ?? 'both';

/** The pattern's codes on a nucleotide target: U as T (the target's letters are folded so). */
const asNucleotides = (codes: Uint8Array) => codes.map((c) => (c === 85 ? 84 : c));

export const ExactMatcher: SearchMatcher<CodesCompiled> = kernelMatcher<CodesCompiled>({
  name: 'exact',
  title: 'Exact',
  compile(pattern, options, bounds) {
    const codes = patternCodes(pattern, bounds);
    return codes instanceof Uint8Array ? { pattern, strands: strandsOf(options), codes } : codes;
  },
  kernelOf(c, target) {
    const nuc = isNucleotide(target.alphabet);
    const issue = alphabetIssue(c.codes, target, nuc ? NucleotideLetters : ProteinLetters);
    return issue ?? exactKernel(nuc ? asNucleotides(c.codes) : c.codes);
  },
});

export const DegenerateMatcher: SearchMatcher<CodesCompiled> = kernelMatcher<CodesCompiled>({
  name: 'degenerate',
  title: 'IUPAC',
  compile(pattern, options, bounds) {
    const codes = patternCodes(pattern, bounds);
    return codes instanceof Uint8Array ? { pattern, strands: strandsOf(options), codes } : codes;
  },
  kernelOf(c, target) {
    const nuc = isNucleotide(target.alphabet);
    const issue = alphabetIssue(c.codes, target, nuc ? NucleotideLetters : ProteinDegenerate);
    if (issue) return issue;
    const sets = nuc ? NucleotideSets : ProteinSets;
    const m = c.codes.length;
    // position i takes a sequence code whose set lies inside the pattern letter's
    return classKernel(
      classMasks(m, (i, code) => (sets[code] & ~sets[c.codes[i]]) === 0),
      m,
    );
  },
});

interface RegexCompiled extends SearchCompiled {
  readonly parsed: RegexParsed;
  readonly tables: RegexTables;
}

export const RegexMatcher: SearchMatcher<RegexCompiled> = kernelMatcher<RegexCompiled>({
  name: 'regex',
  title: 'Regular expression',
  compile(pattern, options, bounds) {
    const parsed = parseRegex(pattern, bounds);
    if (Array.isArray(parsed)) return parsed;
    const p = parsed as RegexParsed;
    // the tables only once every bound held
    return { pattern, strands: strandsOf(options), parsed: p, tables: regexTables(glushkov(p)) };
  },
  kernelOf(c, target, bounds) {
    const { parsed } = c;
    if ((parsed.anchorStart || parsed.anchorEnd) && target.cyclic && !target.within)
      return issueOf(
        SearchIssues,
        'seq.search.pattern.anchor',
        `${target.system} is cyclic: it has no ends to anchor at`,
        {
          address: target.system,
        },
      );
    const nuc = isNucleotide(target.alphabet);
    const issue = alphabetIssue(parsed.letters, target, nuc ? NucleotideLetters : ProteinLetters);
    if (issue) return issue;
    const longest = parsed.maxMatch < 0 ? bounds.maxHitLength : Math.min(parsed.maxMatch, bounds.maxHitLength);
    return regexKernel(
      c.tables,
      parsed,
      longest,
      nuc ? (set, code) => set[code] === 1 || (code === 84 && set[85] === 1) : (set, code) => set[code] === 1,
    );
  },
});

/** The built-in matchers by name. */
export const BuiltInSearchMatchers: Readonly<Record<'exact' | 'degenerate' | 'regex', SearchMatcher<any>>> = {
  exact: ExactMatcher,
  degenerate: DegenerateMatcher,
  regex: RegexMatcher,
};
