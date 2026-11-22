/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A GenBank feature's location (_spec/README.md "Data", GenBank): INSDC's grammar as a feature table writes it, into
// pieces of one item (0-based, end exclusive), its strand from `complement`, its partial ends. `a..b`, a single base
// `a`, `complement(…)`, `join(…)` and `order(…)` nested as INSDC writes them (`complement(join(…))`,
// `join(complement(…),…)`), at most `MaxDepth` deep and `maxPieces` pieces; `<` and `>` kept as partial ends; `a.b` read
// as its hull (said). A `join` over a circular record's origin (`join(2600..2686,1..120)`) is one item across it: the
// pieces after the origin a lap on, so its end passes the length (the cyclic system's wrap, as built for ranges).
// Refused, with the reason: a remote reference (`J00194.1:100..202`), a site between two bases (`a^b`), a position
// past the record's length, a malformed text. One pass over the text, never a throw.

/** How deep `complement`, `join` and `order` may nest. */
export const MaxDepth = 16;

/** The largest position a location may name (int32). */
const MaxPosition = 0x7fffffff;

export type Location =
  | {
      readonly ok: true;
      /** The pieces, `[start, end, …]`, 0-based and end exclusive, sorted, apart (touching ones merged). */
      readonly pieces: number[];
      readonly start: number;
      readonly end: number;
      /** −1 every piece complemented, 1 none, 0 mixed (or a protein's). */
      readonly strand: -1 | 0 | 1;
      /** A `<` (the lower end) or a `>` (the upper end) somewhere in it. */
      readonly partialStart: boolean;
      readonly partialEnd: boolean;
      /** An `a.b` read as its hull. */
      readonly approximate: boolean;
      /** A join across a circular record's origin. */
      readonly wraps: boolean;
    }
  | { readonly ok: false; readonly why: string };

/** A piece as read: its ends and its strand, in reading order. */
interface Leaf {
  start: number;
  end: number;
  strand: 1 | -1;
}

class Refused {
  constructor(readonly why: string) {}
}

const isDigit = (c: number) => c >= 48 && c <= 57;
const isLetter = (c: number) => (c >= 65 && c <= 90) || (c >= 97 && c <= 122) || c === 95;

/**
 * The location `text` (whitespace already dropped) on a record of `length` positions, `circular` or not: its pieces, or
 * why it is refused. `maxPieces` bounds the pieces before they are held.
 */
export function parseLocation(text: string, length: number, circular: boolean, maxPieces: number): Location {
  let i = 0;
  let count = 0;
  let partialStart = false;
  let partialEnd = false;
  let approximate = false;

  const word = (w: string) => {
    if (!text.startsWith(w, i)) return false;
    i += w.length;
    return true;
  };
  const expect = (c: string) => {
    if (text[i] !== c) throw new Refused(i < text.length ? `\`${text[i]}\` where \`${c}\` was expected` : 'cut short');
    i++;
  };
  const number = (): number => {
    const from = i;
    while (i < text.length && isDigit(text.charCodeAt(i)) && i - from <= 10) i++;
    if (i === from) throw new Refused(i < text.length ? `\`${text[i]}\` where a position was expected` : 'cut short');
    const v = Number(text.slice(from, i));
    if (i - from > 10 || v > MaxPosition) throw new Refused('a position past 2^31 − 1');
    if (v < 1) throw new Refused('a position before 1');
    return v;
  };
  const leaf = (out: Leaf[]) => {
    const c = text.charCodeAt(i);
    if (isLetter(c)) {
      // an accession before `:` is another record's; any other word is no location
      let k = i;
      while (
        k < text.length &&
        k - i < 64 &&
        (isLetter(text.charCodeAt(k)) || isDigit(text.charCodeAt(k)) || text[k] === '.')
      )
        k++;
      throw new Refused(text[k] === ':' ? 'a remote reference (another record’s positions)' : 'not a location');
    }
    if (word('<')) partialStart = true;
    const a = number();
    let b = a;
    if (word('..')) {
      if (word('>')) partialEnd = true;
      else if (word('<')) partialStart = true;
      b = number();
    } else if (word('.')) {
      // one base unknown within `a.b`: its hull
      approximate = true;
      b = number();
    } else if (text[i] === '^') throw new Refused('a site between two bases (`a^b`)');
    if (word('>')) partialEnd = true;
    if (b < a) throw new Refused('an end before its start');
    if (b > length) throw new Refused(`a position past the record’s length (${length})`);
    if (++count > maxPieces) throw new Refused(`more than ${maxPieces} pieces`);
    out.push({ start: a - 1, end: b, strand: 1 });
  };
  /** `depth` the operators around this location: at most `MaxDepth` nest. */
  const loc = (depth: number, out: Leaf[]) => {
    if (word('complement(')) {
      if (depth >= MaxDepth) throw new Refused(`nested deeper than ${MaxDepth}`);
      const inner: Leaf[] = [];
      loc(depth + 1, inner);
      expect(')');
      // the other strand, read the other way
      for (let k = inner.length - 1; k >= 0; k--) {
        const p = inner[k];
        p.strand = p.strand === 1 ? -1 : 1;
        out.push(p);
      }
      return;
    }
    if (word('join(') || word('order(')) {
      if (depth >= MaxDepth) throw new Refused(`nested deeper than ${MaxDepth}`);
      loc(depth + 1, out);
      while (word(',')) loc(depth + 1, out);
      expect(')');
      return;
    }
    leaf(out);
  };

  const leaves: Leaf[] = [];
  try {
    if (!text.length) throw new Refused('empty');
    loc(0, leaves);
    if (i !== text.length) throw new Refused(`\`${text.slice(i, i + 8)}\` after the location`);
  } catch (e) {
    if (e instanceof Refused) return { ok: false, why: e.why };
    throw e;
  }

  let plus = 0;
  for (const p of leaves) if (p.strand === 1) plus++;
  const strand: -1 | 0 | 1 = plus === leaves.length ? 1 : plus === 0 ? -1 : 0;
  // in the direction of the sequence's numbering: a minus strand's pieces were read from its far end
  const along = strand === -1 ? leaves.reverse() : leaves;
  let wraps = false;
  if (circular) {
    for (let k = 1; k < along.length; k++) {
      if (along[k].start + (wraps ? length : 0) >= along[k - 1].start) {
        if (wraps) along[k].start += length;
        if (wraps) along[k].end += length;
        continue;
      }
      if (wraps) return { ok: false, why: 'more than one lap of the circular record' };
      wraps = true;
      along[k].start += length;
      along[k].end += length;
    }
  }
  const pieces = merged(along);
  const start = pieces[0];
  const end = pieces[pieces.length - 1];
  if (end - start > length && circular) return { ok: false, why: 'more than one lap of the circular record' };
  // a wrap near 2^31 would pass the int32 positions a table holds
  if (end > MaxPosition) return { ok: false, why: 'a join over the origin past 2^31 − 1' };
  return { ok: true, pieces, start, end, strand, partialStart, partialEnd, approximate, wraps };
}

/** The leaves as pieces: sorted by start, overlapping or touching ones merged (the origin's seam one piece). */
function merged(leaves: readonly Leaf[]): number[] {
  const order = leaves.length > 1 ? [...leaves].sort((a, b) => a.start - b.start) : leaves;
  const out: number[] = [];
  for (const p of order) {
    if (out.length && p.start <= out[out.length - 1]) out[out.length - 1] = Math.max(out[out.length - 1], p.end);
    else out.push(p.start, p.end);
  }
  return out;
}
