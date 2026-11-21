/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A text read whole as lines (_spec/README.md "Data", Genome runtimes' internals): apart from `lines.ts` (by region,
// over io's tabix), so a runtime reading files whole loads none of io's region readers. Not published.

/** One line's default bound: io's `DefaultTabixLimits.maxLineBytes` (16 MiB), not imported so tabix stays unloaded. */
export const MaxLineBytes = 16 * 1024 * 1024;

/** The lines of a file read whole: its header lines (`meta`, the first `skip`) and the rest, with 1-based numbers. */
export interface WholeLines {
  readonly header: LineSpans;
  readonly lines: LineSpans;
}

export interface LineSpans {
  readonly bytes: Uint8Array;
  readonly lineStart: Uint32Array;
  readonly lineEnd: Uint32Array;
  /** Each line's number in the file, 1-based (an issue's place). */
  readonly lineNumber: Uint32Array;
}

export interface LineOptions {
  /** A header line's first character (`#`); empty: none. */
  readonly meta: string;
  /** Lines skipped at the top as header. */
  readonly skip?: number;
  /** One line's bound (`MaxLineBytes`, 16 MiB). */
  readonly maxLineBytes?: number;
  /** Lines a file may hold (`limits.maxRows`). */
  readonly maxLines?: number;
}

export type LinesResult<T> =
  | { readonly ok: true; readonly value: T }
  | { readonly ok: false; readonly kind: 'too-large' | 'malformed'; readonly message: string };

/**
 * A text read whole split into lines (`\n`, a `\r` before it dropped; blank lines left out): the header lines at the top
 * (those starting with `meta`, and the first `skip`), then the rest, a `meta` line among them kept as a comment's
 * place is the format's. A line over `maxLineBytes` is `malformed`, more than `maxLines` lines `too-large`; the spans
 * are counted before they are allocated.
 */
export function wholeLines(bytes: Uint8Array, options: LineOptions): LinesResult<WholeLines> {
  const maxLine = options.maxLineBytes ?? MaxLineBytes;
  const maxLines = options.maxLines ?? Number.POSITIVE_INFINITY;
  const skip = options.skip ?? 0;
  const meta = options.meta ? options.meta.charCodeAt(0) : -1;
  /** Calls `line` for each line kept (the blank ones left out but a skipped one), the header ones first. */
  const scan = (line: (s: number, e: number, number: number, header: boolean) => boolean | undefined) => {
    let inHeader = true;
    for (let s = 0, number = 1; s < bytes.length; number++) {
      let e = bytes.indexOf(10, s);
      if (e < 0) e = bytes.length;
      const next = e + 1;
      if (e > s && bytes[e - 1] === 13) e--;
      if (inHeader && !(number <= skip || (e > s && bytes[s] === meta))) inHeader = false;
      if ((e > s || (inHeader && number <= skip)) && line(s, e, number, inHeader) === false) return;
      s = next;
    }
  };
  // counted first, so the spans are allocated once and a bound is refused before they are
  let headerCount = 0;
  let count = 0;
  let failure: LinesResult<WholeLines> | undefined;
  scan((s, e, number, header) => {
    if (e - s > maxLine) {
      failure = { ok: false, kind: 'malformed', message: `Line ${number} is over ${maxLine} bytes` };
      return false;
    }
    if (header) headerCount++;
    else if (count - headerCount >= maxLines) {
      failure = { ok: false, kind: 'too-large', message: `The file has more than ${maxLines} lines` };
      return false;
    }
    count++;
  });
  if (failure) return failure;
  const starts = new Uint32Array(count);
  const ends = new Uint32Array(count);
  const numbers = new Uint32Array(count);
  let n = 0;
  scan((s, e, number) => {
    starts[n] = s;
    ends[n] = e;
    numbers[n] = number;
    n++;
  });
  const spans = (from: number, to: number): LineSpans => ({
    bytes,
    lineStart: starts.subarray(from, to),
    lineEnd: ends.subarray(from, to),
    lineNumber: numbers.subarray(from, to),
  });
  return { ok: true, value: { header: spans(0, headerCount), lines: spans(headerCount, n) } };
}
