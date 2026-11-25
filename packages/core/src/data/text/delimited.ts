/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * Delimited records, RFC 4180 (https://www.rfc-editor.org/rfc/rfc4180), over text, with plain string cells (written
 * from `@molstar/io`'s CSV parser as the guide; that reader keeps its CIF column model):
 * - a cell is the text between delimiters, spaces included (never trimmed);
 * - a cell starting with the quote is quoted: it may hold delimiters, line breaks and the quote doubled (`""`), read
 *   as one quote; only spaces and tabs may follow its closing quote;
 * - records end at CRLF, LF or CR; the last line break is optional; empty lines are skipped;
 * - a leading byte order mark (U+FEFF) is skipped;
 * - not in the RFC: lines starting with `comment`, when set, are skipped (off by default).
 * An unterminated quoted cell, or text after a closing quote, is an error with its line. The cells are as written:
 * what a record's count of cells means (a header, ragged rows) is the caller's.
 */

export interface DelimitedOptions {
  /** One character; `,` by default (`\t` for TSV). */
  readonly delimiter?: string;
  /** One character; `"` by default. */
  readonly quote?: string;
  /** One character starting a comment line; none by default (comments are not in RFC 4180). */
  readonly comment?: string;
}

export type DelimitedResult =
  | {
      readonly ok: true;
      /** Each record's cells, in the text's order. */
      readonly records: string[][];
      /** Each record's first line, 1-based. */
      readonly lines: number[];
    }
  | { readonly ok: false; readonly message: string; readonly line: number };

class DelimitedError {
  constructor(
    readonly message: string,
    readonly line: number,
  ) {}
}

function charCode(char: string, name: string): number {
  if (char.length !== 1) throw new Error(`The ${name} must be one character, got '${char}'.`);
  return char.charCodeAt(0);
}

const LF = 10;
const CR = 13;

/**
 * The records of `text`. Throws only on bad options (a delimiter, quote or comment not one character, two of them the
 * same, or one a line break); a text that does not read is `{ ok: false }` with its line.
 */
export function readDelimited(text: string, options: DelimitedOptions = {}): DelimitedResult {
  const delimiter = charCode(options.delimiter ?? ',', 'delimiter');
  const quote = charCode(options.quote ?? '"', 'quote');
  const comment = options.comment ? charCode(options.comment, 'comment') : -1;
  const breaks = (c: number) => c === LF || c === CR;
  if (quote === delimiter || breaks(quote) || breaks(delimiter))
    throw new Error('The quote and the delimiter must differ and not be a line break.');
  if (comment === quote || comment === delimiter || breaks(comment))
    throw new Error('The comment must differ from the quote and the delimiter and not be a line break.');

  const length = text.length;
  let position = text.charCodeAt(0) === 0xfeff ? 1 : 0;
  let line = 1;
  const records: string[][] = [];
  const lines: number[] = [];

  const eatLineBreak = () => {
    const c = text.charCodeAt(position);
    if (c === CR) {
      ++position;
      if (text.charCodeAt(position) === LF) ++position;
      ++line;
    } else if (c === LF) {
      ++position;
      ++line;
    }
  };

  const readQuoted = (): string => {
    const open = line;
    let at = position + 1;
    let from = at;
    let value = '';
    while (true) {
      if (at >= length) throw new DelimitedError('An unterminated quoted cell.', open);
      const c = text.charCodeAt(at);
      if (c === quote) {
        if (text.charCodeAt(at + 1) === quote) {
          value += text.slice(from, at + 1);
          at += 2;
          from = at;
          continue;
        }
        break;
      }
      if (c === LF) {
        if (text.charCodeAt(at - 1) !== CR) ++line;
      } else if (c === CR) ++line;
      ++at;
    }
    value += text.slice(from, at);
    ++at;
    // spaces and tabs, unless one is the delimiter
    while (at < length) {
      const c = text.charCodeAt(at);
      if ((c !== 32 && c !== 9) || c === delimiter) break;
      ++at;
    }
    if (at < length) {
      const c = text.charCodeAt(at);
      if (c !== delimiter && !breaks(c)) throw new DelimitedError('Text after a closing quote.', line);
    }
    position = at;
    return value;
  };

  try {
    while (position < length) {
      const c = text.charCodeAt(position);
      if (breaks(c)) {
        eatLineBreak();
        continue;
      }
      if (c === comment) {
        while (position < length && !breaks(text.charCodeAt(position))) ++position;
        continue;
      }
      const cells: string[] = [];
      lines.push(line);
      while (true) {
        if (text.charCodeAt(position) === quote) cells.push(readQuoted());
        else {
          const start = position;
          let at = position;
          // a local per character: the hot loop of a large file
          while (at < length) {
            const d = text.charCodeAt(at);
            if (d === delimiter || d === LF || d === CR) break;
            ++at;
          }
          position = at;
          cells.push(text.slice(start, at));
        }
        if (position < length && text.charCodeAt(position) === delimiter) {
          ++position;
          continue;
        }
        eatLineBreak();
        break;
      }
      records.push(cells);
    }
  } catch (e) {
    if (e instanceof DelimitedError) return { ok: false, message: e.message, line: e.line };
    throw e;
  }
  return { ok: true, records, lines };
}
