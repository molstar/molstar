/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Queries in free text (seqspec:seqql#queries-in-text): the [SeqQL]…[/SeqQL] blocks of a reply, each with its place in
// the whole text, so the reader's, binder's and checker's issues inside a block point into the reply. Tags first; in a
// text with none, ```seqql fences; with neither, the whole text if it reads as one form.
import { queryIssue } from './issue.js';
import { readQuery, type TextPosition } from './read.js';
import type { QueryIssue } from './types.js';

export interface QueryBlock {
  /** The query as written between the tags (or inside the fence, or the whole text). */
  readonly text: string;
  /** Where `text` starts in the whole text; `checkText(block.text, language, block.at)` places issues there. */
  readonly at: TextPosition;
  /** A `[SeqQL]` block, a ```seqql fence, or the whole text read as one query. */
  readonly via: 'tag' | 'fence' | 'text';
}

export interface QueryExtraction {
  /** The queries in order. */
  readonly blocks: readonly QueryBlock[];
  /** Tags that do not pair up, empty blocks, or no query at all (`seqql.reply.*`). */
  readonly issues: readonly QueryIssue[];
}

/** The canonical tags; matched in any case. */
export const SEQQL_OPEN = '[SeqQL]';
export const SEQQL_CLOSE = '[/SeqQL]';
const TAG = /\[(\/?)seqql\]/gi;
const FENCE = /^```[ \t]*seqql[ \t]*\r?\n([\s\S]*?)^```[ \t]*$/gim;
const WRAP = `wrap each query in ${SEQQL_OPEN}…${SEQQL_CLOSE}`;

/** The line and column (1-based) of an offset into a text. */
function positionOf(text: string, offset: number): TextPosition {
  let line = 1;
  let start = 0;
  for (let i = text.indexOf('\n'); i !== -1 && i < offset; i = text.indexOf('\n', i + 1)) {
    line++;
    start = i + 1;
  }
  return { line, column: offset - start + 1, offset };
}

/** Every query block of a text, in order, with where it starts; never throws. */
export function extractQueries(text: string): QueryExtraction {
  const blocks: QueryBlock[] = [];
  const issues: QueryIssue[] = [];
  const where = (o: number) => `${positionOf(text, o).line}:${positionOf(text, o).column}`;
  const say = (code: string, message: string, at: number | null, open: number | null = null) =>
    issues.push(
      queryIssue(code, message, {
        at: at === null ? null : positionOf(text, at),
        open: open === null ? null : positionOf(text, open),
      }),
    );
  const block = (body: string, start: number, via: QueryBlock['via'], opened: number) => {
    if (body.trim() === '')
      say('seqql.reply.empty', `an empty block: put one query in it, as ${SEQQL_OPEN}(all)${SEQQL_CLOSE}`, opened);
    else blocks.push({ text: body, at: positionOf(text, start), via });
  };

  let open: number | null = null;
  let tags = 0;
  for (const m of text.matchAll(TAG)) {
    tags++;
    if (!m[1]) {
      if (open !== null) {
        const message = `${SEQQL_OPEN} inside the block opened at ${where(open)}: close each block with ${SEQQL_CLOSE} before the next`;
        say('seqql.reply.nested', message, m.index, open);
      }
      open = m.index;
    } else if (open === null) {
      say('seqql.reply.unopened', `${SEQQL_CLOSE} closes no block: start each query with ${SEQQL_OPEN}`, m.index);
    } else {
      const start = open + SEQQL_OPEN.length;
      block(text.slice(start, m.index), start, 'tag', open);
      open = null;
    }
  }
  if (open !== null) {
    const message = `unclosed ${SEQQL_OPEN} opened at ${where(open)}: end the query with ${SEQQL_CLOSE}`;
    say('seqql.reply.unclosed', message, text.length, open);
  }
  if (tags > 0) return { blocks, issues };

  for (const m of text.matchAll(FENCE)) {
    const start = m.index + m[0].indexOf('\n') + 1;
    block(m[1] ?? '', start, 'fence', m.index);
  }
  if (blocks.length > 0 || issues.length > 0) return { blocks, issues };
  if (text.trim() !== '' && readQuery(text).issues.length === 0) {
    return { blocks: [{ text, at: positionOf(text, 0), via: 'text' }], issues };
  }
  say('seqql.reply.none', text.trim() === '' ? `the reply is empty: ${WRAP}` : `no query found: ${WRAP}`, null);
  return { blocks, issues };
}
