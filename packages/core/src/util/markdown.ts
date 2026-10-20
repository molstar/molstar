/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { decodeHtmlEntities } from './html-entities.js';

// A small Markdown tokenizer for labels, tooltips and short descriptions: a typed tree, never HTML, so whatever renders
// it builds elements from text (a document's Markdown cannot inject markup). The approach follows snarkdown (Jason
// Miller, MIT): a line pass for the blocks and one scanner for the inline marks. What it reads: ATX headings, paragraphs
// (a single newline is a space; two trailing spaces or a backslash make a break), bullet and ordered lists (nested by
// indentation), fenced code, rules; inline: code spans, strong and emphasis (stars or underscores, underscores not inside
// words), links (balanced parentheses in the destination) and autolinks, image alt text, backslash escapes, the common
// entities. HTML is dropped (as react-markdown's `skipHtml`). Not read: tables, block quotes, setext headings, reference
// links.
//
// The input is untrusted (a loaded document's tooltips, parsed on each hover), so the parser is linear in its input: no
// regular expression with ambiguous adjacent quantifiers, an inline scanner whose closer searches are cached forward
// scans (each mark's search resumes where the last one stopped, never rescans from every opener), and nesting capped
// at `MaxInlineDepth` (deeper strong, emphasis and links stay plain text), which also bounds the recursion.

export type MarkdownInline =
  | { readonly kind: 'text'; readonly text: string }
  | { readonly kind: 'strong' | 'em'; readonly children: readonly MarkdownInline[] }
  | { readonly kind: 'code'; readonly text: string }
  | { readonly kind: 'link'; readonly href: string; readonly children: readonly MarkdownInline[] }
  | { readonly kind: 'break' };

export type MarkdownBlock =
  | { readonly kind: 'heading'; readonly level: number; readonly children: readonly MarkdownInline[] }
  | { readonly kind: 'paragraph'; readonly children: readonly MarkdownInline[] }
  | {
      readonly kind: 'list';
      readonly ordered: boolean;
      readonly start: number;
      readonly items: readonly (readonly MarkdownBlock[])[];
    }
  | { readonly kind: 'code'; readonly lang: string; readonly text: string }
  | { readonly kind: 'rule' };

const Escapable = '\\`*_{}[]()#+-.!<>|~';

const decode = decodeHtmlEntities;

/** A link destination's text: backslash escapes and entities read in one pass. */
function decodeHref(s: string): string {
  return s.replace(
    /\\([\\`*_{}[\]()#+\-.!<>|~])|&(?:#x[0-9a-f]+|#\d+|\w+);/gi,
    (m, escaped: string | undefined) => escaped ?? decodeHtmlEntities(m),
  );
}

/**
 * Whether a link's target is one a renderer may follow: http(s), mailto, or relative (no scheme). The scheme is read as
 * a browser's URL parser reads it: leading and trailing C0 controls and spaces stripped, tabs and newlines removed
 * anywhere (`\x01javascript:` and `java\tscript:` are `javascript:`).
 */
export function isSafeMarkdownHref(href: string): boolean {
  let start = 0;
  let end = href.length;
  while (start < end && href.charCodeAt(start) <= 0x20) start++;
  while (end > start && href.charCodeAt(end - 1) <= 0x20) end--;
  const url = href.slice(start, end).replace(/[\t\n\r]/g, '');
  const scheme = /^([a-z][a-z0-9+.-]*):/i.exec(url);
  return !scheme || ['http', 'https', 'mailto'].includes(scheme[1].toLowerCase());
}

/** The blocks of a Markdown text. */
export function parseMarkdown(source: string): MarkdownBlock[] {
  const lines = source.replace(/\r\n?/g, '\n').split('\n');
  const blocks: MarkdownBlock[] = [];
  let i = 0;
  let paragraph: string[] = [];
  const endParagraph = () => {
    if (paragraph.length) blocks.push({ kind: 'paragraph', children: parseMarkdownInline(paragraph.join('\n')) });
    paragraph = [];
  };
  while (i < lines.length) {
    const line = lines[i];
    const trimmed = line.trim();
    const fence = /^(`{3,}|~{3,})\s*([\w+-]*)/.exec(trimmed);
    if (fence) {
      endParagraph();
      const body: string[] = [];
      i++;
      while (i < lines.length && !lines[i].trim().startsWith(fence[1])) body.push(lines[i++]);
      i++;
      blocks.push({ kind: 'code', lang: fence[2], text: body.join('\n') });
      continue;
    }
    if (!trimmed) {
      endParagraph();
      i++;
      continue;
    }
    const heading = atxHeading(trimmed);
    if (heading) {
      endParagraph();
      blocks.push({ kind: 'heading', level: heading.level, children: parseMarkdownInline(heading.text) });
      i++;
      continue;
    }
    if (/^([-*_])(\s*\1){2,}$/.test(trimmed)) {
      endParagraph();
      blocks.push({ kind: 'rule' });
      i++;
      continue;
    }
    if (listItem(line)) {
      endParagraph();
      const [list, next] = parseList(lines, i, indentOf(line));
      blocks.push(list);
      i = next;
      continue;
    }
    paragraph.push(line);
    i++;
  }
  endParagraph();
  return blocks;
}

/**
 * An ATX heading (`## Title ##`) of a trimmed line: its level and its text, the closing `#`s and spaces stripped by
 * hand (a regular expression with `\s*#*\s*$` after a lazy group is cubic on a run of spaces).
 */
function atxHeading(trimmed: string): { level: number; text: string } | undefined {
  let level = 0;
  while (level < trimmed.length && trimmed[level] === '#') level++;
  if (level < 1 || level > 6 || !/\s/.test(trimmed[level] ?? '')) return;
  let end = trimmed.length;
  while (end > level && trimmed[end - 1] === '#') end--;
  return { level, text: trimmed.slice(level, end).trim() };
}

// `[\s\S]*` rather than `.*`: `.` stops at U+2028, and the failed `$` would backtrack over every split of the spaces.
const ListItem = /^(\s*)([-*+]|(\d+)[.)])\s+([\s\S]*)$/;

function listItem(line: string) {
  return ListItem.exec(line);
}

function indentOf(line: string): number {
  return /^\s*/.exec(line)![0].replace(/\t/g, '    ').length;
}

/** A list from line `i` at `indent`: its items, deeper items as nested lists, continuation lines into the item's text. */
function parseList(lines: readonly string[], i: number, indent: number): [MarkdownBlock, number] {
  const first = listItem(lines[i])!;
  const ordered = first[3] !== undefined;
  const items: MarkdownBlock[][] = [];
  let text: string[] = [];
  let blocks: MarkdownBlock[] = [];
  const endItem = () => {
    if (text.length) blocks.unshift({ kind: 'paragraph', children: parseMarkdownInline(text.join('\n')) });
    if (text.length || blocks.length) items.push(blocks);
    text = [];
    blocks = [];
  };
  while (i < lines.length) {
    const line = lines[i];
    if (!line.trim()) {
      // A blank line ends the list unless the next line goes on with it (an item or an indented line).
      const next = lines[i + 1];
      if (next === undefined || (!listItem(next) && indentOf(next) <= indent)) break;
      i++;
      continue;
    }
    const m = listItem(line);
    const at = indentOf(line);
    if (m && at === indent && (m[3] !== undefined) === ordered) {
      endItem();
      text = [m[4]];
      i++;
    } else if (m && at > indent) {
      const [nested, next] = parseList(lines, i, at);
      blocks.push(nested);
      i = next;
    } else if (!m && at > indent) {
      text.push(line.trim());
      i++;
    } else if (!m && text.length && !blocks.length) {
      // A lazy continuation line of the item's paragraph.
      text.push(line.trim());
      i++;
    } else break;
  }
  endItem();
  return [{ kind: 'list', ordered, start: ordered ? parseInt(first[3], 10) : 1, items }, i];
}

/** How deep strong, emphasis and links nest; deeper, their marks stay plain text (and the recursion stops). */
const MaxInlineDepth = 8;

/** The inline marks of a run of text (a paragraph, a heading, an item). */
export function parseMarkdownInline(text: string): MarkdownInline[] {
  return inline(text, 0);
}

const Space = /\s/;
/** `\s` by character code, the regular expression only past ASCII (the closer tests run on every character). */
const space = (code: number) =>
  code === 32 || (code >= 9 && code <= 13) || (code > 127 && Space.test(String.fromCharCode(code)));
/** A character at `c` that is not a space (`\S`); none outside the text. */
const nonSpace = (text: string, c: number) => c >= 0 && c < text.length && !space(text.charCodeAt(c));
const alnum = (text: string, c: number) => {
  const code = text.charCodeAt(c);
  return (code >= 0x30 && code <= 0x39) || ((code | 0x20) >= 0x61 && (code | 0x20) <= 0x7a);
};
const Star = 42;
const Under = 95;

/**
 * The first position at or after `from` that passes `test`, or -1; cached, so the queries of one left-to-right scan
 * (rising `from`) cost one pass over the text in all, however many openers ask.
 */
function seeker(text: string, test: (c: number) => boolean): (from: number) => number {
  let from = Infinity;
  let found = -1;
  return (p) => {
    if (p >= from && (found < 0 || p <= found)) return found;
    from = p;
    found = -1;
    for (let c = p; c < text.length; c++) {
      if (test(c)) {
        found = c;
        break;
      }
    }
    return found;
  };
}

/**
 * Where a link destination read from each position stops (a space, a `)` outside parentheses, the end), or -1 if a
 * `(` in it is not closed first: one backward pass, so balanced parentheses (`/wiki/Foo_(bar)`) cost no rescans.
 */
function destinationEnds(text: string): Int32Array {
  const n = text.length;
  const ends = new Int32Array(n + 1);
  ends[n] = n;
  for (let p = n - 1; p >= 0; p--) {
    const ch = text[p];
    if (ch === ')' || space(text.charCodeAt(p))) ends[p] = p;
    else if (ch === '\\' && p + 1 < n && Escapable.includes(text[p + 1])) ends[p] = ends[p + 2];
    else if (ch === '(') {
      const inner = ends[p + 1];
      ends[p] = inner >= 0 && text[inner] === ')' ? ends[inner + 1] : -1;
    } else ends[p] = ends[p + 1];
  }
  return ends;
}

function inline(text: string, depth: number): MarkdownInline[] {
  const n = text.length;
  const nest = depth < MaxInlineDepth;
  const out: MarkdownInline[] = [];
  let plain = '';
  let start = 0;
  const flush = () => {
    if (plain) out.push({ kind: 'text', text: decode(plain) });
    plain = '';
  };
  /** Ends the pending text at `i` (it goes to `plain`) and resumes it at `next`. */
  const take = (i: number, next: number) => {
    plain += text.slice(start, i);
    start = next;
  };

  const bracket = seeker(text, (c) => text[c] === ']');
  const quote = seeker(text, (c) => text[c] === '"');
  const gt = seeker(text, (c) => text[c] === '>');
  const spaceOrGt = seeker(text, (c) => text[c] === '>' || space(text.charCodeAt(c)));
  const notSpace = seeker(text, (c) => !space(text.charCodeAt(c)));
  const strong = (mark: number) => (c: number) => {
    if (text.charCodeAt(c) !== mark || text.charCodeAt(c + 1) !== mark || !nonSpace(text, c - 1)) return false;
    const after = text.charCodeAt(c + 2);
    return after !== Star && after !== Under;
  };
  const strongStars = seeker(text, strong(Star));
  const strongUnders = seeker(text, strong(Under));
  const emStar = seeker(text, (c) => text.charCodeAt(c) === Star && nonSpace(text, c - 1));
  const emUnder = seeker(text, (c) => text.charCodeAt(c) === Under && nonSpace(text, c - 1) && !alnum(text, c + 1));
  let destinations: Int32Array | undefined;
  let codeRuns: Map<number, { starts: number[]; next: number }> | undefined;

  /** A `[label](destination "title")` at `b` (the `[`): its parts and where it ends, or undefined. */
  const link = (b: number, minLabel: number) => {
    const close = bracket(b + 1);
    if (close < 0 || close - b - 1 < minLabel || text[close + 1] !== '(') return;
    destinations ??= destinationEnds(text);
    const destEnd = destinations[close + 2];
    if (destEnd < 0) return;
    let end = destEnd;
    if (text[end] !== ')') {
      const q = notSpace(end);
      if (q < 0 || text[q] !== '"') return;
      const r = quote(q + 1);
      if (r < 0 || text[r + 1] !== ')') return;
      end = r + 1;
    }
    return { label: text.slice(b + 1, close), href: text.slice(close + 2, destEnd), end: end + 1 };
  };

  /** The closing run of a code span opened by `k` backticks before `from`: a run of exactly `k`. */
  const codeCloser = (k: number, from: number) => {
    if (!codeRuns) {
      codeRuns = new Map();
      for (let c = text.indexOf('`'); c >= 0; ) {
        let e = c;
        while (text[e] === '`') e++;
        let runs = codeRuns.get(e - c);
        if (!runs) codeRuns.set(e - c, (runs = { starts: [], next: 0 }));
        runs.starts.push(c);
        c = text.indexOf('`', e);
      }
    }
    const runs = codeRuns.get(k);
    if (!runs) return -1;
    while (runs.next < runs.starts.length && runs.starts[runs.next] < from) runs.next++;
    return runs.next < runs.starts.length ? runs.starts[runs.next] : -1;
  };

  let i = 0;
  while (i < n) {
    const ch = text[i];
    if (ch === '\\') {
      if (Escapable.includes(text[i + 1] ?? '\n')) {
        take(i, i + 2);
        plain += text[i + 1];
        i += 2;
        continue;
      }
      if (text[i + 1] === '\n') {
        take(i, i + 2);
        flush();
        out.push({ kind: 'break' });
        i += 2;
        continue;
      }
    } else if (ch === '`') {
      let k = 1;
      while (text[i + k] === '`') k++;
      const c = codeCloser(k, i + k + 1);
      if (c >= 0) {
        take(i, c + k);
        flush();
        out.push({ kind: 'code', text: text.slice(i + k, c).replace(/^ (.*) $/, '$1') });
        i = c + k;
      } else i += k;
      continue;
    } else if (ch === '!' && text[i + 1] === '[') {
      const image = link(i + 1, 0);
      if (image) {
        take(i, image.end);
        plain += image.label;
        i = image.end;
        continue;
      }
    } else if (ch === '[' && nest) {
      const found = link(i, 1);
      if (found) {
        take(i, found.end);
        flush();
        out.push({ kind: 'link', href: decodeHref(found.href), children: inline(found.label, depth + 1) });
        i = found.end;
        continue;
      }
    } else if (ch === '<') {
      const scheme = /^(?:https?|mailto):/.exec(text.slice(i + 1, i + 9));
      const e = scheme ? spaceOrGt(i + 1 + scheme[0].length) : -1;
      if (scheme && e > i + 1 + scheme[0].length && text[e] === '>') {
        const href = text.slice(i + 1, e);
        take(i, e + 1);
        flush();
        out.push({ kind: 'link', href, children: [{ kind: 'text', text: href }] });
        i = e + 1;
        continue;
      }
      const name = text[i + 1] === '/' ? i + 2 : i + 1;
      const g = /[a-zA-Z]/.test(text[name] ?? '') ? gt(name + 1) : -1;
      if (g >= 0) {
        // An HTML tag, dropped.
        take(i, g + 1);
        i = g + 1;
        continue;
      }
    } else if ((ch === '*' || ch === '_') && nest) {
      if (text[i + 1] === ch && nonSpace(text, i + 2)) {
        const c = (ch === '*' ? strongStars : strongUnders)(i + 3);
        if (c >= 0) {
          take(i, c + 2);
          flush();
          out.push({ kind: 'strong', children: inline(text.slice(i + 2, c), depth + 1) });
          i = c + 2;
          continue;
        }
      }
      if (nonSpace(text, i + 1) && (ch === '*' || !alnum(text, i - 1))) {
        const c = (ch === '*' ? emStar : emUnder)(i + 2);
        if (c >= 0) {
          take(i, c + 1);
          flush();
          out.push({ kind: 'em', children: inline(text.slice(i + 1, c), depth + 1) });
          i = c + 1;
          continue;
        }
      }
    } else if (ch === ' ') {
      let e = i + 1;
      while (text[e] === ' ') e++;
      if (e - i >= 2 && text[e] === '\n') {
        take(i, e + 1);
        flush();
        out.push({ kind: 'break' });
        i = e + 1;
      } else i = e;
      continue;
    } else if (ch === '\n') {
      take(i, i + 1);
      plain += ' ';
      i++;
      continue;
    }
    i++;
  }
  take(n, n);
  flush();
  return out;
}

/** The text of inline marks (no emphasis): for labels that show plain text, and for search. */
export function markdownInlineText(nodes: readonly MarkdownInline[]): string {
  let s = '';
  for (const n of nodes) {
    if (n.kind === 'text' || n.kind === 'code') s += n.text;
    else if (n.kind === 'break') s += '\n';
    else s += markdownInlineText(n.children);
  }
  return s;
}
