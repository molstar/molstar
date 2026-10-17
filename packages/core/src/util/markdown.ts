/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A small Markdown tokenizer for labels, tooltips and short descriptions: a typed tree, never HTML, so whatever renders
// it builds elements from text (a document's Markdown cannot inject markup). The approach follows snarkdown (Jason
// Miller, MIT): a line pass for the blocks and one scanner for the inline marks. What it reads: ATX headings, paragraphs
// (a single newline is a space; two trailing spaces or a backslash make a break), bullet and ordered lists (nested by
// indentation), fenced code, rules; inline: code spans, strong and emphasis (stars or underscores, underscores not inside
// words), links and autolinks, image alt text, backslash escapes, the common entities. HTML is dropped (as
// react-markdown's `skipHtml`). Not read: tables, block quotes, setext headings, reference links.

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

const Entities: Record<string, string> = { amp: '&', lt: '<', gt: '>', quot: '"', apos: "'", nbsp: ' ', '#39': "'" };

function decode(s: string): string {
  return s.replace(/&(#x[0-9a-f]+|#\d+|\w+);/gi, (m, name: string) => {
    if (name[0] === '#') {
      const code = name[1] === 'x' || name[1] === 'X' ? parseInt(name.slice(2), 16) : parseInt(name.slice(1), 10);
      return Number.isFinite(code) && code > 0 && code < 0x110000 ? String.fromCodePoint(code) : m;
    }
    return Entities[name] ?? m;
  });
}

/** Whether a link's target is one a renderer may follow: http(s), mailto, or relative (no scheme). */
export function isSafeMarkdownHref(href: string): boolean {
  const scheme = /^([a-z][a-z0-9+.-]*):/i.exec(href.trim());
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
    const heading = /^(#{1,6})\s+(.*?)\s*#*\s*$/.exec(trimmed);
    if (heading) {
      endParagraph();
      blocks.push({ kind: 'heading', level: heading[1].length, children: parseMarkdownInline(heading[2]) });
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

const ListItem = /^(\s*)([-*+]|(\d+)[.)])\s+(.*)$/;

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

/** The inline marks of a run of text (a paragraph, a heading, an item). */
export function parseMarkdownInline(text: string): MarkdownInline[] {
  const out: MarkdownInline[] = [];
  let plain = '';
  const flush = () => {
    if (plain) out.push({ kind: 'text', text: decode(plain) });
    plain = '';
  };
  const re =
    /\\([\\`*_{}[\]()#+\-.!<>|~])|(`+)([\s\S]+?)\2(?!`)|!\[([^\]]*)\]\(([^)\s]*)(?:\s+"[^"]*")?\)|\[([^\]]+)\]\(([^)\s]*)(?:\s+"[^"]*")?\)|<((?:https?|mailto):[^\s>]+)>|<\/?[a-zA-Z][^>]*>|(\*\*|__)(\S(?:[\s\S]*?\S)?)\9(?![*_])|\*(\S(?:[\s\S]*?\S)?)\*|(?<![A-Za-z0-9])_(\S(?:[\s\S]*?\S)?)_(?![A-Za-z0-9])|( {2,}|\\)\n|\n/g;
  let last = 0;
  for (let m = re.exec(text); m; m = re.exec(text)) {
    plain += text.slice(last, m.index);
    last = re.lastIndex;
    if (m[1] !== undefined) plain += m[1];
    else if (m[2] !== undefined) {
      flush();
      out.push({ kind: 'code', text: m[3].replace(/^ (.*) $/, '$1') });
    } else if (m[4] !== undefined) plain += m[4];
    else if (m[6] !== undefined) {
      flush();
      out.push({ kind: 'link', href: m[7], children: parseMarkdownInline(m[6]) });
    } else if (m[8] !== undefined) {
      flush();
      out.push({ kind: 'link', href: m[8], children: [{ kind: 'text', text: m[8] }] });
    } else if (m[10] !== undefined) {
      flush();
      out.push({ kind: 'strong', children: parseMarkdownInline(m[10]) });
    } else if (m[11] !== undefined || m[12] !== undefined) {
      flush();
      out.push({ kind: 'em', children: parseMarkdownInline(m[11] ?? m[12]) });
    } else if (m[13] !== undefined) {
      flush();
      out.push({ kind: 'break' });
    } else if (m[0] === '\n') plain += ' ';
    // else: an HTML tag, dropped.
  }
  plain += text.slice(last);
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
