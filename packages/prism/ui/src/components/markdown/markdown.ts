/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { bind, getScope, onCleanup } from '../../dom/scope.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { classes } from './_generated/markdown.module.css.js';

/** What a host's link does, per event: a click, and the pointer entering and leaving it (a highlight). */
export type MarkdownLinkEvent = 'click' | 'enter' | 'leave';

/**
 * What the page that shows the markdown adds to it (the viewer: Mol*'s markdown extensions), so the renderer itself
 * knows no plugin: links it runs (`!command` links), images it resolves (an archive's files, inline commands such as a
 * color swatch), and media to show above the text (the audio player a command made).
 */
export interface MarkdownHost {
  /** The host's handler for a link, when it is one of its own (a command); otherwise the renderer's own handling. */
  link?(href: string): ((event: MarkdownLinkEvent) => void) | undefined;
  /**
   * An image: its URL (or a promise of one, for a file read from an archive), an element to show in its place (an
   * inline command), or undefined to load `src` as it is.
   */
  image?(src: string, alt: string): string | Promise<string | undefined> | Node | undefined;
  /** Shown above the text while set (one node, moved here from the last description that showed it). */
  media?: ReadonlyState<Node | null | undefined>;
}

/**
 * A small, safe markdown renderer for descriptions (no HTML is parsed; everything becomes DOM nodes): headings,
 * paragraphs, lists, tables, fenced code, rules, **bold**, *italic*, `code`, images and [links](href). A link the host
 * claims goes to it; `#key` and other non-web links go to `onLink`; web links open in a new tab.
 */
export function Markdown(props: { source: string; onLink?: (href: string) => void; host?: MarkdownHost }): HTMLElement {
  const root = document.createElement('div');
  root.className = classes.markdown;
  const ctx: InlineContext = { onLink: props.onLink, host: props.host };
  const media = props.host?.media;
  if (media) {
    const slot = document.createElement('div');
    slot.className = classes.media;
    root.appendChild(slot);
    const show = (node: Node | null | undefined) => {
      slot.replaceChildren(...(node ? [node] : []));
      slot.hidden = !node;
    };
    // Outside a view (a test, a one-off render) the current node is shown once.
    if (getScope()) bind(media, show);
    else show(media.value);
  }
  const lines = props.source.split('\n');
  let list: HTMLUListElement | HTMLOListElement | undefined;
  let paragraph: HTMLParagraphElement | undefined;
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i].trimEnd();
    const heading = /^(#{1,4})\s+(.*)$/.exec(line);
    const bullet = /^\s*[-*]\s+(.*)$/.exec(line);
    const numbered = /^\s*\d+\.\s+(.*)$/.exec(line);
    if (!line.trim()) {
      list = paragraph = undefined;
    } else if (/^\s*```/.test(line)) {
      list = paragraph = undefined;
      const code: string[] = [];
      for (i++; i < lines.length && !/^\s*```/.test(lines[i]); i++) code.push(lines[i]);
      const pre = document.createElement('pre');
      pre.appendChild(element('code', code.join('\n')));
      root.appendChild(pre);
    } else if (isTableStart(lines, i)) {
      list = paragraph = undefined;
      i = table(root, lines, i, ctx) - 1;
    } else if (/^\s*([-*_])(\s*\1){2,}\s*$/.test(line)) {
      list = paragraph = undefined;
      root.appendChild(document.createElement('hr'));
    } else if (heading) {
      list = paragraph = undefined;
      const h = document.createElement(`h${Math.min(6, heading[1].length + 2)}`);
      inline(h, heading[2], ctx);
      root.appendChild(h);
    } else if (bullet || numbered) {
      paragraph = undefined;
      const tag = bullet ? 'UL' : 'OL';
      if (!list || list.tagName !== tag) {
        list = document.createElement(bullet ? 'ul' : 'ol');
        root.appendChild(list);
      }
      const li = document.createElement('li');
      inline(li, (bullet ?? numbered)![1], ctx);
      list.appendChild(li);
    } else {
      list = undefined;
      if (!paragraph) {
        paragraph = document.createElement('p');
        root.appendChild(paragraph);
      } else {
        paragraph.appendChild(document.createTextNode(' '));
      }
      inline(paragraph, line.trim(), ctx);
    }
  }
  return root;
}

interface InlineContext {
  onLink?: (href: string) => void;
  host?: MarkdownHost;
}

// GFM tables: a header row, a separator row (`|---|:-:|`), then rows until a line without a pipe.
const TableSeparator = /^\s*\|?\s*:?-+:?\s*(\|\s*:?-+:?\s*)*\|?\s*$/;

function isTableStart(lines: readonly string[], i: number): boolean {
  return (
    lines[i].includes('|') && i + 1 < lines.length && TableSeparator.test(lines[i + 1]) && lines[i + 1].includes('-')
  );
}

/** The cells of a table row: split on pipes that are not escaped, the outer ones dropped. */
function cells(line: string): string[] {
  let text = line.trim();
  if (text.startsWith('|')) text = text.slice(1);
  if (text.endsWith('|') && !text.endsWith('\\|')) text = text.slice(0, -1);
  return text.split(/(?<!\\)\|/).map((c) => c.trim().replace(/\\\|/g, '|'));
}

/** Renders the table at line `start`; returns the line after it. */
function table(root: HTMLElement, lines: readonly string[], start: number, ctx: InlineContext): number {
  const align = cells(lines[start + 1]).map((c) =>
    c.startsWith(':') && c.endsWith(':') ? 'center' : c.endsWith(':') ? 'right' : '',
  );
  const el = document.createElement('table');
  const row = (parent: HTMLElement, text: string, tag: 'th' | 'td') => {
    const tr = document.createElement('tr');
    const values = cells(text);
    for (let c = 0; c < align.length; c++) {
      const cell = document.createElement(tag);
      if (align[c]) cell.style.textAlign = align[c];
      inline(cell, values[c] ?? '', ctx);
      tr.appendChild(cell);
    }
    parent.appendChild(tr);
  };
  const head = document.createElement('thead');
  row(head, lines[start], 'th');
  el.appendChild(head);
  const body = document.createElement('tbody');
  let i = start + 2;
  for (; i < lines.length && lines[i].trim() && lines[i].includes('|'); i++) row(body, lines[i], 'td');
  if (body.childNodes.length) el.appendChild(body);
  // A wide table scrolls by itself, not the description.
  const wrap = document.createElement('div');
  wrap.className = classes.table;
  wrap.appendChild(el);
  root.appendChild(wrap);
  return i;
}

const InlinePattern = /\*\*([^*]+)\*\*|\*([^*]+)\*|`([^`]+)`|!\[([^\]]*)\]\(([^)\s]+)\)|\[([^\]]+)\]\(([^)\s]+)\)/g;

function inline(parent: Node, text: string, ctx: InlineContext) {
  let last = 0;
  // A fresh pattern per call: a link's text is rendered by a nested call.
  const pattern = new RegExp(InlinePattern.source, 'g');
  for (let m = pattern.exec(text); m; m = pattern.exec(text)) {
    if (m.index > last) parent.appendChild(document.createTextNode(text.slice(last, m.index)));
    if (m[1] !== undefined) parent.appendChild(element('strong', m[1]));
    else if (m[2] !== undefined) parent.appendChild(element('em', m[2]));
    else if (m[3] !== undefined) parent.appendChild(element('code', m[3]));
    else if (m[5] !== undefined) parent.appendChild(image(m[5], m[4], ctx));
    else parent.appendChild(link(m[7], m[6], ctx));
    last = m.index + m[0].length;
  }
  if (last < text.length) parent.appendChild(document.createTextNode(text.slice(last)));
}

function link(href: string, text: string, ctx: InlineContext): HTMLAnchorElement {
  const a = document.createElement('a');
  inline(a, text, ctx);
  const own = ctx.host?.link?.(href);
  if (own) {
    a.href = '#';
    a.addEventListener('click', (e) => {
      e.preventDefault();
      own('click');
    });
    let over = false;
    a.addEventListener('pointerenter', () => {
      over = true;
      own('enter');
    });
    a.addEventListener('pointerleave', () => {
      over = false;
      own('leave');
    });
    // A link removed under the pointer (the next snapshot's description) gets no pointerleave: its highlight would stay.
    if (getScope()) onCleanup(() => over && own('leave'));
    return a;
  }
  const web = /^https?:/i.test(href);
  a.href = href.startsWith('#') || web ? href : '#';
  if (web) {
    a.target = '_blank';
    a.rel = 'noopener noreferrer';
  } else {
    a.addEventListener('click', (e) => {
      e.preventDefault();
      ctx.onLink?.(href);
    });
  }
  return a;
}

function image(src: string, alt: string, ctx: InlineContext): Node {
  const resolved = ctx.host?.image?.(src, alt);
  if (resolved instanceof Node) return resolved;
  const img = document.createElement('img');
  img.alt = alt;
  if (typeof resolved === 'string') img.src = resolved;
  else if (resolved) {
    // An archive's file: shown once it is read; nothing when it cannot be.
    img.hidden = true;
    resolved.then(
      (url) => {
        if (!url) return;
        img.src = url;
        img.hidden = false;
      },
      () => {},
    );
  } else img.src = src;
  return img;
}

function element(tag: string, text: string) {
  const el = document.createElement(tag);
  el.textContent = text;
  return el;
}
