/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type MarkdownBlock, type MarkdownInline, parseMarkdown } from '@molstar/core/util/markdown';
import type { HoverLine, HoverRun } from '@molstar/prism-ui/model/playback';

const Entities: Record<string, string> = { amp: '&', lt: '<', gt: '>', quot: '"', apos: "'", nbsp: ' ', '#39': "'" };

/**
 * Mol*'s loci labels (`plugin.behaviors.labels.highlight`: strings with `</br>` between lines and a `<small>` count)
 * as plain lines for the shell's hover label. Edge rate: once per hovered element, not per pointer move.
 */
export function hoverLines(labels: readonly unknown[]): readonly string[] {
  const out: string[] = [];
  for (const label of labels) {
    if (typeof label !== 'string') continue;
    const text = label
      .replace(/<\s*\/?\s*br\s*\/?\s*>/gi, '\n')
      .replace(/<[^>]*>/g, '')
      .replace(/&(#?\w+);/g, (m, name: string) => Entities[name] ?? m)
      .replace(/\s*\|\|\s*/g, ' ');
    for (const line of text.split('\n')) {
      const trimmed = line.trim();
      if (trimmed) out.push(trimmed);
    }
  }
  return out;
}

/**
 * The hover label's lines with their formatting, as v5 showed them: a label with a newline is Markdown (MolViewSpec's
 * tooltips, extensions), one without is Mol*'s HTML. Both become runs (`HoverRun`): text and its emphasis, never markup,
 * so a label from a loaded document cannot inject anything (the view builds elements from runs, `textContent` only).
 * Mol*'s HTML keeps `b`/`strong`, `i`/`em`, `small`, `sub`, `sup` and `br` (a line); other tags are dropped, their text
 * kept. Markdown keeps headings (a title line), list items (`-`, `*`, `+`, `1.`), paragraphs, `**bold**`, `*italic*`,
 * `` `code` `` and `[links](…)` (their text, shown as a link: a label takes no pointer); HTML in it is skipped, as v5's
 * `skipHtml`. Plain lines (no emphasis but a Mol* title) stay strings.
 */
export function hoverRich(labels: readonly unknown[]): readonly HoverLine[] {
  const out: HoverLine[] = [];
  for (const label of labels) {
    if (typeof label !== 'string' || !label.trim()) continue;
    if (label.includes('\n')) markdownLines(label, out);
    else htmlLines(label, out, out.length === 0);
  }
  return out;
}

const decode = (s: string) => s.replace(/&(#?\w+);/g, (m, name: string) => Entities[name] ?? m);

/** Mol*'s HTML: a line per `<br>`, runs by the allowed tags. The first line of the label is its title. */
function htmlLines(label: string, out: HoverLine[], first: boolean): void {
  const text = label.replace(/\s*\|\|\s*/g, ' ');
  let runs: HoverRun[] = [];
  const style = { bold: 0, italic: 0, small: 0, sub: 0, sup: 0 };
  const flush = () => {
    const merged = trimRuns(runs);
    runs = [];
    if (!merged.length) return;
    const role = first && out.length === 0 ? 'title' : 'text';
    out.push(merged.every(plain) ? merged.map((r) => r.text).join('') : { runs: merged, role });
  };
  const re = /<\s*(\/?)\s*([a-zA-Z]+)[^>]*>/g;
  let last = 0;
  for (let m = re.exec(text); m; m = re.exec(text)) {
    push(runs, decode(text.slice(last, m.index)), style);
    last = re.lastIndex;
    const close = m[1] === '/';
    const tag = m[2].toLowerCase();
    if (tag === 'br') flush();
    else if (tag === 'b' || tag === 'strong') style.bold += close ? -1 : 1;
    else if (tag === 'i' || tag === 'em') style.italic += close ? -1 : 1;
    else if (tag === 'small' || tag === 'sub' || tag === 'sup') style[tag] += close ? -1 : 1;
  }
  push(runs, decode(text.slice(last)), style);
  flush();
}

/**
 * Markdown (core's tokenizer, `@molstar/core/util/markdown`) as lines: a heading a title, a paragraph a line per hard
 * break, list items with their marker (nested ones indented), code blocks a line each; rules dropped.
 */
function markdownLines(label: string, out: HoverLine[]): void {
  const blocks = (list: readonly MarkdownBlock[], depth: number) => {
    for (const b of list) {
      if (b.kind === 'heading') out.push({ runs: runsOf(b.children), role: 'title' });
      else if (b.kind === 'paragraph')
        for (const runs of splitBreaks(runsOf(b.children))) out.push({ runs, role: 'text' });
      else if (b.kind === 'code')
        for (const line of b.text.split('\n')) out.push({ runs: [{ text: line, code: true }], role: 'text' });
      else if (b.kind === 'list') {
        b.items.forEach((item, k) => {
          const marker = `${'  '.repeat(depth)}${b.ordered ? `${b.start + k}.` : depth ? '◦' : '•'} `;
          const [head, ...rest] = item;
          if (head?.kind === 'paragraph') {
            const [first, ...more] = splitBreaks(runsOf(head.children));
            out.push({ runs: [{ text: marker }, ...(first ?? [])], role: 'item' });
            for (const runs of more) out.push({ runs, role: 'item' });
            blocks(rest, depth + 1);
          } else {
            out.push({ runs: [{ text: marker }], role: 'item' });
            blocks(item, depth + 1);
          }
        });
      }
    }
  };
  blocks(parseMarkdown(label), 0);
}

/** Inline marks as runs (emphasis flattened into flags); a break stays a `\n` run for `splitBreaks`. */
function runsOf(nodes: readonly MarkdownInline[], flags: Omit<HoverRun, 'text'> = {}): HoverRun[] {
  const out: HoverRun[] = [];
  for (const n of nodes) {
    if (n.kind === 'text') out.push({ ...flags, text: n.text });
    else if (n.kind === 'code') out.push({ ...flags, text: n.text, code: true });
    else if (n.kind === 'break') out.push({ text: '\n' });
    else if (n.kind === 'strong') out.push(...runsOf(n.children, { ...flags, bold: true }));
    else if (n.kind === 'em') out.push(...runsOf(n.children, { ...flags, italic: true }));
    else out.push(...runsOf(n.children, { ...flags, link: true }));
  }
  return out;
}

function splitBreaks(runs: readonly HoverRun[]): HoverRun[][] {
  const lines: HoverRun[][] = [[]];
  for (const r of runs) {
    if (r.text === '\n') lines.push([]);
    else lines[lines.length - 1].push(r);
  }
  return lines.map(trimRuns).filter((l) => l.length > 0);
}

function push(
  runs: HoverRun[],
  text: string,
  style: { bold: number; italic: number; small: number; sub: number; sup: number },
) {
  if (!text) return;
  const run: HoverRun = {
    text,
    ...(style.bold > 0 ? { bold: true } : {}),
    ...(style.italic > 0 ? { italic: true } : {}),
    ...(style.small > 0 ? { small: true } : {}),
    ...(style.sub > 0 ? { sub: true } : {}),
    ...(style.sup > 0 ? { sup: true } : {}),
  };
  runs.push(run);
}

const plain = (r: HoverRun) => !r.bold && !r.italic && !r.small && !r.sub && !r.sup && !r.code && !r.link;

/** Whitespace collapsed and trimmed at the line's ends; runs of the same style merged. */
function trimRuns(runs: readonly HoverRun[]): HoverRun[] {
  const out: HoverRun[] = [];
  for (const r of runs) {
    let text = r.text.replace(/\s+/g, ' ');
    const prev = out[out.length - 1];
    // One space across a run boundary (HTML collapses it too).
    if (prev?.text.endsWith(' ') && text.startsWith(' ')) text = text.slice(1);
    if (prev && sameStyle(prev, r)) out[out.length - 1] = { ...prev, text: prev.text + text };
    else out.push({ ...r, text });
  }
  if (out.length) {
    out[0] = { ...out[0], text: out[0].text.trimStart() };
    out[out.length - 1] = { ...out[out.length - 1], text: out[out.length - 1].text.trimEnd() };
  }
  return out.filter((r) => r.text);
}

const sameStyle = (a: HoverRun, b: HoverRun) =>
  !!a.bold === !!b.bold &&
  !!a.italic === !!b.italic &&
  !!a.small === !!b.small &&
  !!a.sub === !!b.sub &&
  !!a.sup === !!b.sup &&
  !!a.code === !!b.code &&
  !!a.link === !!b.link;
