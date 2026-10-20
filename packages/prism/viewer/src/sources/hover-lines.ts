/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { decodeHtmlEntities } from '@molstar/core/util/html-entities';
import { type MarkdownBlock, type MarkdownInline, parseMarkdown } from '@molstar/core/util/markdown';
import type { HoverLine, HoverRun } from '@molstar/prism-ui/model/playback';

const isSpace = (c: number) => c === 32 || (c >= 9 && c <= 13) || (c > 127 && /\s/.test(String.fromCharCode(c)));

/** A `<…>` of a label: where it is and, when it reads as a tag (`<b>`, `</ b >`, `<br/>`), its name (lower case). */
interface TagHit {
  readonly start: number;
  readonly end: number;
  readonly close: boolean;
  readonly name: string;
}

/**
 * The tags of `text` in one pass (labels come from loaded documents: no regex with adjacent quantifiers over runs of
 * `<` or spaces, which are quadratic). A tag starts at the last `<` before its `>`; a `<` with no `>` after it is text.
 */
function* scanTags(text: string): Generator<TagHit> {
  let lt = -1;
  for (let p = 0; p < text.length; p++) {
    const c = text.charCodeAt(p);
    if (c === 60) lt = p;
    else if (c === 62 && lt >= 0) {
      let i = lt + 1;
      while (i < p && isSpace(text.charCodeAt(i))) i++;
      const close = text.charCodeAt(i) === 47;
      if (close) i++;
      while (i < p && isSpace(text.charCodeAt(i))) i++;
      const from = i;
      for (; i < p; i++) {
        const l = text.charCodeAt(i) | 32;
        if (l < 97 || l > 122) break;
      }
      yield { start: lt, end: p + 1, close, name: text.slice(from, i).toLowerCase() };
      lt = -1;
    }
  }
}

/** `\s*\|\|\s*` → one space, in one pass. */
function collapseBars(text: string): string {
  if (!text.includes('||')) return text;
  let out = '';
  let i = 0;
  for (let bar = text.indexOf('||'); bar >= 0; bar = text.indexOf('||', i)) {
    let end = bar;
    while (end > i && isSpace(text.charCodeAt(end - 1))) end--;
    out += `${text.slice(i, end)} `;
    i = bar + 2;
    while (i < text.length && isSpace(text.charCodeAt(i))) i++;
  }
  return out + text.slice(i);
}

/**
 * Mol*'s loci labels (`plugin.behaviors.labels.highlight`: strings with `</br>` between lines and a `<small>` count)
 * as plain lines for the shell's hover label. Edge rate: once per hovered element, not per pointer move.
 */
export function hoverLines(labels: readonly unknown[]): readonly string[] {
  const out: string[] = [];
  for (const label of labels) {
    if (typeof label !== 'string') continue;
    let stripped = '';
    let last = 0;
    for (const tag of scanTags(label)) {
      stripped += label.slice(last, tag.start);
      if (tag.name === 'br') stripped += '\n';
      last = tag.end;
    }
    const text = collapseBars(decode(stripped + label.slice(last)));
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

/** A label's references as characters (core's decoder); a no-break space a plain one, as before. */
const decode = (s: string) => decodeHtmlEntities(s).replace(/\u00a0/g, ' ');

/** Mol*'s HTML: a line per `<br>`, runs by the allowed tags. The first line of the label is its title. */
function htmlLines(label: string, out: HoverLine[], first: boolean): void {
  const text = collapseBars(label);
  let runs: HoverRun[] = [];
  const style = { bold: 0, italic: 0, small: 0, sub: 0, sup: 0 };
  const flush = () => {
    const merged = trimRuns(runs);
    runs = [];
    if (!merged.length) return;
    const role = first && out.length === 0 ? 'title' : 'text';
    out.push(merged.every(plain) ? merged.map((r) => r.text).join('') : { runs: merged, role });
  };
  let last = 0;
  for (const hit of scanTags(text)) {
    // Not a tag (`<3 >`): text.
    if (!hit.name) continue;
    push(runs, decode(text.slice(last, hit.start)), style);
    last = hit.end;
    const close = hit.close;
    const tag = hit.name;
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

/**
 * Which label the canvas shows (sequence.md §10): Mol*'s whenever there is one (a sequence hover the link maps into the
 * structure is lit in 3D, and Mol* labels it), the sequence's own for a hover that does not map, else none (a mapped
 * hover whose Mol* label has not come yet: the label holds a frame meanwhile).
 */
export function pickHoverLabel(
  molstar: readonly HoverLine[],
  sequence: readonly HoverLine[],
  sequenceMapped: boolean,
): readonly HoverLine[] {
  if (molstar.length) return molstar;
  return sequence.length && !sequenceMapped ? sequence : [];
}
