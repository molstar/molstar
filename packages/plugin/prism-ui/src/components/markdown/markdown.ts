/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { classes } from './_generated/markdown.module.css.js';

/**
 * A small, safe markdown renderer for descriptions (no HTML is parsed; everything becomes DOM nodes): headings,
 * paragraphs, lists, **bold**, *italic*, `code` and [links](href). Links go to `onLink`, which applies `#key` snapshot
 * links and commands. A stand-in until the plugin's markdown extensions are ported.
 */
export function Markdown(props: { source: string; onLink?: (href: string) => void }): HTMLElement {
  const root = document.createElement('div');
  root.className = classes.markdown;
  let list: HTMLUListElement | HTMLOListElement | undefined;
  let paragraph: HTMLParagraphElement | undefined;
  for (const raw of props.source.split('\n')) {
    const line = raw.trimEnd();
    const heading = /^(#{1,4})\s+(.*)$/.exec(line);
    const bullet = /^\s*[-*]\s+(.*)$/.exec(line);
    const numbered = /^\s*\d+\.\s+(.*)$/.exec(line);
    if (!line.trim()) {
      list = undefined;
      paragraph = undefined;
    } else if (heading) {
      list = paragraph = undefined;
      const h = document.createElement(`h${Math.min(6, heading[1].length + 2)}`);
      inline(h, heading[2], props.onLink);
      root.appendChild(h);
    } else if (bullet || numbered) {
      paragraph = undefined;
      const tag = bullet ? 'UL' : 'OL';
      if (!list || list.tagName !== tag) {
        list = document.createElement(bullet ? 'ul' : 'ol');
        root.appendChild(list);
      }
      const li = document.createElement('li');
      inline(li, (bullet ?? numbered)![1], props.onLink);
      list.appendChild(li);
    } else {
      list = undefined;
      if (!paragraph) {
        paragraph = document.createElement('p');
        root.appendChild(paragraph);
      } else {
        paragraph.appendChild(document.createTextNode(' '));
      }
      inline(paragraph, line.trim(), props.onLink);
    }
  }
  return root;
}

const InlinePattern = /\*\*([^*]+)\*\*|\*([^*]+)\*|`([^`]+)`|\[([^\]]+)\]\(([^)\s]+)\)/g;

function inline(parent: Node, text: string, onLink?: (href: string) => void) {
  let last = 0;
  InlinePattern.lastIndex = 0;
  for (let m = InlinePattern.exec(text); m; m = InlinePattern.exec(text)) {
    if (m.index > last) parent.appendChild(document.createTextNode(text.slice(last, m.index)));
    if (m[1] !== undefined) parent.appendChild(element('strong', m[1]));
    else if (m[2] !== undefined) parent.appendChild(element('em', m[2]));
    else if (m[3] !== undefined) parent.appendChild(element('code', m[3]));
    else {
      const a = element('a', m[4]) as HTMLAnchorElement;
      const href = m[5];
      a.href = href.startsWith('#') || /^https?:/.test(href) ? href : '#';
      if (/^https?:/.test(href)) {
        a.target = '_blank';
        a.rel = 'noopener noreferrer';
      } else {
        a.addEventListener('click', (e) => {
          e.preventDefault();
          onLink?.(href);
        });
      }
      parent.appendChild(a);
    }
    last = m.index + m[0].length;
  }
  if (last < text.length) parent.appendChild(document.createTextNode(text.slice(last)));
}

function element(tag: string, text: string) {
  const el = document.createElement(tag);
  el.textContent = text;
  return el;
}
