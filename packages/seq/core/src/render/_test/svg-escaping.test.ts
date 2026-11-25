// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The SVG backend's escaping (_spec/README.md, Render, Backends; the security bounds): labels, ids, titles, fonts and
// metadata from a spec or data, hostile as they can be, never become markup: the document parses as XML, holds only
// the elements the backend writes, and no attribute or text is anything but what was given (a character XML forbids
// replaced one for one).

import { describe, expect, it } from 'vitest';
import { Anchor, DrawLists, drawListsToSvg, escapeXml, PassWriter } from '../index.js';

const Hostile: readonly string[] = [
  '"><script>alert(1)</script>',
  "' onload='alert(1)",
  '" onload="alert(1)',
  ']]>',
  '<![CDATA[x]]>',
  '<!-- c -->',
  '-->',
  '<?pi x?>',
  '&amp; &lt; &#x41; &#65; &unknown;',
  '</text><text>',
  '\u0000\u0001\u0008\u000b\u000c\u000e\u001f￾￿',
  '\ud800',
  '\udc00',
  'a\ud800b',
  '\udc00\ud800',
  '😀',
  '  \u0085',
  '<svg onload="x">',
];

const Elements = new Set([
  'svg',
  'title',
  'defs',
  'g',
  'rect',
  'path',
  'text',
  'tspan',
  'clipPath',
  'mask',
  'linearGradient',
  'stop',
]);

/** What the text should read back as: forbidden characters and lone surrogates replaced one for one. */
const cleaned = (s: string) =>
  s
    .replace(/[\u0000-\u0008\u000b\u000c\u000e-\u001f￾￿]/g, '�')
    .replace(/[\ud800-\udbff](?![\udc00-\udfff])|(?<![\ud800-\udbff])[\udc00-\udfff]/g, '�');

function parse(svg: string): Document {
  const doc = new DOMParser().parseFromString(svg, 'image/svg+xml');
  expect(doc.getElementsByTagName('parsererror')).toHaveLength(0);
  return doc;
}

function onlyOurMarkup(doc: Document) {
  for (const el of Array.from(doc.getElementsByTagName('*'))) expect(Elements.has(el.localName)).toBe(true);
  for (const el of Array.from(doc.getElementsByTagName('*')))
    for (const a of Array.from(el.attributes)) expect(a.name).not.toMatch(/^on/i);
}

describe('escapeXml', () => {
  it('writes no raw markup character, and replaces what XML forbids one for one', () => {
    for (const s of Hostile) {
      const out = escapeXml(s);
      expect(out).not.toMatch(/[<>"']/);
      expect(out).not.toMatch(/&(?!amp;|lt;|gt;|quot;|#39;)/);
      expect(out.length).toBeGreaterThanOrEqual(s.length);
      expect(out).toBe(
        cleaned(s).replace(
          /[&<>"']/g,
          (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' })[c]!,
        ),
      );
    }
  });

  it('keeps a surrogate pair and replaces a lone half', () => {
    expect(escapeXml('😀')).toBe('😀');
    expect(escapeXml('\ud800')).toBe('�');
    expect(escapeXml('\udc00')).toBe('�');
    expect(escapeXml('a\ud800b')).toBe('a�b');
    expect(escapeXml('\udc00\ud800')).toBe('��');
    expect(escapeXml(']]>')).toBe(']]&gt;');
  });
});

describe('a document with hostile strings', () => {
  function lists(label: string, font = '12px sans'): DrawLists {
    const l = new DrawLists();
    l.width = 120;
    l.height = 40;
    const labels = new PassWriter(l, l.passes.labels);
    // wide text, so a fade (a mask with an id) is written too
    labels.text(label, 2, 12, Anchor.start, font, 0x000000ff, 400, 40);
    labels.beginGlyphs(font, 20, 0x000000ff);
    for (let i = 0; i < label.length; i++) labels.letter(label.charCodeAt(i), 10 + i);
    labels.endGlyphs();
    return l;
  }

  it('parses as XML and holds only the backend’s elements, whatever the label', () => {
    for (const s of Hostile) {
      const doc = parse(drawListsToSvg(lists(s)));
      onlyOurMarkup(doc);
      const texts = Array.from(doc.getElementsByTagName('text'));
      expect(texts[0].textContent).toBe(cleaned(s));
      // the glyph run: one character a glyph, each replaced as the label's
      expect(texts[1].textContent).toBe(cleaned(s));
    }
  });

  it('escapes a title, metadata values and a font that is none of the roles', () => {
    for (const s of Hostile) {
      const svg = drawListsToSvg(lists('x', `12px ${s}`), { title: s, metadata: { 'seq-x': s } });
      const doc = parse(svg);
      onlyOurMarkup(doc);
      expect(doc.getElementsByTagName('title')[0].textContent).toBe(cleaned(s));
      expect(doc.documentElement.getAttribute('data-seq-x')).toBe(cleaned(s));
      expect(doc.getElementsByTagName('text')[0].getAttribute('font-family')).toBe(cleaned(`12px ${s}`));
    }
  });

  it('writes a metadata key only when it is a plain name', () => {
    const svg = drawListsToSvg(lists('x'), {
      metadata: { 'seq-ok': '1', 'a"b': '2', 'x onload': '3', 'a>b': '4', 'ok\n': '5', '': '6' },
    });
    const doc = parse(svg);
    expect(Array.from(doc.documentElement.attributes).map((a) => a.name)).toEqual([
      'xmlns',
      'width',
      'height',
      'viewBox',
      'data-seq-ok',
    ]);
  });

  it('keeps an id prefix to name characters: no quote or bracket gets into an id or a url()', () => {
    const svg = drawListsToSvg(lists('a long label'), { idPrefix: '"/><x onload=y>)(#' });
    const doc = parse(svg);
    onlyOurMarkup(doc);
    for (const el of Array.from(doc.getElementsByTagName('*'))) {
      const id = el.getAttribute('id');
      if (id !== null) expect(id).toMatch(/^[A-Za-z0-9_-]+$/);
      const mask = el.getAttribute('mask');
      if (mask !== null) expect(mask).toMatch(/^url\(#[A-Za-z0-9_-]+\)$/);
    }
  });
});
