/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { decodeHtmlEntities } from '../html-entities.js';

describe('decodeHtmlEntities', () => {
  it('reads decimal, hexadecimal and named references', () => {
    expect(decodeHtmlEntities('&#x2B05; Back')).toBe('⬅ Back');
    expect(decodeHtmlEntities('a&nbsp;b &amp; c &lt;b&gt; &#65;&#X42; &mdash; &Aring; &alpha;')).toBe(
      'a b & c <b> AB — Å α',
    );
    expect(decodeHtmlEntities('&#x1F600;')).toBe('\u{1F600}');
  });

  it('keeps unknown names as written; invalid numbers are U+FFFD; no prototype lookups', () => {
    expect(decodeHtmlEntities('&unknownname; &constructor; &toString;')).toBe('&unknownname; &constructor; &toString;');
    expect(decodeHtmlEntities('&#0; &#xD800; &#x110000;')).toBe('� � �');
    expect(decodeHtmlEntities('AT&T & co; &amp')).toBe('AT&T & co; &amp');
    expect(decodeHtmlEntities('plain')).toBe('plain');
  });
});
