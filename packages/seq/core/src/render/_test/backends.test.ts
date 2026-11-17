/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The backends: the recording format line by line, the SVG backend's groups, clips, squircles, images and fades, and
// escaping (labels from data never become markup); the token table, the text port and the squircle's curve.

import { categoryColumn } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import {
  Anchor,
  bundledTextMetrics,
  bundledWidth,
  DrawLists,
  drawListsToSvg,
  escapeXml,
  inkOn,
  LightTokens,
  luminance,
  PassWriter,
  parseColor,
  Round,
  recordDrawLists,
  renderSvg,
  resolveColor,
  SquircleHandle,
  squircleRadius,
} from '../index.js';

function sample(): DrawLists {
  const lists = new DrawLists();
  lists.width = 200;
  lists.height = 40;
  const data = new PassWriter(lists, lists.passes.data);
  data.clip(10, 0, 180, 40);
  data.rect(10, 4, 5.5, 16, 0xff0000ff);
  data.rect(20, 4, 60, 12, 0x7a8fa6ff, 0, 0, 4, Round.left);
  data.beginGlyphs('10px mono', 15.333, 0x2d271eff);
  data.letter(77, 12.75);
  data.letter(69, 18.25);
  data.endGlyphs();
  data.image(10, 30, 180, 4, {
    width: 3,
    height: 1,
    data: Uint8Array.from([255, 0, 0, 255, 255, 0, 0, 255, 0, 0, 0, 0]),
  });
  data.unclip();
  const labels = new PassWriter(lists, lists.passes.labels);
  labels.text('A <long> label', 2, 12, Anchor.start, '12px sans', 0x6b655cff, 80, 40);
  return lists;
}

describe('the recording backend', () => {
  it('writes a line a primitive, pass by pass, numbers to two decimals', () => {
    expect(recordDrawLists(sample(), ['data', 'labels'])).toBe(
      [
        'size 200 40',
        'pass data',
        'clip 10 0 180 40',
        'rect 10 4 5.5 16 #ff0000ff',
        'rect 20 4 60 12 #7a8fa6ff squircle 4 left',
        'glyphs "10px mono" 15.33 #2d271eff M 12.75, E 18.25',
        'image 10 30 180 4 3x1 ' + recordDrawLists(sample(), ['data']).split('3x1 ')[1].split('\n')[0],
        'unclip',
        'pass labels',
        'text "A <long> label" 2 12 start "12px sans" #6b655cff width 80 max 40',
        '',
      ].join('\n'),
    );
  });
});

describe('the SVG backend', () => {
  it('writes passes as groups, clips, squircles as the core’s path, images as runs, fades past a max width', () => {
    const svg = drawListsToSvg(sample(), { metadata: { 'seq-settled': 'true', 'Bad Key': 'x' } });
    expect(svg.startsWith('<svg xmlns="http://www.w3.org/2000/svg" width="200" height="40" viewBox="0 0 200 40"')).toBe(
      true,
    );
    expect(svg).toContain('data-seq-settled="true"');
    expect(svg).not.toContain('Bad Key');
    expect(svg).toContain('<g data-pass="data"><g clip-path="url(#c0)">');
    expect(svg).toContain('<rect x="10" y="4" width="5.5" height="16" fill="#ff0000"/>');
    // the left end round (a cubic), the right square
    expect(svg).toMatch(/<path d="M 24 4 L 80 4 L 80 16 L 24 16 C [^"]+ Z" fill="#7a8fa6"\/>/);
    expect(svg).toContain('<text x="12.75 18.25" y="15.33" text-anchor="middle"');
    // the image's two red pixels are one run; the clear one is not drawn
    expect(svg).toContain('<rect x="10" y="30" width="120" height="4" fill="#ff0000"/>');
    expect(svg).toContain('mask="url(#f1)"');
    // ids take a prefix, so two documents inlined in one page do not collide
    expect(drawListsToSvg(sample(), { idPrefix: 'a-' })).toContain('clip-path="url(#a-c0)"');
    expect(svg).toContain('A &lt;long&gt; label');
    expect(svg.endsWith('</svg>')).toBe(true);
  });

  it('escapes text and attributes; a character XML forbids becomes U+FFFD, one for one', () => {
    expect(escapeXml(`<script>alert("x")</script> & 'q'\u0001`)).toBe(
      '&lt;script&gt;alert(&quot;x&quot;)&lt;/script&gt; &amp; &#39;q&#39;�',
    );
  });

  it('a spec’s labels and data never become markup', async () => {
    const b = createSeqSpec({ title: '<svg onload="x">' });
    const s = b.data.sequence('s', { alphabet: 'protein', letters: 'MEEPQSDPSV' });
    const r = b.data.ranges('r', {
      target: s,
      rows: {
        start: new Int32Array([0]),
        end: new Int32Array([10]),
        id: categoryColumn(['<img src=x onerror=alert(1)>']),
        label: ['</text><script>alert(1)</script>'],
      },
    });
    b.views.view('v', (v) => {
      const sec = v.section('s', { of: s });
      sec.track('t', { label: '<b onclick="x">&</b>' }).layer({ type: 'ranges', data: r });
    });
    const out = await renderSvg(b.build(), { size: { width: 600, height: 60 } });
    expect(out.svg).not.toMatch(/<script|<img|<b |onclick="/);
    expect(out.svg).toContain('&lt;/text&gt;&lt;script&gt;alert(1)&lt;/script&gt;');
    expect(out.svg).toContain('&lt;b onclick=&quot;x&quot;&gt;&amp;&lt;/b&gt;');
  });
});

describe('tokens, text and shapes', () => {
  it('reads hex colors and tokens; ink follows its fill', () => {
    expect(parseColor('#f00')).toBe(0xff0000ff);
    expect(parseColor('#ff000080')).toBe(0xff000080);
    expect(parseColor('red')).toBeUndefined();
    expect(resolveColor('token:accent', LightTokens)).toBe(LightTokens.colors['msui-accent']);
    expect(resolveColor('token:range', LightTokens)).toBe(LightTokens.colors['seq-range']);
    expect(resolveColor('token:nope', LightTokens)).toBeUndefined();
    expect(inkOn(0xffffffff, LightTokens)).toBe(LightTokens.colors['seq-ink-dark']);
    expect(inkOn(0x000000ff, LightTokens)).toBe(LightTokens.colors['seq-ink-light']);
    expect(luminance(0xffffffff)).toBeCloseTo(1);
  });

  it('measures with the bundled table: an Arial-compatible sans, a fixed-advance mono', () => {
    expect(bundledWidth('10px mono', 'MEEP')).toBeCloseTo(24);
    expect(bundledWidth('12px sans', 'Domains')).toBeCloseTo((722 + 556 + 833 + 556 + 222 + 556 + 500) * 0.012);
    const m = bundledTextMetrics();
    expect(m.measure('12px sans', 'Domains')).toBe(bundledWidth('12px sans', 'Domains'));
  });

  it('squircles: a third of the height, capped, square when small; the handle of a quartic superellipse', () => {
    expect(squircleRadius(100, 12, 1 / 3, 8, 6)).toBe(4);
    expect(squircleRadius(100, 30, 1 / 3, 8, 6)).toBe(8);
    expect(squircleRadius(5, 12, 1 / 3, 8, 6)).toBe(0);
    // the cubic's midpoint on |x|⁴ + |y|⁴ = 1
    const a = SquircleHandle;
    const mid = (0 + 3 * a + 3 + 1) / 8;
    expect(2 * mid ** 4).toBeCloseTo(1, 6);
  });
});
