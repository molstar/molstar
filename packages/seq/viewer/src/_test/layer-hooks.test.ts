// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A layer type's words (core's `LayerType.describe`, `summary`, `announce`; _spec/README.md "Hover card", Providers,
// `seq:layer`): an extension's layer says what it draws: its lines head the hover card (the viewer's own text after),
// its summary is its layer's ARIA group, and what it adds at the keyboard cursor's position is said by the live region.
// A binding scoped to its type (`GestureBinding.layer`) wins over the grammar's double-click on its hits.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { type LayerType, RangesType } from '@molstar/seqx-core/render';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { Announcer } from '../announcer.js';
import { bindingTable } from '../bindings.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { describerFor } from '../describe.js';
import { defineViewerExtension } from '../extension.js';
import { layerAria, setLayerAria } from '../layer-aria.js';
import { createSeqViewer, type SeqViewer, type ViewerCommand } from '../viewer.js';
import { stubCanvasContext } from './canvas-stub.js';
import { type Frames, installFrames } from './frames.js';

let frames: Frames;
const viewers: SeqViewer[] = [];

beforeEach(() => {
  frames = installFrames();
  vi.spyOn(DomHost, 'now').mockImplementation(() => 1000);
  stubCanvasContext();
});

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames.restore();
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
  document.body.innerHTML = '';
});

const Section = '/views/v/s';

/** Ranges drawn as bars that say what they are. */
const Bars: LayerType = {
  ...(RangesType as LayerType),
  name: 'demo:bars',
  shape: 'ranges',
  fallback: 'ranges',
  describe: (hit) => (hit.element ? [`Bar ${hit.element.key}`, 'Links: 3'] : undefined),
  summary: () => '1 bar in view',
  announce: (_lay, _system, position) => (position >= 2 && position < 8 ? 'stem of 3 links' : undefined),
};

function spec() {
  const b = createSeqSpec({ title: 'words', extensions: [{ id: 'demo', version: '0.1.0' }] });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFS' });
  const d = b.data.ranges('d', { target: p, rows: { start: [2], end: [8], id: ['A'], label: ['Domain A'] } });
  b.views.view('v', (v) => {
    const s = v.section('s', { of: p });
    s.track('letters', { label: 'p' }).layer('letters', { type: 'letters' });
    s.track('bars', { label: 'Bars' }).layer('b', { type: 'demo:bars', fallback: 'ranges', data: d } as never);
  });
  return b.build();
}

async function make() {
  const opened: unknown[] = [];
  const open = {
    id: 'demo:open',
    title: 'Open',
    category: 'Demo',
    run: (_args: unknown, context: unknown) => (opened.push(context), { message: 'opened' }),
  } as unknown as ViewerCommand;
  const ext = defineViewerExtension({
    id: 'demo',
    version: '0.1.0',
    requires: '0.10',
    registry: [{ kind: 'layer-type', type: Bars }],
    commands: [open],
    bindings: [{ gesture: 'double-click', target: 'item', action: 'demo:open', layer: 'demo:bars' }],
  });
  const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined, extensions: [ext] });
  viewers.push(v);
  await v.load(spec());
  v.resize(600);
  await v.settle();
  frames.flush();
  v.draw();
  return { v, opened };
}

const bars = (v: SeqViewer) => v.frame!.layout.sections[0].tracks[1];

describe('a layer type’s words', () => {
  it('its `describe` heads the hover card, the viewer’s own text after', async () => {
    const { v } = await make();
    const l = bars(v).layers[0];
    const d = describerFor({ ...v, layout: () => v.frame?.layout });
    const card = d.describe({
      section: Section,
      layer: l.key,
      system: '/data/p',
      position: 4,
      element: { owner: '/data/d', key: 'A' },
    })!;
    expect(card.title).toBe('Bar A');
    expect(card.rows[0]).toEqual({ label: 'Links', value: '3' });
    // the viewer's item rows after the layer's
    expect(card.rows.some((r) => r.label === 'Extent')).toBe(true);
    expect(card.text.startsWith('Bar A, Links: 3, ')).toBe(true);
  });

  it('its `announce` is said at the keyboard cursor, in the live region', async () => {
    const { v } = await make();
    const a = new Announcer(v, { minGapMs: 0 });
    for (let i = 0; i < 4; i++) await v.execute('cursor.next', { section: Section });
    expect(a.message.value.text).toContain('stem of 3 links');
    a.dispose();
  });

  it('its `summary` is its layer’s ARIA group, set when the frame settles', async () => {
    const { v } = await make();
    expect(layerAria(v.frame!.layout)).toEqual([{ label: 'Bars', description: '1 bar in view' }]);
    const host = document.createElement('div');
    setLayerAria(host, v.frame!.layout);
    const g = host.querySelector('[role="group"]')!;
    expect([g.getAttribute('aria-label'), g.textContent]).toEqual(['Bars', '1 bar in view']);
  });
});

describe('a binding scoped to a layer type', () => {
  it('resolves on that type’s hits only, beside the grammar’s', () => {
    const t = bindingTable(
      undefined,
      [[{ gesture: 'double-click', target: 'item', action: 'view.fit', layer: 'demo:bars' }]],
      {},
    );
    expect(t.resolveIn('demo:bars', 'double-click', {}, 'item')).toBe('view.fit');
    expect(t.resolveIn('ranges', 'double-click', {}, 'item')).toBeUndefined();
    // the grammar's own binding stays for every other layer
    expect(t.resolve('double-click', {}, 'item')).toBe('view.focus');
  });

  it('a double-click on its hit runs its command with the hit, never the focus', async () => {
    const { v, opened } = await make();
    const t = bars(v);
    const x = t.box.x + ((4.5 - 0) / 20) * t.box.width;
    const y = t.box.y + t.box.height / 2;
    const hit = v.hit(x, y);
    expect(hit?.element?.key).toBe('A');
    const before = v.view.get(Section);
    expect(v.input.doubleClick(x, y)).toBe(true);
    await v.settle();
    expect(opened.length).toBe(1);
    expect((opened[0] as { hit: { element: { key: string } } }).hit.element.key).toBe('A');
    expect(v.view.get(Section)).toEqual(before);
  });
});
