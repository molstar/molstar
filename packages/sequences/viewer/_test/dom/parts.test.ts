/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The DOM parts' content (UX2; seqspec:interaction#the-default-ui, seqspec:interaction#the-legend): the legend's "+n
// more" and ramps with their captions (UX1's helpers drawn), the default Description as plain text only and a host's
// part for markdown, and the element on the real canvas surface with real pointer and keyboard events.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { afterEach, describe, expect, it } from 'vitest';
import * as D from '../../src/dom/index.js';
import '../../src/skin.css';
import {
  createSeqViewer,
  KernelRenderers,
  type LayerRenderer,
  manualScheduler,
  type SeqViewerContext,
} from '../../src/index.js';
import { fixtureMeasurer } from '../../src/testing/index.js';
import { flush, headlessSpec, stubSurface, uiDoc, WIDTH, xOf } from './fixtures.js';

const cleanups: (() => void)[] = [];
afterEach(() => {
  for (const c of cleanups.splice(0)) c();
});
const A = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ';

/** A context showing `doc`, the default spec with some kernel renderers replaced by name, mounted in a host. */
async function mount(
  doc: ReturnType<ReturnType<typeof createSeqViewSpec>['getDocument']>,
  replace: Record<string, LayerRenderer> = {},
) {
  const spec = headlessSpec();
  const [first, ...rest] = spec.registry ?? [];
  const registry = [
    { ...first, renderers: (KernelRenderers as LayerRenderer[]).map((r) => replace[r.name] ?? r) },
    ...rest,
  ];
  const ctx = createSeqViewer({ ...spec, registry }, { scheduler: manualScheduler(), text: fixtureMeasurer() });
  const done = ctx.setDocument(doc);
  for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(WIDTH, 220, 1);
  await done;
  const host = document.createElement('div');
  host.style.width = `${WIDTH}px`;
  document.body.appendChild(host);
  const m = D.mountSeqViewer(host, ctx, { surface: stubSurface });
  cleanups.push(() => {
    m.dispose();
    host.remove();
    ctx.dispose();
  });
  await flush();
  return { ctx, host };
}

/** `n` ranges, a category each. */
function manyDoc(n: number) {
  const b = createSeqViewSpec({ title: 'Many' });
  const a = b.data.sequence({ ref: 'a', label: 'A', value: A });
  const start = Array.from({ length: n }, (_, k) => k);
  const d = b.data.annotation({
    ref: 'd',
    geometry: 'ranges',
    target: a,
    value: { start, stop: start.map((x) => x + 1), label: start.map((k) => `cat${k}`) },
  });
  const v = b.view({ ref: 'v', label: 'V' });
  v.block({ sequence: a });
  v.track({ ref: 'v.t', label: 'Many' }).ranges({ data: d });
  return b.getDocument();
}
function heatDoc() {
  const b = createSeqViewSpec({ title: 'Heat' });
  const a = b.data.sequence({ ref: 'a', label: 'A', value: A });
  const d = b.data.annotation({
    ref: 'd',
    geometry: 'values',
    target: a,
    value: { index: [0, 5, 9], value: [0, 5, 20] },
  });
  b.view({ ref: 'v', label: 'V' }).block({ sequence: a }).track({ label: 'Heat' }).color({ data: d });
  return b.getDocument();
}

const inputs = (ramp: HTMLElement) =>
  [...ramp.querySelectorAll<HTMLInputElement>('.seqv-range-input')].map((i) => i.value);

describe('the legend', () => {
  it('a group shows `UI.LegendItems` items and "+n more"; the control expands and collapses it and keeps focus', async () => {
    const { host, ctx } = await mount(manyDoc(30));
    const items = () => host.querySelectorAll('.seqv-legend-item').length;
    const more = () => host.querySelector('.seqv-legend-more') as HTMLButtonElement;
    expect([items(), more().textContent, more().getAttribute('aria-expanded')]).toEqual([24, '+6 more', 'false']);
    more().focus();
    more().click();
    expect([items(), more().textContent, more().getAttribute('aria-expanded'), [...ctx.ui.legendOpen.value]]).toEqual([
      30,
      'fewer',
      'true',
      ['v.main/v.t'],
    ]);
    expect(document.activeElement).toBe(more());
    more().click();
    expect([items(), more().textContent, document.activeElement === more()]).toEqual([24, '+6 more', true]);
    // The open state is the context's: another UI on it agrees, and the details element's own state is kept.
    ctx.ui.toggleLegend('v.main/v.t');
    expect(items()).toBe(30);
    const legend = host.querySelector('.seqv-legend') as HTMLDetailsElement;
    legend.open = false;
    await flush();
    ctx.ui.toggleLegend('v.main/v.t');
    expect((host.querySelector('.seqv-legend') as HTMLDetailsElement).open).toBe(false);
  });

  it('a group at the cap has no control', async () => {
    const { host } = await mount(manyDoc(24));
    expect(host.querySelectorAll('.seqv-legend-item').length).toBe(24);
    expect(host.querySelector('.seqv-legend-more')).toBeNull();
  });

  it('a continuous palette is a gradient bar with its two ends, not a list; a caption stands beside it', async () => {
    const { host } = await mount(heatDoc());
    const ramp = host.querySelector('.seqv-ramp') as HTMLElement;
    expect(inputs(ramp)).toEqual(['0', '20']);
    expect(host.querySelectorAll('.seqv-legend-item').length).toBe(0);
    expect((ramp.querySelector('.seqv-ramp-bar') as HTMLElement).style.background).toContain('linear-gradient(');
    const captioned: LayerRenderer = {
      ...(KernelRenderers.find((r) => r.name === 'color') as LayerRenderer),
      legend: (_p, _t, sink) =>
        sink.ramp(Float64Array.of(0, 0x112233, 1, 0xffeeaa), 0, 12.3456, 'contacts (number of kinds)'),
    };
    const t = await mount(heatDoc(), { color: captioned });
    const cap = t.host.querySelector('.seqv-ramp') as HTMLElement;
    expect([inputs(cap), cap.textContent]).toEqual([['0', '12.3'], 'contacts (number of kinds)']);
    expect([...cap.children].map((c) => c.className)).toEqual([
      'seqv-range-input',
      'seqv-ramp-track',
      'seqv-range-input',
      'seqv-quiet',
    ]);
  });
});

describe('the color range editor', () => {
  it("sets the ends from the fields, marks the data's span, refuses a bad number, resets", async () => {
    const { host, ctx } = await mount(heatDoc());
    const field = (end: string) =>
      host.querySelector(`.seqv-range-input[data-seqv-key^="range-${end}"]`) as HTMLInputElement;
    expect(host.querySelector('.seqv-ramp-data')).not.toBeNull();
    expect(host.querySelector('[data-seqv-key^="range-reset"]')).toBeNull();
    field('hi').value = '40';
    field('hi').dispatchEvent(new Event('change'));
    await flush();
    expect([...ctx.ui.colorRanges.value.values()]).toEqual([[0, 40]]);
    // Wider than the data: the mark covers half of the bar.
    const mark = host.querySelector('.seqv-ramp-data') as HTMLElement;
    expect([mark.style.left, mark.style.width]).toEqual(['0%', '50%']);
    field('hi').value = '-3';
    field('hi').dispatchEvent(new Event('change'));
    expect(field('hi').value).toBe('40');
    (host.querySelector('[data-seqv-key^="range-reset"]') as HTMLButtonElement).click();
    await flush();
    expect([ctx.ui.colorRanges.value.size, field('hi').value]).toEqual([0, '20']);
  });
});

describe('the default Description', () => {
  it("shows a plaintext description as paragraphs of text only, a snapshot's over the document's", async () => {
    const { host, ctx } = await mount(uiDoc());
    const paragraphs = () => [...host.querySelectorAll('.seqv-description p')].map((p) => p.textContent);
    expect(paragraphs()).toEqual(['The kinase and its domains, as given.', 'A second paragraph.']);
    await ctx.snapshots.apply('s1');
    await flush();
    expect(paragraphs()).toEqual(['Only the first domain is hidden.']);
  });

  it('puts nothing of the text into markup', async () => {
    const doc = {
      ...uiDoc(),
      metadata: {
        description: '<img src=x onerror=alert(1)> **not bold** [x](!seq-select=a)',
        description_format: 'plaintext' as const,
      },
    };
    const { host } = await mount(doc);
    const d = host.querySelector('.seqv-description') as HTMLElement;
    expect(d.children.length).toBe(1);
    expect(d.querySelectorAll('img, b, strong, a, button').length).toBe(0);
    expect(d.textContent).toBe('<img src=x onerror=alert(1)> **not bold** [x](!seq-select=a)');
  });

  it("shows a markdown description not at all; a host's part through the slot shows it with its own renderer", async () => {
    const doc = {
      ...uiDoc(),
      metadata: { description: '# Title\n\nSome **bold** words.', description_format: 'markdown' as const },
    };
    const t = await mount(doc);
    expect(t.host.querySelector('.seqv-description')).toBeNull();
    // The host's part reads the description's source from the context and draws it as it likes.
    const Markdown: D.DomPart<D.PartProps> = (host, ctx) => {
      const at = D.spot(host);
      const subs = D.derived(
        [ctx.document.shown],
        () => ctx.document.shown.value?.metadata?.description ?? '',
        (text) =>
          at.set(
            D.h(
              'div',
              { class: 'host-markdown', 'data-text': text },
              D.h('h1', null, 'Title'),
              'Some ',
              D.h('b', null, 'bold'),
              ' words.',
            ),
          ),
      );
      return { dispose: () => D.off([...subs, at]) };
    };
    const host2 = document.createElement('div');
    document.body.append(host2);
    const m = D.mountSeqViewer(host2, t.ctx, { surface: stubSurface, slots: { Description: Markdown } });
    await flush();
    expect(host2.querySelector('.host-markdown h1')?.textContent).toBe('Title');
    expect(host2.querySelector('.seqv-description')).toBeNull();
    m.dispose();
    host2.remove();
  });
});

describe('the element on the real canvas surface', () => {
  it('mounts canvases, hovers from real pointer events with a tooltip at the pointer, and tells the host once a frame', async () => {
    D.defineSeqViewer();
    const el = document.createElement('seq-viewer');
    el.spec = headlessSpec();
    el.ports = { text: fixtureMeasurer() };
    el.document = uiDoc();
    const host = document.createElement('div');
    host.style.width = `${WIDTH}px`;
    document.body.append(host);
    cleanups.push(() => host.remove());
    host.append(el);
    const hover = new Promise<D.SeqHoverDetail>((r) =>
      el.addEventListener('seq-hover', (e) => r((e as CustomEvent).detail), { once: true }),
    );
    for (let i = 0; i < 20 && !el.shadowRoot?.querySelector('.seqv-surface canvas'); i++)
      await new Promise((r) => setTimeout(r, 25));
    const ctx = el.context as SeqViewerContext;
    const surface = el.shadowRoot?.querySelector('.seqv-surface [tabindex]') as HTMLElement;
    expect(surface.querySelectorAll('canvas').length).toBeGreaterThan(0);
    await flush();
    const layout = ctx.layout.get('chains.main')?.current;
    const rect = surface.getBoundingClientRect();
    const view = ctx.graph.assembled.value?.view;
    expect([!!layout, !!view, rect.width > 100]).toEqual([true, true, true]);
    // The residues track's row, at one of the first cells.
    let y = 0;
    const l = layout as NonNullable<typeof layout>;
    for (let e = 0; e < l.entryCount; e++) {
      const t = l.entryTrack[e];
      if (t >= 0 && view?.set.strings.get(view.resolved.tracks.ref[t]) === 't.res')
        y = l.scrolledTop + l.entryTop[e] + l.entryHeight[e] / 2;
    }
    const x = xOf(ctx, 'chains.main', 10);
    surface.dispatchEvent(
      new PointerEvent('pointermove', {
        clientX: rect.left + x,
        clientY: rect.top + y,
        pointerId: 1,
        isPrimary: true,
        bubbles: true,
      }),
    );
    const detail = await hover;
    expect([detail.owner, detail.track, detail.index]).toEqual(['pointer', 't.res', 10]);
    const tip = el.shadowRoot?.querySelector('.seqv-tooltip') as HTMLElement;
    expect(tip).not.toBeNull();
    const root = el.shadowRoot?.querySelector('.seqv-root') as HTMLElement;
    const at = tip.getBoundingClientRect();
    // At the pointer: right of and below it, in the viewer.
    expect(at.left).toBeGreaterThanOrEqual(rect.left + x);
    expect(at.top).toBeGreaterThanOrEqual(rect.top + y);
    expect(at.right).toBeLessThanOrEqual(root.getBoundingClientRect().right + 1);
    el.remove();
  });
});
