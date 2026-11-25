// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The card's overflow and the enterable card (_spec/README.md "Ports", Hover card; "Behaviors", hoverCard): the default
// adapter's box is capped and bounds its value column, a long value wraps there (break opportunities after `,` `;` `|`
// `:`) to `hover.tooltipLines`, then an ellipsis; a long VCF REF is titled by its first letters and length, not a hash;
// a card with a link takes the pointer and is entered across the gap (EC-261), one without never takes it. Canvas-free:
// the card's box as styled (jsdom lays nothing out).

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { BridgeMs, HoverCardProvider } from '../behaviors/hover-card.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { type Allele, alleleTitle } from '../describe-allele.js';
import type { HoverCard, HoverCardPlace, HoverCardShowOptions } from '../hover-card.js';
import { breakable, CardMaxPx, hoverCardOverlay } from '../mount-hover-card.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

let frames: Frames;
const viewers: SeqViewer[] = [];
let clock = 1000;

beforeEach(() => {
  frames = installFrames();
  clock = 1000;
  vi.spyOn(DomHost, 'now').mockImplementation(() => clock);
});

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames.restore();
  vi.restoreAllMocks();
  document.body.innerHTML = '';
});

const place: HoverCardPlace = {
  anchor: { x: 10, y: 10, width: 8, height: 8 },
  pointer: { x: 14, y: 14 },
  view: { x: 0, y: 0, width: 800, height: 400 },
  side: 'below',
  x: 14,
  y: 26,
};

const LongValue =
  'MONDO:MONDO:0018875,MedGen:C0085390,OMIM:151623,Orphanet:524,SNOMED_CT:428850001,Human_Phenotype_Ontology';

describe('the card’s overflow', () => {
  it('a capped box whose value column wraps a long value to `lines`, then an ellipsis', () => {
    const host = document.createElement('div');
    document.body.append(host);
    const port = hoverCardOverlay()(host);
    const card: HoverCard = { title: 'p.R175H', rows: [{ label: 'CLNDISDB', value: LongValue }], text: '' };
    port.show(card, place, { lines: 3 });
    const box = host.querySelector<HTMLElement>('.seqx-hover-card')!;
    expect(box.style.maxWidth).toBe(`${CardMaxPx}px`);
    const form = box.querySelector<HTMLElement>('.msui-form')!;
    expect(form.style.gridTemplateColumns).toBe('max-content minmax(0, 1fr)');
    const value = box.querySelector<HTMLElement>('.msui-value')!;
    expect(value.style.getPropertyValue('-webkit-line-clamp')).toBe('3');
    expect(value.style.overflow).toBe('hidden');
    expect(value.style.overflowWrap).toBe('anywhere');
    // a break opportunity after each `,` and `:` (the text read aloud is the card's data, not this)
    expect(value.textContent).toBe(breakable(LongValue));
    expect(value.textContent!.split('​').length).toBeGreaterThan(8);
    port.show(card, place, { lines: 5 });
    expect(value.style.getPropertyValue('-webkit-line-clamp')).toBe('5');
    // shown as a clamped box: the line clamp holds only with its display
    expect(value.style.display).toBe('-webkit-box');
    port.dispose?.();
  });

  it('a long text with no label spans the row and wraps as a value does', () => {
    const host = document.createElement('div');
    document.body.append(host);
    const port = hoverCardOverlay()(host);
    port.show({ title: 'p 4', rows: [{ label: LongValue, value: '' }], text: '' }, place, { lines: 3 });
    const label = host.querySelector<HTMLElement>('.msui-label')!;
    const value = host.querySelector<HTMLElement>('.msui-value')!;
    expect(label.style.gridColumn).toBe('1 / -1');
    expect(label.style.getPropertyValue('-webkit-line-clamp')).toBe('3');
    expect(label.textContent).toBe(breakable(LongValue));
    expect(value.style.display).toBe('none');
    // the next card's row with a value: the columns again
    port.show({ title: 'p 4', rows: [{ label: 'Note', value: 'x' }], text: '' }, place, { lines: 3 });
    expect(label.style.gridColumn).toBe('');
    expect(value.style.display).toBe('-webkit-box');
    port.dispose?.();
  });

  it('a long VCF REF titled by its first letters and its length, not its hash', () => {
    const ref = `G${'AGGCTT'.repeat(78)}`;
    const allele = { ref, alt: 'G' } as Allele;
    expect(alleleTitle('7674964:469:edaeae47>G', allele)).toBe(`7674964:GAGGCT… ${ref.length} bp>G`);
    expect(alleleTitle('7674964:A>G', { ref: 'A', alt: 'G' } as Allele)).toBe('7674964:A>G');
  });
});

describe('the enterable card (EC-261)', () => {
  it('the adapter: a linked card takes the pointer and runs its link; a plain one never takes it', () => {
    const host = document.createElement('div');
    document.body.append(host);
    const port = hoverCardOverlay()(host);
    const plain: HoverCard = { title: 'x', rows: [{ label: 'Note', value: 'text' }], text: '' };
    port.show(plain, place, { lines: 3 });
    const box = host.querySelector<HTMLElement>('.seqx-hover-card')!;
    expect(box.style.pointerEvents).toBe('none');
    const ran: string[] = [];
    const entered: string[] = [];
    const linked: HoverCard = {
      title: 'x',
      rows: [
        {
          label: 'Note',
          value: 'see',
          rich: [{ kind: 'link', href: '!seq-focus=/refs/a', children: [{ kind: 'text', text: 'see' }] }],
        },
      ],
      text: '',
    };
    const opts: HoverCardShowOptions = {
      lines: 3,
      link: (h) => void ran.push(h),
      enter: () => void entered.push('enter'),
      leave: () => void entered.push('leave'),
    };
    port.show(linked, place, opts);
    expect(box.style.pointerEvents).toBe('auto');
    box.dispatchEvent(new Event('pointerenter'));
    box.querySelector<HTMLAnchorElement>('a[data-seq-link]')!.click();
    box.dispatchEvent(new Event('pointerleave'));
    expect(ran).toEqual(['!seq-focus=/refs/a']);
    expect(entered).toEqual(['enter', 'leave']);
    port.hide();
    expect(box.style.pointerEvents).toBe('none');
    port.dispose?.();
  });

  it('the behavior: a linked card waits the grace for the pointer, stays while entered, goes on its leave', async () => {
    const calls: ({ card: HoverCard; options?: HoverCardShowOptions } | 'hide')[] = [];
    const base = DefaultSeqViewerSpec();
    const v = createSeqViewer({
      ...base,
      ui: undefined,
      config: [[viewerItem('viewport.maxPxPerPosition'), 1000]] as never,
      behaviors: base.behaviors.map((b) =>
        b.provider === HoverCardProvider
          ? {
              provider: HoverCardProvider,
              params: {
                port: () => ({
                  show: (card: HoverCard, _: HoverCardPlace, options?: HoverCardShowOptions) =>
                    void calls.push({ card, ...(options ? { options } : {}) }),
                  hide: () => void calls.push('hide'),
                }),
              },
            }
          : b,
      ),
    });
    viewers.push(v);
    const b = createSeqSpec();
    const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLL' });
    b.refs.ref('a', { to: `${p}/(range :from 1 :to 3)` });
    b.views.tooltip('t', {
      target: p,
      rows: { start: [2], end: [6], tooltip: ['x'] },
      text: '{} [see](!seq-focus=/refs/a)',
    });
    b.views.view('v', (vw) => vw.section('s', { of: p }).track('t').layer('l', { type: 'letters' }));
    await v.load(b.build());
    v.resize(900);
    await v.settle();
    frames.flush();
    v.draw();
    const s = v.frame!.layout.sections[0];
    const t = s.tracks[0].box;
    const at = (pos: number) => [(s.axis!.x(pos) + s.axis!.x(pos + 1)) / 2, t.y + t.height / 2] as const;
    v.input.move(...at(3));
    clock += 500;
    frames.timers.fire();
    const shown = calls.at(-1);
    if (!shown || shown === 'hide') throw new Error('no card');
    expect(shown.options?.link).toBeDefined();
    // the pointer leaves the anchor: the card stays the grace
    v.input.move(...at(20));
    expect(calls.at(-1)).not.toBe('hide');
    // entered within it: it stays past the grace, the hover frozen
    shown.options!.enter!();
    clock += BridgeMs + 10;
    frames.timers.fire();
    expect(calls.at(-1)).not.toBe('hide');
    // its leave hides it
    shown.options!.leave!();
    expect(calls.at(-1)).toBe('hide');
    // not entered: the grace ends and the card goes
    v.input.move(...at(3));
    clock += 500;
    frames.timers.fire();
    expect(calls.at(-1)).not.toBe('hide');
    v.input.move(...at(20));
    expect(calls.at(-1)).not.toBe('hide');
    clock += BridgeMs + 10;
    frames.timers.fire();
    expect(calls.at(-1)).toBe('hide');
  });
});
