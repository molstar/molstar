/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { DomHost, render } from '../../dom/index.js';
import { BandPlaceholder, StagePlaceholder } from '../placeholders.js';

// The stand-ins (stage §6): the canvas's atoms and the band's residues feed hover, the band's column once a frame.

const host = DomHost.requestAnimationFrame;
let queue: ((t: number) => void)[] = [];
const frame = () => {
  const run = queue;
  queue = [];
  for (const f of run) f(0);
};
beforeEach(() => {
  DomHost.requestAnimationFrame = (f) => {
    queue.push(f);
    return 0;
  };
});
afterEach(() => {
  document.body.innerHTML = '';
  frame();
  DomHost.requestAnimationFrame = host;
  queue = [];
});

function mount(content: () => Node) {
  const container = document.createElement('div');
  document.body.append(container);
  render(content, container);
  return container;
}

describe('placeholders', () => {
  it('the canvas: the residue under the pointer, once per residue, none when it leaves', () => {
    const hover = vi.fn();
    const el = mount(() => <StagePlaceholder onHover={hover} />);
    // Two atoms of residue 7: one call.
    const [a, b] = el.querySelectorAll('[data-index="7"]');
    a.dispatchEvent(new PointerEvent('pointerover', { bubbles: true }));
    b.dispatchEvent(new PointerEvent('pointerover', { bubbles: true }));
    expect(hover.mock.calls).toEqual([[7]]);
    el.firstElementChild!.dispatchEvent(new PointerEvent('pointerleave'));
    expect(hover).toHaveBeenLastCalledWith(undefined);
    // Three chains of 70 residues, four atoms each.
    const atoms = [...el.querySelectorAll('[data-index]')];
    expect(atoms.length).toBe(840);
    expect(new Set(atoms.map((x) => x.getAttribute('data-index'))).size).toBe(210);
  });

  it('an empty canvas draws nothing', () => {
    const el = mount(() => <StagePlaceholder empty />);
    expect(el.querySelector('svg')).toBeNull();
  });

  it('the band: rows of residues under a ruler, a selection, the hovered column in a frame', () => {
    const hover = vi.fn();
    const el = mount(() => <BandPlaceholder rows={3} length={40} selection={[2, 5]} onHover={hover} />);
    const svg = el.querySelector('svg')!;
    vi.spyOn(svg, 'getBoundingClientRect').mockReturnValue({ left: 0, top: 0 } as DOMRect);
    expect(el.querySelectorAll('rect[data-kind]').length).toBe(120);
    expect(el.querySelector('[data-msui-hover-source]')).not.toBeNull();
    svg.dispatchEvent(new PointerEvent('pointerenter'));
    // Column 3 (12 px a residue), the second row (a 16 px ruler, 16 px rows).
    svg.dispatchEvent(new PointerEvent('pointermove', { clientX: 40, clientY: 40 }));
    svg.dispatchEvent(new PointerEvent('pointermove', { clientX: 41, clientY: 40 }));
    expect(hover).not.toHaveBeenCalled();
    frame();
    expect(hover.mock.calls).toEqual([[{ row: 1, index: 3 }]]);
    const column = el.querySelector('[data-shown]')!;
    expect(column.getAttribute('transform')).toBe('translate(36 0)');
    svg.dispatchEvent(new PointerEvent('pointerleave'));
    frame();
    expect(hover).toHaveBeenLastCalledWith(undefined);
    expect(el.querySelector('[data-shown]')).toBeNull();
  });
});
