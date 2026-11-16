/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { describe, expect, it } from 'vitest';
import { flush, render } from '../../../dom/index.js';
import { Loaded, Loading, Skeleton, SkeletonRows, SkeletonText } from '../skeleton.js';

function mount(content: () => any) {
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(content, container);
  return { container, dispose };
}

describe('Skeleton', () => {
  it('a loading region says so once; its shapes are hidden from assistive technology', () => {
    const { container, dispose } = mount(() => (
      <Loading label="Loading components">
        <SkeletonRows count={3} />
        <SkeletonText lines={2} />
        <Skeleton height="64px" />
      </Loading>
    ));
    const status = container.querySelector('[role="status"]')!;
    expect(status.getAttribute('aria-label')).toBe('Loading components');
    expect(status.getAttribute('aria-busy')).toBe('true');
    // Every child of the region is hidden: 3 rows, 2 lines, 1 block.
    expect([...status.children].map((c) => c.getAttribute('aria-hidden'))).toEqual(Array(6).fill('true'));
    dispose();
  });

  it('rows sit on the row primitive, with or without the icon slot; text lines end with a shorter one', () => {
    const { container, dispose } = mount(() => [
      <SkeletonRows count={2} />,
      <SkeletonRows count={1} icon={false} />,
      <SkeletonText lines={3} size="md" />,
    ]);
    const rows = [...container.querySelectorAll('.msui-row')];
    expect(rows).toHaveLength(3);
    expect(rows.map((r) => r.querySelector('.msui-slot') !== null)).toEqual([true, true, false]);
    expect(rows[2].classList.contains('msui-row-plain')).toBe(true);
    const lines = [...container.querySelectorAll('[data-size="md"]')];
    expect(lines).toHaveLength(3);
    expect((lines[2].firstElementChild as HTMLElement).style.width).toBe('40%');
    dispose();
  });

  it('a block takes its size and corners', () => {
    const { container, dispose } = mount(() => [
      <Skeleton width="32px" height="32px" round />,
      <Skeleton height="80px" radius="lg" />,
    ]);
    const [dot, block] = [...container.querySelectorAll<HTMLElement>('[aria-hidden="true"]')];
    expect(dot.style.borderRadius).toBe('50%');
    expect(block.style.height).toBe('80px');
    expect(block.style.borderRadius).toBe('var(--msui-radius-lg)');
    dispose();
  });
});

describe('Loaded', () => {
  it('shows the placeholder, then the content in the same place; content there at first shows at once', async () => {
    const loaded = new State<string | undefined>(undefined);
    const { container, dispose } = mount(() => (
      <Loaded when={loaded} label="Loading the summary" placeholder={() => <SkeletonText lines={2} />}>
        {(text) => <p>{text}</p>}
      </Loaded>
    ));
    expect(container.querySelector('[role="status"]')?.getAttribute('aria-label')).toBe('Loading the summary');
    expect(container.querySelector('p')).toBeNull();
    loaded.set('1TQN');
    await flush();
    // No Web Animations here: the placeholder goes at once; the content enters with its fade.
    expect(container.querySelector('[role="status"]')).toBeNull();
    const layer = container.querySelector('p')!.parentElement!;
    expect(layer.getAttribute('data-reveal')).toBe('fade');
    expect(layer.parentElement!.children).toHaveLength(1);
    dispose();

    const ready = mount(() => (
      <Loaded when="ready" placeholder={() => <SkeletonText />}>
        {(text) => <p>{text}</p>}
      </Loaded>
    ));
    expect(ready.container.querySelector('[role="status"]')).toBeNull();
    expect(ready.container.querySelector('p')!.parentElement!.hasAttribute('data-reveal')).toBe(false);
    ready.dispose();
  });
});
