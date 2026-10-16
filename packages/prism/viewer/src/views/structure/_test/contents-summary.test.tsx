/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { render } from '@molstar/prism-ui/dom/render';
import { afterEach, describe, expect, it } from 'vitest';
import type { ContentsSummary } from '../../../model/scene.js';
import { ContentsSummaryView } from '../contents-summary.js';

// A structure's contents in its inspector: a mesoscale structure's many entities are not drawn all at once.

afterEach(() => {
  document.body.innerHTML = '';
});

function mount(count: number) {
  const summary: ContentsSummary = {
    counts: [{ label: 'Entities', value: String(count) }],
    entities: Array.from({ length: count }, (_, i) => ({
      id: `${i + 1}`,
      name: `PROTEIN ${i + 1}`,
      type: 'polymer',
      chains: ['A'],
      atoms: 100,
    })),
    reveal: () => {},
  };
  const container = document.createElement('div');
  document.body.append(container);
  render(() => <ContentsSummaryView summary={summary} />, container);
  const entities = () => container.querySelectorAll('[role=listitem]').length;
  const more = () => [...container.querySelectorAll('button')].find((b) => b.textContent?.startsWith('Show all'));
  return { entities, more };
}

describe('the contents summary', () => {
  it('shows the first few entities, then all on Show all', () => {
    const t = mount(400);
    expect(t.entities()).toBe(5);
    expect(t.more()?.textContent).toContain('Show all 400');
    expect(t.more()?.textContent).toContain('395 more');
    t.more()!.click();
    expect(t.entities()).toBe(400);
    expect(t.more()).toBeUndefined();
  });

  it('shows a few entities all at once', () => {
    const t = mount(6);
    expect(t.entities()).toBe(6);
    expect(t.more()).toBeUndefined();
  });
});
