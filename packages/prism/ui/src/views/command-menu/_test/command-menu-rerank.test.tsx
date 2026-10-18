/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { CommandManager } from '../../../commands/manager.js';
import { render } from '../../../dom/render.js';
import { flushFrame } from '../../../dom/scope.js';
import { observableState } from '../../../reactive/observable-state.js';
import { CommandMenu } from '../command-menu.js';

const tick = () => new Promise((r) => setTimeout(r, 0));

// A re-ranked search with the same length and ends must show the new middle rows (the highlight and Enter use them).
describe('command menu re-ranking', () => {
  it('shows the middle rows in the order the query ranks them', async () => {
    const manager = new CommandManager<object>({
      context: observableState<object>({}).asReadonly,
      formats: { get: () => undefined },
      preview: { set() {}, clear() {} },
      browseDelay: 0,
    });
    const titles = ['Color', 'Collapse', 'Outline color', 'Cool', 'Clear color', 'Recolor'];
    manager.register(titles.map((t) => ({ id: t, title: t, category: 'X', run: () => ({}) })));
    const container = document.createElement('div');
    document.body.append(container);
    const dispose = render(() => <CommandMenu manager={manager} onClose={() => {}} />, container);
    const input = container.querySelector('input') as HTMLInputElement;
    const shown = () => [...container.querySelectorAll('[role=option]')].map((r) => r.textContent ?? '');
    const order = (t: string) => manager.query(t).matches.map((m) => m.command.title);
    let reordered = 0;
    let prev: string[] = [];
    for (const t of ['c', 'co', 'col', 'colo', 'color']) {
      input.value = t;
      input.dispatchEvent(new Event('input'));
      flushFrame();
      await tick();
      const want = order(t);
      const got = shown().map((s) => want.find((w) => s.startsWith(w)) ?? s);
      expect(got).toEqual(want);
      if (
        prev.length === want.length &&
        prev[0] === want[0] &&
        prev.at(-1) === want.at(-1) &&
        prev.join() !== want.join()
      )
        reordered++;
      prev = want;
    }
    expect(reordered).toBeGreaterThan(0);
    dispose();
  });
});
