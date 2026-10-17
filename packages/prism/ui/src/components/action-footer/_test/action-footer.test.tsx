/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @vitest-environment jsdom

import { afterEach, describe, expect, it, vi } from 'vitest';
import { render } from '../../../dom/render.js';
import { observableState } from '../../../reactive/observable-state.js';
import { ActionFooter } from '../action-footer.js';

// Enter in a field of an action's block applies it, as Mol* v5's actions did (a PDB id typed, Enter: it loads).

afterEach(() => {
  document.body.innerHTML = '';
});

function mount(disabled = false) {
  const host = document.createElement('div');
  document.body.append(host);
  const apply = vi.fn();
  const committed: string[] = [];
  const off = observableState(disabled);
  render(
    () => (
      <div>
        <input type="text" class="id" onChange={(e) => committed.push((e.currentTarget as HTMLInputElement).value)} />
        <input type="number" class="n" />
        <textarea class="notes" />
        <input type="checkbox" class="check" />
        <ActionFooter icon="download" label="Download" disabled={off.asReadonly} onApply={apply} />
      </div>
    ),
    host,
  );
  const key = (sel: string, init: KeyboardEventInit = {}) => {
    const el = host.querySelector<HTMLElement>(sel)!;
    const e = new KeyboardEvent('keydown', { key: 'Enter', bubbles: true, cancelable: true, ...init });
    el.dispatchEvent(e);
    return e;
  };
  return { host, apply, committed, key, off };
}

const tick = () => new Promise((r) => setTimeout(r, 0));

describe('the action footer', () => {
  it('applies on Enter in a text or number field, after the field commits its value', async () => {
    const t = mount();
    await tick();
    (t.host.querySelector('.id') as HTMLInputElement).value = '1grm';
    const e = t.key('.id');
    expect(e.defaultPrevented).toBe(true);
    expect(t.committed).toEqual(['1grm']);
    await tick();
    expect(t.apply).toHaveBeenCalledOnce();
    t.key('.n');
    await tick();
    expect(t.apply).toHaveBeenCalledTimes(2);
  });

  it('leaves Enter alone in a text area, a check box, while composing, with a modifier, or when it cannot apply', async () => {
    const t = mount();
    await tick();
    t.key('.notes');
    t.key('.check');
    t.key('.id', { isComposing: true });
    t.key('.id', { ctrlKey: true });
    await tick();
    expect(t.apply).not.toHaveBeenCalled();
    t.off.set(true);
    t.key('.id');
    await tick();
    expect(t.apply).not.toHaveBeenCalled();
  });
});
