/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { render } from '../../../dom/render.js';
import { Show } from '../../../dom/control-flow.js';
import { observableState } from '../../../reactive/observable-state.js';
import { Popover } from '../popover.js';

const flush = () => new Promise((r) => setTimeout(r, 0));
const press = (el: Element) => el.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, composed: true }));

function mount() {
  const open = observableState(true);
  const onClose = vi.fn(() => open.set(false));
  const anchor = document.createElement('button');
  // another floating layer (a select's options, opened from inside the panel)
  const layer = document.createElement('div');
  layer.setAttribute('popover', 'manual');
  const option = document.createElement('button');
  layer.append(option);
  const outside = document.createElement('div');
  const container = document.createElement('div');
  document.body.append(anchor, container, layer, outside);
  render(
    () => (
      <Show when={open.asReadonly}>
        {() => (
          <Popover anchor={anchor} onClose={onClose} title="Scene defaults">
            <button type="button">Inside</button>
          </Popover>
        )}
      </Show>
    ),
    container,
  );
  return { open, onClose, anchor, option, outside, container };
}

describe('Popover', () => {
  it('is a dialog named by its title; a press inside it, on its button or in another layer keeps it', async () => {
    const { onClose, anchor, option } = mount();
    await flush();
    const dialog = document.querySelector('[role=dialog]')!;
    expect(dialog.getAttribute('aria-label')).toBe('Scene defaults');
    press(dialog.querySelector('button')!);
    press(anchor);
    press(option);
    expect(onClose).not.toHaveBeenCalled();
  });

  it('closes on a press outside, and on Escape with the focus back', async () => {
    const a = mount();
    await flush();
    press(a.outside);
    expect(a.onClose).toHaveBeenCalledWith();

    const b = mount();
    await flush();
    document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
    expect(b.onClose).toHaveBeenCalledWith(true);
  });

  it('leaves Escape to a layer it opened while the focus is there', async () => {
    const { onClose, option } = mount();
    await flush();
    option.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true, composed: true }));
    expect(onClose).not.toHaveBeenCalled();
  });
});
