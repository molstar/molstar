/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { render } from '@molstar/prism-ui/dom/render';
import { flush } from '@molstar/prism-ui/dom/scope';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { derive } from '@molstar/prism-ui/reactive/state';
import { afterEach, describe, expect, it, vi } from 'vitest';
import type { LookOverrideSummary, NodeLookModel, SceneLookModel } from '../../../model/look.js';
import { classes } from '../_generated/look-popover.module.css.js';
import { LookButton, NodeLookPopover, SceneLookPopover, sceneOverridesText } from '../look-popover.js';

// The look's popovers and their button (look-popover.tsx) over stub models: the dot and the name telling of overrides,
// what keeps its own values listed, Apply to all and Use default calling the model.

(globalThis as { ResizeObserver?: unknown }).ResizeObserver ??= class {
  observe() {}
  unobserve() {}
  disconnect() {}
};

const Material: LookOverrideSummary = { key: 'material', label: 'Material', names: ['Cartoon of Polymer (1CBS)'] };
const Light: LookOverrideSummary = { key: 'ignoreLight', label: 'Ignore light', names: ['A', 'B'] };
const form = () => new ParamsFormModel({ ignoreLight: PD.Boolean(false, { label: 'Ignore Light' }) });
const settle = () => new Promise((r) => setTimeout(r, 0)).then(() => flush());

let dispose: (() => void) | undefined;
afterEach(() => {
  dispose?.();
  dispose = undefined;
  document.body.innerHTML = '';
});

function mount(view: () => unknown) {
  const box = document.createElement('div');
  document.body.append(box);
  dispose = render(view as never, box);
  flush();
  return box;
}

describe('LookButton', () => {
  it('shows a dot and says so in its name while marked; opens its popover', async () => {
    const marked = observableState<string | undefined>(undefined);
    const box = mount(() => (
      <LookButton
        title="Scene defaults"
        marked={marked.asReadonly}
        popover={(_, close) => (
          <button type="button" data-testid="inside" onClick={() => close()}>
            inside
          </button>
        )}
      />
    ));
    const button = box.querySelector('button')!;
    expect(button.getAttribute('aria-haspopup')).toBe('dialog');
    expect(button.getAttribute('aria-expanded')).toBe('false');
    expect(box.querySelector(`.${classes.dot}`)).toBeNull();
    marked.set('2 representations have their own look');
    flush();
    expect(box.querySelector(`.${classes.dot}`)).not.toBeNull();
    expect(button.getAttribute('aria-label')).toBe('Scene defaults (2 representations have their own look)');

    button.click();
    await settle();
    expect(button.getAttribute('aria-expanded')).toBe('true');
    expect(document.querySelector('[data-testid=inside]')).not.toBeNull();
  });

  it('counts the representations with their own look once each', () => {
    expect(sceneOverridesText([])).toBeUndefined();
    expect(sceneOverridesText([Material])).toBe('a representation has its own look');
    expect(sceneOverridesText([Material, Light])).toBe('3 representations have their own look');
  });
});

describe('the popovers', () => {
  it('Scene defaults lists what is kept as their own; Apply to all writes the default everywhere', async () => {
    const overrides = observableState<readonly LookOverrideSummary[]>([Material]);
    const model: SceneLookModel = { form: form(), overrides: overrides.asReadonly, applyToAll: vi.fn() };
    const anchor = document.createElement('button');
    document.body.append(anchor);
    mount(() => <SceneLookPopover model={model} anchor={anchor} onClose={() => {}} />);
    await settle();
    const dialog = document.querySelector('[role=dialog]')!;
    expect(dialog.getAttribute('aria-label')).toBe('Scene defaults');
    expect(dialog.textContent).toContain('Kept as their own');
    expect(dialog.textContent).toContain('Material');
    [...dialog.querySelectorAll('button')].find((b) => b.textContent === 'Apply to all')!.click();
    expect(model.applyToAll).toHaveBeenCalledWith('material');

    overrides.set([]);
    flush();
    expect(dialog.textContent).not.toContain('Kept as their own');
  });

  it("a node's look offers Use default per value, and all of them back at once", async () => {
    const own = observableState<readonly LookOverrideSummary[]>([Material, Light]);
    const model: NodeLookModel = { form: form(), own: own.asReadonly, reset: vi.fn(), refresh: () => {} };
    const anchor = document.createElement('button');
    document.body.append(anchor);
    const label = observableState('Ligand');
    mount(() => (
      <NodeLookPopover
        model={model}
        title={derive(label.asReadonly, (l) => `Look of ${l}`)}
        anchor={anchor}
        onClose={() => {}}
      />
    ));
    await settle();
    const dialog = document.querySelector('[role=dialog]')!;
    expect(dialog.textContent).toContain('Look of Ligand');
    const buttons = () => [...dialog.querySelectorAll('button')];
    buttons()
      .filter((b) => b.textContent === 'Use default')[0]
      .click();
    expect(model.reset).toHaveBeenCalledWith('material');
    buttons()
      .find((b) => b.textContent?.includes('Use all the scene defaults'))!
      .click();
    expect(model.reset).toHaveBeenCalledWith();
    // a rename follows
    label.set('RA');
    flush();
    expect(dialog.textContent).toContain('Look of RA');
  });
});
