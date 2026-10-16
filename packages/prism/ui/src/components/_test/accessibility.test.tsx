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
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { render } from '../../dom/render.js';
import { ParamsFormModel } from '../../forms/params-form.js';
import { observableState } from '../../reactive/observable-state.js';
import { expectNoAxeViolations } from '../../testing/axe.js';
import { Choices } from '../choices/choices.js';
import { ColorPicker } from '../color/color-picker.js';
import { openContextMenu } from '../context-menu/context-menu.js';
import { ParamsForm } from '../params/params-form.js';
import { Select } from '../select/select.js';

// Automated checks (axe-core, jsdom) of the shared components (accessibility.md, Automated checks): a params form with
// one control of each kind, the open select list, the color picker, a choice row and a context menu.

const flush = () => new Promise((r) => setTimeout(r, 0));

const onOff = (on: boolean, params: PD.Params = {}, info: PD.Info = {}) =>
  PD.MappedStatic(on ? 'on' : 'off', { on: PD.Group(params), off: PD.Group({}) }, info);

function mount(view: () => any) {
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(view, container);
  return { container, dispose };
}

describe('components have no accessibility violations (axe)', () => {
  beforeEach(() => {
    vi.stubGlobal(
      'ResizeObserver',
      class {
        observe() {}
        disconnect() {}
        unobserve() {}
      },
    );
  });
  afterEach(() => {
    document.body.innerHTML = '';
    vi.unstubAllGlobals();
  });

  it('a params form with one control of each kind', async () => {
    const form = new ParamsFormModel({
      boolean: PD.Boolean(true),
      number: PD.Numeric(5, { min: 0, max: 20, step: 0.5 }),
      plain: PD.Numeric(1920, undefined, { label: 'Width (px)' }),
      text: PD.Text('1tqn', { placeholder: 'PDB ID' }),
      notes: PD.Text('a\nb', { multiline: true, label: 'Notes' }),
      select: PD.Select('cartoon', [
        ['cartoon', 'Cartoon'],
        ['spacefill', 'Spacefill'],
      ]),
      color: PD.Color(0x3362b2 as any),
      file: PD.File({ accept: '.cif,.pdb' }),
      files: PD.FileList({ multiple: true }),
      outline: onOff(true, { scale: PD.Numeric(1, { min: 1, max: 5, step: 1 }) }, { label: 'Outline' }),
      mode: PD.MappedStatic(
        'auto',
        {
          auto: PD.Group({}, { isFlat: true }),
          custom: PD.Group({ width: PD.Numeric(1920), height: PD.Numeric(1080) }, { isFlat: true }),
        },
        {
          options: [
            ['auto', 'Auto'],
            ['custom', 'Custom'],
          ],
        },
      ),
      group: PD.Group({ samples: PD.Numeric(32, { min: 1, max: 256, step: 1 }) }, { label: 'Group', isExpanded: true }),
      folded: PD.Group({ bias: PD.Numeric(0.8) }, { label: 'Folded' }),
      size: PD.Group(
        { width: PD.Numeric(1920), height: PD.Numeric(1080) },
        {
          isFlat: true,
          presets: [
            [{ width: 1920, height: 1080 }, 'Full HD'],
            [{ width: 3840, height: 2160 }, '4K'],
          ],
        },
      ),
      channels: PD.ObjectList(
        { name: PD.Text('channel'), level: PD.Numeric(1.5, { min: 0, max: 5, step: 0.1 }) },
        (e) => `${e.name} at ${e.level}σ`,
        { defaultValue: [{ name: '2fo-fc', level: 1.5 }] },
      ),
    });
    const t = mount(() => <ParamsForm form={form} />);
    await flush();
    expect(t.container.querySelectorAll('input').length).toBeGreaterThan(5);
    await expectNoAxeViolations(t.container);
    // Its select and a presets list open (popovers: in the top layer, named by their buttons).
    const openers = t.container.querySelectorAll<HTMLElement>('button[aria-haspopup=listbox]');
    expect(openers.length).toBeGreaterThan(1);
    for (const opener of openers) {
      opener.click();
      await flush();
      expect(document.querySelectorAll('[role=listbox]').length).toBe(1);
      await expectNoAxeViolations(document.body);
      opener.click();
      await flush();
    }
    t.dispose();
  });

  it('a select with its list open', async () => {
    const value = observableState('b');
    const t = mount(() => (
      <Select
        aria-label="Letter"
        options={[
          ['a', 'Alpha'],
          ['b', 'Beta'],
        ]}
        value={value}
        onChange={(v) => value.set(v)}
      />
    ));
    t.container.querySelector<HTMLElement>('button[aria-haspopup=listbox]')!.click();
    await flush();
    expect(t.container.querySelector('[role=listbox]')).not.toBeNull();
    await expectNoAxeViolations(t.container);
    t.dispose();
  });

  it('the color picker', async () => {
    const color = observableState(0x3362b2);
    const t = mount(() => <ColorPicker value={color.asReadonly} onChange={(c) => color.set(c)} />);
    await flush();
    await expectNoAxeViolations(t.container);
    t.dispose();
  });

  it('a row of choices', async () => {
    const t = mount(() => (
      <Choices
        label="Style"
        choices={[
          { id: 'a', label: 'Cartoon' },
          { id: 'b', label: 'Surface' },
        ]}
        onPick={() => {}}
      />
    ));
    await expectNoAxeViolations(t.container);
    t.dispose();
  });

  it('a context menu', async () => {
    const root = document.createElement('div');
    document.body.append(root);
    openContextMenu(root, { x: 10, y: 20 }, [
      { id: 'a', label: 'Color by pLDDT', icon: 'gauge', shortcut: '⌘P', run: () => {} },
      { id: 'b', label: 'Hide', disabled: 'Nothing to hide', run: () => {} },
    ]);
    expect(root.querySelector('[role=menu]')).not.toBeNull();
    await expectNoAxeViolations(root);
  });
});
