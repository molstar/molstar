/** @jest-environment jsdom */
/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { render } from '../../../dom/render.js';
import { ParamsFormModel } from '../../../forms/params-form.js';
import { ParamsForm } from '../params-form.js';

const flush = () => new Promise((r) => setTimeout(r, 0));

function mount(form: ParamsFormModel<any>) {
  const container = document.createElement('div');
  document.body.append(container);
  render(() => <ParamsForm form={form} />, container);
  return container;
}

describe('ParamsForm', () => {
  it('takes dropped files into a file list, names them, clears', async () => {
    const form = new ParamsFormModel({ files: PD.FileList({ multiple: true }) });
    const container = mount(form);
    const pick = container.querySelector('button') as HTMLButtonElement;
    expect(pick.textContent).toContain('Select files…');

    const files = [new File(['a'], '1tqn.cif'), new File(['b'], 'emd_2984.map')];
    const drop = new Event('drop', { bubbles: true, cancelable: true });
    Object.defineProperty(drop, 'dataTransfer', { value: { files } });
    pick.dispatchEvent(drop);
    await flush();
    expect((form.values.value as any).files.map((f: any) => f.name)).toEqual(['1tqn.cif', 'emd_2984.map']);
    expect(pick.textContent).toContain('2 files · 1tqn.cif, emd_2984.map');

    (container.querySelector('button[aria-label=Clear]') as HTMLButtonElement).click();
    await flush();
    expect((form.values.value as any).files).toBeNull();
    expect(pick.textContent).toContain('Select files…');
  });

  it('shows an on/off mapped param as a switch; on, its settings open behind the sliders button', async () => {
    const form = new ParamsFormModel({
      outline: PD.MappedStatic('off', {
        on: PD.Group({ scale: PD.Numeric(1, { min: 0, max: 5, step: 1 }) }),
        off: PD.Group({}),
      }),
    });
    const container = mount(form);
    const toggle = container.querySelector('[role=switch]') as HTMLButtonElement;
    expect(toggle.getAttribute('aria-checked')).toBe('false');
    expect(container.querySelector('[aria-label="Outline settings"]')).toBeNull();

    toggle.click();
    await flush();
    expect((form.values.value as any).outline.name).toBe('on');
    expect(toggle.getAttribute('aria-checked')).toBe('true');
    const settings = container.querySelector('[aria-label="Outline settings"]') as HTMLButtonElement;
    expect(settings).not.toBeNull();
    expect(container.querySelector('input[type=range]')).toBeNull();
    settings.click();
    await flush();
    expect(container.querySelector('input[type=range]')).not.toBeNull();
  });
});
