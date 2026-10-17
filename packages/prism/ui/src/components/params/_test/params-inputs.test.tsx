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
import { afterEach, describe, expect, it } from 'vitest';
import { render } from '../../../dom/render.js';
import { ParamsFormModel } from '../../../forms/params-form.js';
import { ParamsForm } from '../params-form.js';

afterEach(() => {
  document.body.innerHTML = '';
});

const flush = () => new Promise((r) => setTimeout(r, 0));

function mount(form: ParamsFormModel<any>) {
  const container = document.createElement('div');
  document.body.append(container);
  render(() => <ParamsForm form={form} />, container);
  return container;
}

/** Types into a field and commits it, as the person does (change fires on blur and Enter). */
function type(input: HTMLInputElement, text: string) {
  input.value = text;
  input.dispatchEvent(new Event('change', { bubbles: true }));
}

describe('typed numbers are rewritten to the value that holds', () => {
  it('a slider’s box: an out-of-range value that clamps to the current one, and text that is no number', async () => {
    const form = new ParamsFormModel({ level: PD.Numeric(100, { min: 0, max: 100, step: 1 }) });
    const box = mount(form).querySelector<HTMLInputElement>('input:not([type=range])')!;
    await flush();
    expect(box.value).toBe('100');
    type(box, '500');
    await flush();
    expect(box.value).toBe('100');
    expect((form.values.value as any).level).toBe(100);
    type(box, 'abc');
    await flush();
    expect(box.value).toBe('100');
    type(box, '40');
    await flush();
    expect(box.value).toBe('40');
    expect((form.values.value as any).level).toBe(40);
  });

  it('a number field with no range: "1.50" for a value that is 1.5', async () => {
    const form = new ParamsFormModel({ size: PD.Numeric(1.5) });
    const box = mount(form).querySelector<HTMLInputElement>('input[type=text]')!;
    await flush();
    type(box, '1.50');
    await flush();
    expect(box.value).toBe('1.5');
  });

  it('a vector’s fields: "1.50" is shown as 1.5, and text that is no number goes back to the value', async () => {
    const form = new ParamsFormModel({ offset: PD.Vec3([1.5, 0, 0] as any) });
    const container = mount(form);
    await flush();
    const x = container.querySelector<HTMLInputElement>('input[aria-label$="x"], input[aria-label="x"]')!;
    expect(x.value).toBe('1.5');
    type(x, '1.50');
    await flush();
    expect(x.value).toBe('1.5');
    type(x, 'zzz');
    await flush();
    expect(x.value).toBe('1.5');
    type(x, '2.0');
    await flush();
    expect(x.value).toBe('2');
    expect((form.values.value as any).offset[0]).toBe(2);
  });
});
