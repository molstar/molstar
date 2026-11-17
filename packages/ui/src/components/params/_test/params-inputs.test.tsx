/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { afterEach, describe, expect, it } from 'vitest';
import { flushFrame, render } from '../../../dom/index.js';
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

/** Moves a range input, as a drag does (one input event). */
function slide(input: HTMLInputElement, value: number) {
  input.value = String(value);
  input.dispatchEvent(new Event('input', { bubbles: true }));
}

const values = (form: ParamsFormModel<any>) => form.values.value as any;

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

describe('typed numbers keep to the bounds there are', () => {
  it('a number with only a min, or only a max, clamps to it', async () => {
    const form = new ParamsFormModel({ low: PD.Numeric(5, { min: 0 }), high: PD.Numeric(5, { max: 10 }) });
    const [low, high] = mount(form).querySelectorAll<HTMLInputElement>('input[type=text]');
    await flush();
    type(low, '-3');
    type(high, '20');
    await flush();
    expect(values(form)).toMatchObject({ low: 0, high: 10 });
    expect([low.value, high.value]).toEqual(['0', '10']);
  });
});

describe("a number's decimals are its step's", () => {
  it('shows as many decimals as the step has (0.25: two, 2.5: one)', async () => {
    const form = new ParamsFormModel({ a: PD.Numeric(1, { step: 0.25 }), b: PD.Numeric(2.5, { step: 2.5 }) });
    const [a, b] = mount(form).querySelectorAll<HTMLInputElement>('input[type=text]');
    await flush();
    expect([a.value, b.value]).toEqual(['1.00', '2.5']);
  });

  it("keeps a typed value to the step's decimals (shown as held); a value set elsewhere shows as it is", async () => {
    const form = new ParamsFormModel({ a: PD.Numeric(0, { min: 0, max: 1, step: 0.1 }) });
    const box = mount(form).querySelector<HTMLInputElement>('input[type=text]')!;
    type(box, '0.123');
    await flush();
    expect(values(form).a).toBe(0.1);
    expect(box.value).toBe('0.1');
    form.set('a', 0.123);
    await flush();
    expect(box.value).toBe('0.123');
  });
});

describe('a slider in the frame lane', () => {
  it('a drag back to the last dragged value after a change elsewhere is set again', async () => {
    const form = new ParamsFormModel({ level: PD.Numeric(50, { min: 0, max: 100, step: 1 }) });
    const range = mount(form).querySelector<HTMLInputElement>('input[type=range]')!;
    slide(range, 70);
    flushFrame();
    expect(values(form).level).toBe(70);
    form.set('level', 30);
    await flush();
    slide(range, 70);
    flushFrame();
    expect(values(form).level).toBe(70);
  });
});

describe('an interval', () => {
  it("a knob's drag is set once a frame (the frame lane), the ends not crossing", async () => {
    const form = new ParamsFormModel({ i: PD.Interval([0.2, 0.6], { min: 0, max: 1, step: 0.1 }) });
    const [from, to] = mount(form).querySelectorAll<HTMLInputElement>('input[type=range]');
    slide(from, 0.3);
    slide(from, 0.4);
    expect(values(form).i).toEqual([0.2, 0.6]);
    flushFrame();
    expect(values(form).i).toEqual([0.4, 0.6]);
    slide(to, 0.1);
    flushFrame();
    expect(values(form).i).toEqual([0.4, 0.4]);
    await flush();
  });

  it('without bounds: two number fields, not clamped to 0–1, each up to the other end', async () => {
    const form = new ParamsFormModel({ i: PD.Interval([5, 20]) });
    const container = mount(form);
    expect(container.querySelector('input[type=range]')).toBeNull();
    const [from, to] = container.querySelectorAll<HTMLInputElement>('input[type=text]');
    await flush();
    expect([from.value, to.value]).toEqual(['5', '20']);
    type(from, '-100');
    type(to, '1,5e3');
    await flush();
    expect(values(form).i).toEqual([-100, 1500]);
    type(from, '2000');
    await flush();
    expect(values(form).i).toEqual([1500, 1500]);
    expect(from.value).toBe('1500');
  });
});

describe('vector and matrix cells', () => {
  it('a vector reads a comma, goes back on an empty field, keeps to its range and step', async () => {
    const form = new ParamsFormModel({ v: PD.Vec3([0.5, 0, 0] as any, { min: -1, max: 1, step: 0.5 }) });
    const [x, y] = mount(form).querySelectorAll<HTMLInputElement>('input');
    await flush();
    expect(x.value).toBe('0.5');
    type(x, '');
    await flush();
    expect(values(form).v[0]).toBe(0.5);
    expect(x.value).toBe('0.5');
    type(x, '-0,5');
    type(y, '7');
    await flush();
    expect(values(form).v).toEqual([-0.5, 1, 0]);
    expect(y.value).toBe('1.0');
  });

  it('a matrix cell reads a comma and goes back on an empty field', async () => {
    const form = new ParamsFormModel({ m: PD.Mat4([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1] as any) });
    const container = mount(form);
    container.querySelector<HTMLButtonElement>('button[aria-haspopup=dialog]')!.click();
    await flush();
    const cells = container.querySelectorAll<HTMLInputElement>('[role=dialog] input');
    type(cells[0], '');
    await flush();
    expect(values(form).m[0]).toBe(1);
    expect(cells[0].value).toBe('1');
    type(cells[1], '0,5');
    await flush();
    // Row 1, column 2: stored by columns.
    expect(values(form).m[4]).toBe(0.5);
  });
});
