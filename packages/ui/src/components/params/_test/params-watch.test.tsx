/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { afterEach, beforeEach, describe, expect, it } from 'vitest';
import { DomProbe, flush, render } from '../../../dom/index.js';
import { ParamsFormModel } from '../../../forms/params-form.js';
import type { State } from '@molstar/core/reactive';
import { state } from '../../../state.js';
import { ParamsForm } from '../params-form.js';
import { fieldNamed } from './named.js';

// The form's rows are views of their field (dom §2.2): the DOM made once, one `watch` a row writing what
// differs, a row listening to what decides its output.

// Set before rendering (`DomProbe`): told every `watch` render.
let rendered = 0;
beforeEach(() => {
  rendered = 0;
  DomProbe.on = (kind) => {
    if (kind === 'render') rendered++;
  };
});
afterEach(() => {
  DomProbe.on = undefined;
});

function mount(form: ParamsFormModel<any>, filter?: State<string>) {
  const container = document.createElement('div');
  document.body.append(container);
  render(() => <ParamsForm form={form} {...(filter ? { filter: filter } : {})} />, container);
  return container;
}

/** The `watch` renders a step causes (the first renders at the mount are not counted). */
function renders(step: () => void): number {
  flush();
  const before = rendered;
  step();
  flush();
  return rendered - before;
}

const params = {
  size: PD.Numeric(1, { min: 0, max: 2, step: 0.1 }, { label: 'Size' }),
  mode: PD.Select('a', PD.arrayToOptions(['a', 'b']), { label: 'Mode' }),
  advanced: PD.Boolean(false, { label: 'Advanced mode' }),
  detail: PD.Numeric(1, { min: 0, max: 3, step: 1 }, { label: 'Detail', hideIf: (p) => !p.advanced }),
  name: PD.Text('x', { label: 'Name' }),
};

describe('the params form rows as views of their field', () => {
  it('an edit of one field renders that field’s row only; a select’s label follows it', () => {
    const form = new ParamsFormModel(params);
    const container = mount(form);
    const label = () => fieldNamed(container, 'Mode')!.querySelector('span')!.textContent;
    expect(label()).toBe('A');

    // The slider's row renders; the others (and the hideIf row, whose condition did not change) do not.
    expect(renders(() => form.set('size', 1.5))).toBe(1);
    expect(renders(() => form.set('name', 'y'))).toBe(1);

    expect(renders(() => form.set('mode', 'b'))).toBe(1);
    expect(label()).toBe('B');
    // A burst in one task is one render a row.
    expect(
      renders(() => {
        form.set('mode', 'a');
        form.set('mode', 'b');
        form.set('mode', 'a');
      }),
    ).toBe(1);
    expect(label()).toBe('A');
  });

  it('a hideIf row hides and shows with its condition only (its row’s watch follows the gate)', () => {
    const form = new ParamsFormModel(params);
    const container = mount(form);
    const detail = () => [...container.querySelectorAll('[hidden]')].some((el) => el.textContent?.includes('Detail'));
    expect(detail()).toBe(true);

    // The toggle's own row and the hideIf row render; the others do not.
    expect(renders(() => form.set('advanced', true))).toBe(2);
    expect(detail()).toBe(false);
    expect(renders(() => form.set('advanced', false))).toBe(2);
    expect(detail()).toBe(true);
  });

  it('a search hides rows that do not match, and a section with none showing', () => {
    const grouped = {
      first: PD.Numeric(1, { min: 0, max: 2, step: 0.1 }, { label: 'Alpha', category: 'One' }),
      second: PD.Text('x', { label: 'Beta', category: 'Two' }),
    };
    const form = new ParamsFormModel(grouped);
    const filter = state('');
    const container = mount(form, filter);
    const hidden = () => [...container.querySelectorAll('[hidden]')].map((el) => el.textContent?.trim());
    expect(hidden()).toEqual([]);

    filter.set('bet');
    flush();
    // Alpha's row, and its category with it.
    expect(hidden().some((t) => t?.includes('Alpha'))).toBe(true);
    expect(hidden().some((t) => t?.includes('Beta'))).toBe(false);

    filter.set('');
    flush();
    expect(hidden()).toEqual([]);
  });

  it('writes a vector, a matrix and an interval from the value, one watch each', () => {
    const forms = {
      v: PD.Vec3([1, 2, 3] as any, {}, { label: 'Vector' }),
      m: PD.Mat4([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1] as any, { label: 'Matrix' }),
      i: PD.Interval([0.2, 0.6], { min: 0, max: 1, step: 0.1 }, { label: 'Range' }),
    };
    const form = new ParamsFormModel(forms);
    const container = mount(form);
    const vec = () =>
      [...container.querySelectorAll<HTMLInputElement>('input[aria-label^="Vector"]')].map((i) => i.value);
    expect(vec()).toEqual(['1', '2', '3']);
    expect(container.textContent).toContain('0.2–0.6');
    expect(container.textContent).toContain('Identity');

    // Three numbers, one render; the matrix's summary follows its value.
    const m = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 5, 0, 0, 1];
    expect(
      renders(() => {
        form.set('v', [4, 5, 6] as any);
        form.set('m', m as any);
        form.set('i', [0.1, 0.9]);
      }),
    ).toBe(3);
    expect(vec()).toEqual(['4', '5', '6']);
    expect(container.textContent).toContain('4 × 4 matrix');
    expect(container.textContent).toContain('0.1–0.9');
  });
});
