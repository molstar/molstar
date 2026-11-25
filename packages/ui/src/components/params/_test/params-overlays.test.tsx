/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { afterEach, describe, expect, it } from 'vitest';
import { render } from '../../../dom/index.js';
import { ParamsFormModel } from '../../../forms/params-form.js';
import { ParamsForm } from '../params-form.js';
import { fieldNamed } from './named.js';

const flush = () => new Promise((r) => setTimeout(r, 0));

afterEach(() => {
  document.body.innerHTML = '';
});

function mount(form: ParamsFormModel<any>) {
  const container = document.createElement('div');
  document.body.append(container);
  render(() => <ParamsForm form={form} />, container);
  return container;
}

/** The rows' count of elements, in order (a list or panel that opened inline would add some between them). */
const rowCount = (container: HTMLElement) => container.querySelectorAll('[data-field]').length;

describe('ParamsForm overlays: nothing opens inline', () => {
  it('opens a select, a multi-select and a color list over their field', async () => {
    const form = new ParamsFormModel({
      mode: PD.Select('a', PD.arrayToOptions(['a', 'b'])),
      chains: PD.MultiSelect(['A'], PD.arrayToOptions(['A', 'B', 'C'])),
      colors: PD.ColorList('viridis' as any),
    });
    const container = mount(form);
    const rows = rowCount(container);
    for (const label of ['Mode', 'Chains', 'Colors']) {
      const field = fieldNamed(container, label)!;
      field.click();
      await flush();
      const list = container.querySelector('[role=listbox]')!;
      expect(list.closest('[popover]')!.getAttribute('data-placement'), label).toBe('over');
      expect(rowCount(container)).toBe(rows);
      field.click();
      await flush();
      expect(container.querySelector('[role=listbox]'), label).toBeNull();
    }
  });

  it('keeps a multi-select open while its options toggle', async () => {
    const form = new ParamsFormModel({ chains: PD.MultiSelect(['A'], PD.arrayToOptions(['A', 'B', 'C'])) });
    const container = mount(form);
    fieldNamed(container, 'Chains')!.click();
    await flush();
    const options = () => [...container.querySelectorAll<HTMLButtonElement>('[role=option]')];
    expect(options().map((o) => o.getAttribute('aria-selected'))).toEqual(['true', 'false', 'false']);
    options()[2].click();
    await flush();
    expect((form.values.value as any).chains).toEqual(['A', 'C']);
    expect(options().map((o) => o.getAttribute('aria-selected'))).toEqual(['true', 'false', 'true']);
    expect(fieldNamed(container, 'Chains')!.textContent).toContain('A, C');
  });

  it('opens a color picker and a matrix in a panel under their field', async () => {
    const form = new ParamsFormModel({ tint: PD.Color(0xff0000 as any), transform: PD.Mat4(Mat4.identity()) });
    const container = mount(form);
    const rows = rowCount(container);
    const [color, matrix] = [...container.querySelectorAll<HTMLButtonElement>('button[aria-haspopup=dialog]')];
    color.click();
    await flush();
    expect(container.querySelector('[role=dialog]')!.closest('[popover]')).not.toBeNull();
    expect(container.querySelector('[role=dialog] [aria-label=Palette]')).not.toBeNull();
    expect(rowCount(container)).toBe(rows);
    color.click();
    await flush();
    matrix.click();
    await flush();
    const cells = container.querySelectorAll<HTMLInputElement>('[role=dialog] input');
    expect(cells.length).toBe(16);
    cells[3].value = '5';
    cells[3].dispatchEvent(new Event('change', { bubbles: true }));
    await flush();
    // Row 1, column 4: the translation's x (stored by columns).
    expect((form.values.value as any).transform[12]).toBe(5);
    expect(matrix.textContent).toContain('4 × 4 matrix');
  });
});

describe('ParamsForm levels', () => {
  it('sets a level on each nested block: a group, its group, Advanced, a choice’s params', async () => {
    const form = new ParamsFormModel({
      mode: PD.Select('a', PD.arrayToOptions(['a', 'b']), { isEssential: true }),
      outer: PD.Group(
        { inner: PD.Group({ depth: PD.Numeric(2) }, { isExpanded: true }), after: PD.Boolean(true) },
        { isExpanded: true },
      ),
    });
    const container = mount(form);
    [...container.querySelectorAll<HTMLButtonElement>('button[aria-expanded]')]
      .find((b) => b.textContent === 'Advanced')!
      .click();
    await flush();
    const levelOf = (label: string) => {
      const text = [...container.querySelectorAll('[id^=msui-field]')].find((el) => el.textContent === label)!;
      return (text.closest('[style*="--msui-params-level"]') as HTMLElement | null)?.style.getPropertyValue(
        '--msui-params-level',
      );
    };
    expect(levelOf('Mode')).toBeUndefined();
    // Advanced is level 1; the group's rows level 2; its group's rows level 3.
    expect(levelOf('After')).toBe('2');
    expect(levelOf('Depth')).toBe('3');
  });
});

describe('ParamsForm object list', () => {
  it('adds, names, edits and removes items', async () => {
    const form = new ParamsFormModel({
      channels: PD.ObjectList({ name: PD.Text('c'), level: PD.Numeric(1) }, (e) => `${e.name} at ${e.level}`, {
        defaultValue: [{ name: 'a', level: 1 }],
        ctor: () => ({ name: 'new', level: 5 }),
      }),
    });
    const container = mount(form);
    expect(container.textContent).toContain('1 item');
    container.querySelector<HTMLButtonElement>('button[aria-label="Add to Channels"]')!.click();
    await flush();
    expect(container.textContent).toContain('2 items');
    expect(container.textContent).toContain('new at 5');
    container.querySelectorAll<HTMLButtonElement>('button[aria-label=Edit]')[1].click();
    await flush();
    const name = [...container.querySelectorAll<HTMLInputElement>('input')].find((i) => i.value === 'new')!;
    name.value = 'renamed';
    name.dispatchEvent(new Event('change', { bubbles: true }));
    await flush();
    expect(container.textContent).toContain('renamed at 5');
    container.querySelectorAll<HTMLButtonElement>('button[aria-label=Remove]')[0].click();
    await flush();
    expect((form.values.value as any).channels).toEqual([{ name: 'renamed', level: 5 }]);
  });

  it("an open item's form follows the item at its place: a remove, a load elsewhere; its edits replace that item", async () => {
    const form = new ParamsFormModel({
      channels: PD.ObjectList({ name: PD.Text('c'), level: PD.Numeric(1) }, (e) => `${e.name} at ${e.level}`, {
        defaultValue: [
          { name: 'a', level: 1 },
          { name: 'b', level: 2 },
        ],
        ctor: () => ({ name: 'new', level: 5 }),
      }),
    });
    const container = mount(form);
    container.querySelectorAll<HTMLButtonElement>('button[aria-label=Edit]')[0].click();
    await flush();
    const nameField = () => container.querySelector<HTMLInputElement>('[data-open] input:not([inputmode])')!;
    expect(nameField().value).toBe('a');
    // Removed: the open form is now the item that took its place.
    container.querySelectorAll<HTMLButtonElement>('button[aria-label=Remove]')[0].click();
    await flush();
    expect(nameField().value).toBe('b');
    nameField().value = 'b2';
    nameField().dispatchEvent(new Event('change', { bubbles: true }));
    await flush();
    expect((form.values.value as any).channels).toEqual([{ name: 'b2', level: 2 }]);
    // Loaded elsewhere (a preset, an update): the form shows the new item, and edits it.
    form.load({ channels: [{ name: 'x', level: 9 }] } as any);
    await flush();
    expect(nameField().value).toBe('x');
    nameField().value = 'x2';
    nameField().dispatchEvent(new Event('change', { bubbles: true }));
    await flush();
    expect((form.values.value as any).channels).toEqual([{ name: 'x2', level: 9 }]);
    expect(container.textContent).toContain('x2 at 9');
  });
});
