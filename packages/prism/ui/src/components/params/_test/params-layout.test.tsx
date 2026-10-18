/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { flushFrame } from '../../../reactive/frame-state.js';
import { describe, expect, it } from 'vitest';
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

/** The labels of the rows that show (not under a hidden element). */
function shown(container: HTMLElement) {
  return [...container.querySelectorAll('[data-tip-overflow]')]
    .filter((el) => !el.closest('[hidden]'))
    .map((el) => el.textContent);
}

/** A fold's header (Advanced, a category, a group). */
function fold(container: HTMLElement, title: string) {
  return [...container.querySelectorAll<HTMLButtonElement>('button[aria-expanded]')].find(
    (b) => b.textContent?.trim() === title,
  );
}

/** A representation's params, as Mol*'s: a type with its own params (essential and not), a color theme. */
const cartoon = PD.Group({
  alpha: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }, { label: 'Opacity', isEssential: true }),
  quality: PD.Select('auto', PD.arrayToOptions(['auto', 'custom']), { isEssential: true }),
  detail: PD.Numeric(1, { min: 0, max: 3, step: 1 }, { hideIf: (p) => p.quality !== 'custom', isEssential: true }),
  emissive: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),
  density: PD.Numeric(0.2, { min: 0, max: 1, step: 0.01 }),
  celShaded: PD.Boolean(false, { category: 'Shading' }),
  secret: PD.Text('x', { isHidden: true }),
});
const spacefill = PD.Group({ sizeFactor: PD.Numeric(1, { min: 0, max: 2, step: 0.1 }, { isEssential: true }) });
const params = {
  type: PD.MappedStatic('cartoon', { cartoon, spacefill }),
  colorTheme: PD.MappedStatic('uniform', {
    uniform: PD.Group({ value: PD.Color(0xff0000 as any) }),
    'chain-id': PD.Group({}),
  }),
};

describe('ParamsForm layout (v5 ParameterControls)', () => {
  it("shows a choice's select and a sliders button that opens the chosen case's params", async () => {
    const form = new ParamsFormModel(params);
    const container = mount(form);
    expect(shown(container)).toEqual(['Type', 'Color Theme']);
    const settings = container.querySelector('button[aria-label="Type settings"]') as HTMLButtonElement;
    expect(settings.getAttribute('aria-expanded')).toBe('false');

    settings.click();
    await flush();
    expect(settings.getAttribute('aria-expanded')).toBe('true');
    // Essentials up front (Detail hidden by its hideIf), the rest folded under Advanced.
    expect(shown(container)).toEqual(['Type', 'Opacity', 'Quality', 'Color Theme']);
    expect(fold(container, 'Advanced')).toBeDefined();
    expect(container.textContent).not.toContain('Emissive');
  });

  it('opens Advanced with the rest, a category under its own fold; keeps the folds per form across choices', async () => {
    const form = new ParamsFormModel(params);
    const container = mount(form);
    (container.querySelector('button[aria-label="Type settings"]') as HTMLButtonElement).click();
    await flush();
    fold(container, 'Advanced')!.click();
    await flush();
    expect(shown(container)).toEqual(['Type', 'Opacity', 'Quality', 'Emissive', 'Density', 'Color Theme']);
    fold(container, 'Shading')!.click();
    await flush();
    expect(shown(container)).toContain('Cel Shaded');

    // Another type and back: the button stays open, Advanced and Shading as they were left.
    form.set('type', { name: 'spacefill', params: { sizeFactor: 1 } } as any);
    await flush();
    expect(shown(container)).toEqual(['Type', 'Size Factor', 'Color Theme']);
    expect(fold(container, 'Advanced')).toBeUndefined();
    form.set('type', params.type.defaultValue as any);
    await flush();
    expect(shown(container)).toContain('Cel Shaded');
    expect(shown(container)).toContain('Emissive');
  });

  it("follows a case's hideIf on its values; produces the same values as before", async () => {
    const form = new ParamsFormModel(params);
    const container = mount(form);
    (container.querySelector('button[aria-label="Type settings"]') as HTMLButtonElement).click();
    await flush();
    expect(shown(container)).not.toContain('Detail');
    const type = (form.values.value as any).type;
    form.set('type', { name: 'cartoon', params: { ...type.params, quality: 'custom' } } as any);
    await flush();
    expect(shown(container)).toEqual(['Type', 'Opacity', 'Quality', 'Detail', 'Color Theme']);

    // A slider in the case's params sets the mapped value, keeping its name and the other params.
    const range = container.querySelector('input[type=range]') as HTMLInputElement;
    range.value = '0.5';
    range.dispatchEvent(new Event('input', { bubbles: true }));
    // A drag's value is set in the frame (the frame lane, rendering.md §2.3).
    flushFrame();
    await flush();
    expect((form.values.value as any).type).toEqual({
      name: 'cartoon',
      params: { ...type.params, quality: 'custom', alpha: 0.5 },
    });
  });

  it('has no button for a case without params; a flat case shows its rows at the same level', async () => {
    const form = new ParamsFormModel({
      ...params,
      size: PD.MappedStatic('physical', {
        physical: PD.Group({ scale: PD.Numeric(1, { min: 0, max: 2, step: 0.1 }) }, { isFlat: true }),
        uniform: PD.Group({}),
      }),
    });
    const container = mount(form);
    expect(container.querySelector('button[aria-label="Color Theme settings"]')).not.toBeNull();
    expect(container.querySelector('button[aria-label="Size settings"]')).toBeNull();
    expect(shown(container)).toEqual(['Type', 'Color Theme', 'Size', 'Scale']);
    form.set('colorTheme', { name: 'chain-id', params: {} } as any);
    await flush();
    expect(container.querySelector('button[aria-label="Color Theme settings"]')).toBeNull();
  });

  it('splits a top-level form only when it has both essentials and the rest; hides a fold with nothing to show', async () => {
    const plain = mount(new ParamsFormModel({ a: PD.Boolean(true), b: PD.Text('') }));
    expect(fold(plain, 'Advanced')).toBeUndefined();
    expect(shown(plain)).toEqual(['A', 'B']);

    const form = new ParamsFormModel({
      mode: PD.Select('a', PD.arrayToOptions(['a', 'b']), { isEssential: true }),
      onlyB: PD.Numeric(1, {}, { hideIf: (p) => p.mode !== 'b' }),
    });
    const container = mount(form);
    expect(fold(container, 'Advanced')!.closest('[hidden]')).not.toBeNull();
    form.set('mode', 'b');
    await flush();
    expect(fold(container, 'Advanced')!.closest('[hidden]')).toBeNull();
  });

  it("opens a group as its isExpanded says and keeps a flat group's rows at its level", () => {
    const container = mount(
      new ParamsFormModel({
        open: PD.Group({ inside: PD.Boolean(true) }, { isExpanded: true }),
        closed: PD.Group({ hidden: PD.Boolean(true) }),
        flat: PD.Group({ level: PD.Boolean(true) }, { isFlat: true }),
      }),
    );
    expect(shown(container)).toEqual(['Inside', 'Level']);
    expect(fold(container, 'Open')!.getAttribute('aria-expanded')).toBe('true');
    expect(fold(container, 'Closed')!.getAttribute('aria-expanded')).toBe('false');
  });
});
