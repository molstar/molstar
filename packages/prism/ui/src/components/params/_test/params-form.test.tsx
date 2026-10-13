/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { render } from '../../../dom/render.js';
import { ParamsFormModel } from '../../../forms/params-form.js';
import { observableState } from '../../../reactive/observable-state.js';
import { ParamsForm, paramsMatch } from '../params-form.js';
import { presetIndex } from '../presets.js';

const flush = () => new Promise((r) => setTimeout(r, 0));

function mount(form: ParamsFormModel<any>) {
  const container = document.createElement('div');
  document.body.append(container);
  render(() => <ParamsForm form={form} />, container);
  return container;
}

describe('ParamsForm', () => {
  it("applies a group's preset from its header, names the matching one, and a flat group's from its row", async () => {
    const material = PD.Group(
      { metalness: PD.Numeric(0), roughness: PD.Numeric(1) },
      {
        presets: [
          [{ metalness: 0, roughness: 1 }, 'Matte'],
          [{ metalness: 1, roughness: 0.6 }, 'Metallic'],
        ],
      },
    );
    const size = PD.Group(
      { width: PD.Numeric(1920), height: PD.Numeric(1080) },
      {
        isFlat: true,
        presets: [
          [{ width: 1920, height: 1080 }, 'Full HD'],
          [{ width: 3840, height: 2160 }, '4K'],
        ],
      },
    );
    const form = new ParamsFormModel({ material, size });
    const container = mount(form);
    const button = (label: string) => container.querySelector(`button[aria-label="${label}"]`) as HTMLButtonElement;
    const option = (label: string) =>
      [...container.querySelectorAll('[role=listbox] button')].find((b) => b.textContent?.trim() === label) as
        | HTMLButtonElement
        | undefined;

    expect(container.textContent).toContain('Matte');
    button('Material presets').click();
    await flush();
    option('Metallic')!.click();
    await flush();
    expect((form.values.value as any).material).toEqual({ metalness: 1, roughness: 0.6 });
    expect(container.textContent).toContain('Metallic');
    expect(container.querySelector('[role=listbox]')).toBeNull();

    expect(button('Size presets').textContent).toContain('Full HD');
    button('Size presets').click();
    await flush();
    option('4K')!.click();
    await flush();
    expect((form.values.value as any).size).toEqual({ width: 3840, height: 2160 });
    expect(button('Size presets').textContent).toContain('4K');
  });

  it('finds the preset a value matches: all of it, or the fields a group preset sets', () => {
    const presets: PD.SelectOption<any>[] = [
      [{ a: 1 }, 'One'],
      [{ a: 2, b: 3 }, 'Two'],
    ];
    expect(presetIndex(presets, { a: 1, c: 9 }, true)).toBe(0);
    expect(presetIndex(presets, { a: 1, c: 9 })).toBeUndefined();
    expect(presetIndex(presets, { a: 2, b: 3 })).toBe(1);
    expect(presetIndex(presets, { a: 5 }, true)).toBeUndefined();
  });

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

  it('filters rows by a search: matches show (groups open), a group whose name matches shows all its rows', async () => {
    const params = {
      pickPadding: PD.Numeric(3, { min: 0, max: 10, step: 1 }),
      trackball: PD.Group({ rotateSpeed: PD.Numeric(5), zoomSpeed: PD.Numeric(7) }),
      tooltips: PD.Boolean(true),
    };
    const form = new ParamsFormModel(params);
    const filter = observableState('');
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <ParamsForm form={form} filter={filter.asReadonly} />, container);
    // The visible rows' labels (a group's header is not a row).
    const shown = () =>
      [...container.querySelectorAll('[data-tip-overflow]')]
        .filter((el) => !el.closest('[hidden]'))
        .map((el) => el.textContent);

    filter.set('zoom');
    await flush();
    expect(shown()).toEqual(['Zoom Speed']);

    filter.set('trackball');
    await flush();
    expect(shown()).toEqual(['Rotate Speed', 'Zoom Speed']);

    filter.set('');
    await flush();
    expect(shown()).toEqual(['Pick Padding', 'Tooltips']);

    expect(paramsMatch(params, 'zoom')).toBe(true);
    expect(paramsMatch(params, 'nothing')).toBe(false);
  });

  it('marks a label with help (dotted underline); its text and the ? after it open the help; a long label keeps its ?', async () => {
    const form = new ParamsFormModel({
      backgroundOcclusionRadiusScale: PD.Numeric(
        0.5,
        { min: 0, max: 2 },
        { description: 'How far occlusion samples.' },
      ),
      plain: PD.Boolean(true),
    });
    const container = mount(form);
    await flush();
    const [helped, plain] = [...container.querySelectorAll<HTMLElement>('[data-tip-overflow]')];
    expect(helped.hasAttribute('data-help')).toBe(true);
    expect(plain.hasAttribute('data-help')).toBe(false);
    // The ? is the text's sibling, outside the part that is cut off.
    const toggle = helped.nextElementSibling as HTMLButtonElement;
    expect(toggle.getAttribute('aria-label')).toBe('About Background Occlusion Radius Scale');
    const help = () => container.textContent?.includes('How far occlusion samples.');
    expect(help()).toBe(false);
    helped.click();
    await flush();
    expect(help()).toBe(true);
    expect(toggle.getAttribute('aria-expanded')).toBe('true');
    toggle.click();
    await flush();
    expect(help()).toBe(false);
  });
});
