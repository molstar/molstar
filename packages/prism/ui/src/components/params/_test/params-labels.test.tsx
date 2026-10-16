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
import { Color } from '@molstar/core/util/color/color';
import { render } from '../../../dom/render.js';
import { ParamsFormModel } from '../../../forms/params-form.js';
import { ParamsForm } from '../params-form.js';

/** An element's accessible name as a screen reader computes it here: `aria-labelledby`'s texts, else `aria-label`. */
function nameOf(el: Element): string {
  const ids = el.getAttribute('aria-labelledby');
  if (ids) {
    return ids
      .split(/\s+/)
      .map((id) => el.ownerDocument.getElementById(id)?.textContent?.trim() ?? '')
      .join(' ')
      .trim();
  }
  return el.getAttribute('aria-label')?.trim() ?? '';
}

describe('ParamsForm labels (accessibility.md, WCAG 4.1.2)', () => {
  it('names every control by its row label', () => {
    const form = new ParamsFormModel({
      pdbId: PD.Text('1cbs', { label: 'PDB Id(s)', description: 'The entry to load' }),
      count: PD.Numeric(3),
      opacity: PD.Numeric(0.5, { min: 0, max: 1, step: 0.1 }),
      shadows: PD.Boolean(true),
      tint: PD.Color(Color(0xff0000)),
      notes: PD.Text('', { multiline: true }),
      occlusion: PD.MappedStatic('on', {
        on: PD.Group({ samples: PD.Numeric(32) }),
        off: PD.Group({}),
      }),
      offset: PD.Vec3([0, 0, 0] as any),
      range: PD.Interval([0, 1], { min: 0, max: 1 }),
    });
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <ParamsForm form={form} />, container);

    const controls = [
      ...container.querySelectorAll('input, textarea, [role=switch], button[aria-expanded]:not([aria-label^="About"])'),
    ];
    expect(controls.length).toBeGreaterThan(8);
    for (const c of controls) expect(nameOf(c), c.outerHTML).not.toBe('');

    const byName = (name: string) => controls.filter((c) => nameOf(c) === name);
    // The label's text, not its ? (help) button: "PDB Id(s)", not "PDB Id(s) About PDB Id(s)".
    expect(byName('PDB Id(s)').map((c) => c.tagName)).toEqual(['INPUT']);
    // An unbounded number is text (a native number field shows the locale's decimal comma), with a decimal keyboard.
    expect(byName('Count').map((c) => [(c as HTMLInputElement).type, c.getAttribute('inputmode')])).toEqual([
      ['text', 'decimal'],
    ]);
    // A slider: both its range and its number field.
    expect(byName('Opacity').map((c) => (c as HTMLInputElement).type)).toEqual(['range', 'text']);
    expect(byName('Shadows').map((c) => c.getAttribute('role'))).toEqual(['switch']);
    expect(byName('Occlusion').map((c) => c.getAttribute('role'))).toEqual(['switch']);
    expect(byName('Notes').map((c) => c.tagName)).toEqual(['TEXTAREA']);
    // The color button: its label and its value.
    expect(controls.some((c) => nameOf(c) === 'Tint 255 0 0' || nameOf(c).startsWith('Tint '))).toBe(true);
    expect(controls.filter((c) => nameOf(c).startsWith('Offset ')).length).toBe(3);
    expect(controls.filter((c) => nameOf(c).startsWith('Range ')).map(nameOf)).toEqual(['Range from', 'Range to']);
    container.remove();
  });
});
