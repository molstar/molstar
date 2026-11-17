/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { afterEach, describe, expect, it } from 'vitest';
import { flush, render } from '../../dom/index.js';
import { axeViolations } from '../../testing/axe.js';
import { ParamsFormModel } from '../params-form.js';
import { matchingGroups, type SettingsGroup, SettingsView } from '../settings-view.js';

// Settings: groups of forms with one search over all (matching groups open, the rest hide, "No setting matches"), on
// `Expandable`, the first open.

const tick = () => new Promise((r) => setTimeout(r, 0));

let dispose: (() => void) | undefined;
afterEach(() => {
  dispose?.();
  dispose = undefined;
  document.body.innerHTML = '';
});

function groups(): SettingsGroup[] {
  const group = (id: string, title: string, params: PD.Params): SettingsGroup => ({
    id,
    title,
    icon: 'settings',
    form: new ParamsFormModel(params),
  });
  return [
    group('look', 'Appearance', {
      background: PD.Color(0xffffff as any, { label: 'Background' }),
      outline: PD.Boolean(false, { label: 'Outline' }),
    }),
    group('camera', 'Camera', {
      fov: PD.Numeric(45, { min: 10, max: 90 }, { label: 'Field of view' }),
      mode: PD.Select('perspective', [
        ['perspective', 'Perspective'],
        ['orthographic', 'Orthographic'],
      ]),
    }),
    group('perf', 'Performance', { quality: PD.Select('auto', [['auto', 'Auto']], { label: 'Quality' }) }),
  ];
}

function mount(query = new State('')) {
  const container = document.createElement('div');
  document.body.append(container);
  const gs = groups();
  dispose = render(() => <SettingsView groups={gs} query={query} />, container);
  const group = (id: string) => container.querySelector<HTMLElement>(`[data-settings-group="${id}"]`)!;
  const expanded = (id: string) => group(id).querySelector('button[aria-expanded]')!.getAttribute('aria-expanded');
  return { container, query, gs, group, expanded };
}

describe('the settings view', () => {
  it('lists the groups as panels, the first open', () => {
    const v = mount();
    expect(v.expanded('look')).toBe('true');
    expect(v.expanded('camera')).toBe('false');
    expect(v.container.textContent).toContain('Background');
    expect(v.container.querySelector('input[aria-label="Search settings"]')).not.toBeNull();
  });

  it('searches every group: matching ones open with only their matching rows, the rest hide', async () => {
    const v = mount();
    const input = v.container.querySelector<HTMLInputElement>('input[aria-label="Search settings"]')!;
    input.value = 'field';
    input.dispatchEvent(new Event('input', { bubbles: true }));
    flush();
    await tick();
    expect(v.group('look').hidden).toBe(true);
    expect(v.group('perf').hidden).toBe(true);
    expect(v.group('camera').hidden).toBe(false);
    expect(v.expanded('camera')).toBe('true');
    expect(v.group('camera').textContent).toContain('Field of view');
    v.query.set('');
    flush();
    await tick();
    // The search let go: back as the person left them.
    expect(v.group('look').hidden).toBe(false);
    expect(v.expanded('camera')).toBe('false');
  });

  it('a group whose title matches shows whole: open, every setting in it', async () => {
    const v = mount(new State('cam'));
    flush();
    await tick();
    expect(v.group('camera').hidden).toBe(false);
    expect(v.group('look').hidden).toBe(true);
    expect(v.expanded('camera')).toBe('true');
    // "cam" is in no setting's label: the whole group shows anyway.
    expect(v.group('camera').textContent).toContain('Field of view');
    expect(v.group('camera').textContent).toContain('Perspective');
    expect(v.container.textContent).not.toContain('No setting matches');
    expect(matchingGroups(v.gs, 'perf').map((g) => g.id)).toEqual(['perf']);
  });

  it('words match across levels: a group title and a setting, a setting and its option', async () => {
    const v = mount(new State('perf qua'));
    flush();
    await tick();
    expect(v.group('perf').hidden).toBe(false);
    expect(v.group('perf').textContent).toContain('Quality');
    expect(matchingGroups(v.gs, 'perf qua').map((g) => g.id)).toEqual(['perf']);
    expect(matchingGroups(v.gs, 'camera field').map((g) => g.id)).toEqual(['camera']);
    // Every word must match somewhere on the path.
    expect(matchingGroups(v.gs, 'perf field')).toEqual([]);
    expect(matchingGroups(v.gs, '  QUALITY   perf ')).toHaveLength(1);
  });

  it('says when nothing matches', () => {
    const v = mount(new State('zzz'));
    expect(v.container.textContent).toContain('No setting matches “zzz”.');
    expect(v.gs.every((g) => v.group(g.id).hidden)).toBe(true);
    expect(matchingGroups(v.gs, 'quality').map((g) => g.id)).toEqual(['perf']);
    expect(matchingGroups(v.gs, ' ')).toHaveLength(3);
  });

  it('is accessible (axe)', async () => {
    const v = mount();
    expect(await axeViolations(v.container)).toBe('');
  });
});
