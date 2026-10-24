/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { flush } from '@molstar/prism-ui/dom/scope';
import { type MockScenePreset, mockScene } from '../../../testing/scene.js';
import { render } from '@molstar/prism-ui/dom/render';
import { ShellModel } from '@molstar/prism-ui/shell/shell';
import { StructureTab } from '../structure-tab.js';
import { domView } from '@molstar/prism-ui/dom/view';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { SceneRow } from '../scene-row.js';
import { SceneRowSlot } from '../scene-tree.js';
import type { SceneInspection, SceneModel, SceneNode } from '../../../model/scene.js';
import { shownInOrder } from '@molstar/prism-ui/testing/virtual-list';

function mount(preset?: MockScenePreset, withScene?: (scene: SceneModel) => SceneModel) {
  const base = mockScene(preset ? { preset } : {});
  const scene = withScene ? withScene(base) : base;
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(() => <StructureTab shell={new ShellModel()} scene={scene} />, container);
  const rows = () => shownInOrder(container, '[role=treeitem] > div:first-child');
  const row = (text: string) => rows().find((r) => r.querySelector('span')?.textContent === text)!;
  return { scene, container, dispose, rows, row };
}

afterEach(() => {
  document.body.innerHTML = '';
});

describe('StructureTab', () => {
  it("opens a row's menu as a popover from its ⋯ button; Escape or a press outside closes it", async () => {
    const { container, row } = mount();
    const more = row('Water').querySelector<HTMLButtonElement>('button[data-tip=More]')!;
    more.click();
    const menu = () => container.querySelector('[role=menu] [data-popover]');
    expect(menu()).not.toBeNull();
    expect(more.getAttribute('aria-expanded')).toBe('true');
    document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
    expect(menu()).toBeNull();
    more.click();
    expect(menu()).not.toBeNull();
    document.body.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, composed: true }));
    expect(menu()).toBeNull();
  });

  it("gives a component's representations in its inspector their row's ⋯ menu (Remove among it)", () => {
    const { container, row } = mount();
    row('Polymer').click();
    const reprs = container.querySelector('[aria-label=Representations]')!;
    const more = reprs.querySelector<HTMLButtonElement>('button[aria-haspopup=menu]')!;
    expect(more).not.toBeNull();
    more.click();
    const menu = container.querySelector('[role=menu] [data-popover]');
    expect(menu?.textContent).toContain('Remove');
    expect(more.getAttribute('aria-expanded')).toBe('true');
    document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
    expect(container.querySelector('[role=menu] [data-popover]')).toBeNull();
  });

  it('renders the scene as a tree; hidden nodes keep their eye; collapsing hides children', () => {
    const { rows, row, scene } = mount();
    expect(row('1TQN')).toBeDefined();
    expect(row('Water').querySelector('button[data-tip=Show]')).not.toBeNull();
    const before = rows().length;
    scene.roots.value[0].toggleExpanded();
    flush();
    expect(rows().length).toBeLessThan(before);
  });

  it("shows and hides a node from its inspector's header, not only by its eye", () => {
    const { container, row } = mount();
    row('Water').click();
    const eye = () => container.querySelector<HTMLButtonElement>('[aria-label="Inspector"] header button')!;
    expect(eye().dataset.tip).toBe('Show');
    expect(eye().getAttribute('aria-pressed')).toBe('true');
    eye().click();
    expect(row('Water').querySelector('button[data-tip=Hide]')).not.toBeNull();
    expect(eye().dataset.tip).toBe('Hide');
    eye().click();
    expect(row('Water').querySelector('button[data-tip=Show]')).not.toBeNull();
  });

  it("shows an extension's row in the inspector when picked (its controls, turned on and off at its header), not under it", () => {
    const scene = mockScene();
    const shell = new ShellModel();
    const shown = observableState(false);
    shell.registerSections([
      {
        id: 'x.row.membrane',
        slot: SceneRowSlot,
        title: 'Membrane',
        when: (ctx) => (ctx.subject as { kind: string }).kind === 'structure',
        render: domView(() => (
          <SceneRow
            label="Membrane"
            description={observableState('Off').asReadonly}
            shown={shown.asReadonly}
            onToggle={() => shown.set(!shown.value)}
            body={() => <button type="button">Fit…</button>}
          />
        )),
      },
    ]);
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <StructureTab shell={shell} scene={scene} />, container);
    const row = shownInOrder(container, '[data-tree-row]').find((r) => r.textContent?.startsWith('Membrane'))!;
    row.click();
    const inspector = container.querySelector('[aria-label="Inspector"]')!;
    expect(inspector.textContent).toContain('Membrane');
    expect(row.parentElement!.textContent).not.toContain('Fit…');
    expect(inspector.textContent).toContain('Fit…');
    const eye = inspector.querySelector<HTMLButtonElement>('header button')!;
    expect(eye.dataset.tip).toBe('Turn on');
    eye.click();
    expect(shown.value).toBe(true);
    expect(eye.dataset.tip).toBe('Turn off');
    eye.click();
    expect(shown.value).toBe(false);
    // A node picked takes its place.
    scene.pick(scene.roots.value[0]);
    expect(container.querySelector('[aria-label="Inspector"]')!.textContent).not.toContain('Fit…');
  });

  it("lets an extension row's inspector go with the row (its structure removed)", () => {
    const scene = mockScene();
    const shell = new ShellModel();
    shell.registerSections([
      {
        id: 'x.row.membrane',
        slot: SceneRowSlot,
        title: 'Membrane',
        when: (ctx) => (ctx.subject as { kind: string }).kind === 'structure',
        render: domView(() => <SceneRow label="Membrane" body={() => <button type="button">Fit…</button>} />),
      },
    ]);
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <StructureTab shell={shell} scene={scene} />, container);
    const rows = () => shownInOrder(container, '[data-tree-row]').filter((r) => r.textContent?.startsWith('Membrane'));
    const before = rows().length;
    rows()[0].click();
    const inspector = () => container.querySelector('[aria-label="Inspector"]')!;
    expect(inspector().textContent).toContain('Fit…');
    // Its structure goes, and the row with it: its controls would act on nothing.
    scene.roots.value[0].remove();
    flush();
    expect(rows().length).toBe(before - 1);
    expect(inspector().textContent).not.toContain('Fit…');
    expect(inspector().textContent).toContain('Select an item in the tree');
  });

  it('follows an inspection built anew (its facts), keeping its form as it is', () => {
    const base = mockScene();
    const updates = observableState<SceneInspection | undefined>(undefined);
    const scene = Object.assign(base, { inspection: updates.asReadonly });
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <StructureTab shell={new ShellModel()} scene={scene} />, container);
    const row = shownInOrder(container, '[data-tree-row]').find((r) => r.textContent?.startsWith('Cartoon'))!;
    row.click();
    const picked = scene.picked.value!;
    const inspector = () => container.querySelector('[aria-label="Inspector"]')!;
    const type = inspector().querySelector('button[aria-haspopup=listbox][aria-label=Type]');
    expect(type).not.toBeNull();
    updates.set({ ...scene.inspect(picked), facts: [['Atoms', '4,242']] });
    expect(inspector().textContent).toContain('4,242');
    expect(inspector().querySelector('button[aria-haspopup=listbox][aria-label=Type]')).toBe(type);
  });

  it('picks a node into the inspector, and form edits change the tree', () => {
    const { container, row } = mount();
    row('Cartoon').click();
    // A choice expands inline: open it, pick an option.
    const choose = (label: string, option: string) => {
      const button = container.querySelector(
        `button[aria-haspopup=listbox][aria-label="${label}"]`,
      ) as HTMLButtonElement;
      button.click();
      const pick = [...container.querySelectorAll('[role=listbox] [role=option]')].find(
        (o) => o.textContent === option,
      ) as HTMLButtonElement;
      pick.click();
    };
    const type = container.querySelector('button[aria-haspopup=listbox][aria-label=Type]') as HTMLButtonElement;
    expect(type.textContent).toBe('Cartoon');
    choose('Type', 'Spacefill');
    expect(row('Spacefill')).toBeDefined();
    expect(container.querySelector('[role=listbox]')).toBeNull();
    // The mapped color param: switching to Uniform, its sliders button opens its own params (a color; v5's).
    choose('Color', 'Uniform');
    (container.querySelector('button[aria-label="Color settings"]') as HTMLButtonElement).click();
    // A color row: the swatch button opens the palette under it.
    const swatch = [...container.querySelectorAll('[aria-label=Inspector] button[aria-expanded]')].find(
      (b) => b.textContent?.trim() === '46, 134, 193',
    ) as HTMLButtonElement;
    expect(swatch.textContent?.trim()).toBe('46, 134, 193');
    swatch.click();
    expect(container.querySelectorAll('[role=listbox][aria-label=Palette] [role=option]').length).toBe(8);
  });

  it('toggles visibility inline without picking', () => {
    const { row, scene } = mount();
    const polymer = row('Polymer');
    (polymer.querySelector('button[data-tip=Hide]') as HTMLButtonElement).click();
    expect(row('Polymer').querySelector('button[data-tip=Show]')).not.toBeNull();
    expect(scene.picked.value).toBeUndefined();
  });

  it("shows particle lists as the tree's own rows: their representations and unit cell under them", () => {
    const { rows, dispose } = mount('particles');
    const labels = rows().map((r) => r.querySelector('span')?.textContent);
    expect(labels).toEqual(
      expect.arrayContaining(['cryoet-ribosomes.ndjson', 'Spacefill', 'Orientation', 'Unit cell', 'small.simularium']),
    );
    dispose();
  });

  it('filtered by the selection: only the rows that hold it, a chip in the search; its × or Backspace shows all again', async () => {
    const on = observableState(true);
    const held = new Set(['Ligand', 'Ball & stick']);
    const { rows, container, dispose } = mount(undefined, (scene) => ({
      ...scene,
      selectionFilter: {
        on: on.asReadonly,
        version: observableState(0).asReadonly,
        set: (v: boolean) => on.set(v),
        holds: (n: SceneNode) => held.has(n.label.value),
      },
    }));
    const labels = () => rows().map((r) => r.querySelector('span')?.textContent);
    await new Promise((r) => setTimeout(r, 0));
    expect(labels()).toContain('Ligand');
    expect(labels()).not.toContain('Water');
    const chip = () => container.querySelector('[aria-label="Show every row, not only those in the selection"]');
    expect(chip()).not.toBeNull();

    (chip() as HTMLButtonElement).click();
    await new Promise((r) => setTimeout(r, 0));
    expect(on.value).toBe(false);
    expect(labels()).toContain('Water');

    on.set(true);
    await new Promise((r) => setTimeout(r, 0));
    const input = container.querySelector('input')!;
    input.dispatchEvent(new KeyboardEvent('keydown', { key: 'Backspace', bubbles: true }));
    expect(on.value).toBe(false);
    dispose();
  });
});
