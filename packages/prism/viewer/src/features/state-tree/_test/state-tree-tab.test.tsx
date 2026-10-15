/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { render } from '@molstar/prism-ui/dom/render';
import { afterEach, describe, expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer, until } from '../../../_test/harness.js';
import { StateTreeTab } from '../../../views/state-tree/state-tree-tab.js';
import { StateTreeFeature } from '../index.js';
import { LiveStateTreeModel } from '../model.js';
import type { StateCellNode } from '../source.js';

// The State tree tab over a real (headless) viewer: the tree, the eye, the update forms and the actions.

afterEach(() => {
  document.body.innerHTML = '';
});

async function mount() {
  const v = await headlessViewer({ extensions: [StateTreeFeature], files: [cbs()] });
  const model = new LiveStateTreeModel(v.api.state!);
  const container = document.createElement('div');
  document.body.append(container);
  // As `ui.tsx` wires it while the tab shows.
  const follow = v.api.state!.tree.changed.subscribe((refs) => model.refresh(refs));
  const unmount = render(() => <StateTreeTab model={model} />, container);
  const dispose = () => {
    follow.unsubscribe();
    unmount();
  };
  const rows = () => [...container.querySelectorAll('[role=treeitem] > div:first-child')] as HTMLElement[];
  const row = (text: string) => rows().find((r) => r.querySelector('span')?.textContent === text)!;
  return { v, model, container, dispose, rows, row };
}

const flat = (nodes: readonly StateCellNode[]): StateCellNode[] => nodes.flatMap((n) => [n, ...flat(n.children ?? [])]);

describe('StateTreeTab over the viewer', () => {
  it('renders the cells, and the eye hides a representation', async () => {
    const { v, row, model, dispose } = await mount();
    expect(row('Cartoon')).toBeDefined();
    expect(row('Water')).toBeDefined();
    const water = flat(model.roots.value as StateCellNode[]).find((n) => n.label === 'Water')!;
    const rep = water.children![0];
    const eye = row('Ball & Stick').querySelector('button[title=Hide]');
    expect(eye).not.toBeNull();
    // The first Ball & Stick is the ligand's; use the model for the water's.
    model.toggleVisible(rep.ref);
    await until(() => v.plugin.state.data.cells.get(rep.ref)?.state.isHidden === true);
    await until(() => document.querySelectorAll('button[title=Show]').length > 0);
    dispose();
  });

  it('selects a cell into the inspector: its update form edits the cell, its actions apply', async () => {
    const { v, model, row, container, dispose } = await mount();
    const structure = flat(model.roots.value as StateCellNode[]).find((n) => n.transform === 'Structure')!;
    row(structure.label).click();
    expect(container.textContent).toContain('Update Structure');
    expect(container.textContent).toContain('Actions');

    const inspection = model.inspect(structure.ref)!;
    expect(inspection.updates[0].form!.params).toHaveProperty('type');
    expect(inspection.actions.map((a) => a.id)).toContain('selection');

    const cells = v.plugin.state.data.cells.size;
    await inspection.actions.find((a) => a.id === 'selection')!.apply();
    expect(v.plugin.state.data.cells.size).toBe(cells + 1);
    await until(() => flat(model.roots.value as StateCellNode[]).some((n) => n.transform === 'Selection'));
    dispose();
  });

  it('applies a form edit through state.update, debounced, and follows an Undo', async () => {
    const { v, model, dispose } = await mount();
    const nodes = flat(model.roots.value as StateCellNode[]);
    const rep = nodes.find((n) => n.label === 'Cartoon')!;
    const inspection = model.inspect(rep.ref)!;
    const form = inspection.updates[0].form!;
    const params = () => v.plugin.state.data.cells.get(rep.ref)!.transform.params as any;
    const colorTheme = params().colorTheme;
    expect(form.values.value.colorTheme).toEqual(colorTheme);

    form.set('colorTheme', { name: 'secondary-structure', params: {} });
    // Not at once: the edit waits for the next one.
    expect(params().colorTheme.name).toBe('chain-id');
    await until(() => params().colorTheme.name === 'secondary-structure');
    // The cell changes before its run ends: Undo waits for the step.
    await until(() => v.managers.history.canUndo.value);

    expectOk(await v.api.history.undo());
    expect(params().colorTheme.name).toBe('chain-id');
    // The form was told: it shows the cell again.
    await until(() => (form.values.value.colorTheme as any).name === 'chain-id');
    dispose();
  });

  it('removes a cell through the model, and the selection goes with it', async () => {
    const { v, model, dispose } = await mount();
    const ligand = flat(model.roots.value as StateCellNode[]).find((n) => n.label === 'Ligand')!;
    model.select(ligand.ref);
    model.remove(ligand.ref);
    await until(() => !v.plugin.state.data.cells.has(ligand.ref));
    model.prune();
    expect(model.selected.value).toBeUndefined();
    dispose();
  });
});
