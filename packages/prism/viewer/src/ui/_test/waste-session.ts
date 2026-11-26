/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { EmptyLoci } from '@molstar/model/model/loci';
import { StructureElement } from '@molstar/model/model/structure';
import { expectOk } from '../../_test/harness.js';
import type { ViewerCore } from '../../viewer.js';
import type { MountedShell } from '../shell.js';

// A person's session in the real shell, for the waste audit (waste.test.tsx): a pick, hover across the structure,
// select and clear, a representation's color, tabs opened and closed, the command menu typed in and closed, the theme,
// trajectory frames, notifications. Each step waits for what follows from it to settle: at least its time, then until
// no state has changed for a while (`settle`), so a loaded machine's late work does not land in the next step.

const wait = (ms = 30) => new Promise((r) => setTimeout(r, ms));

export async function personSession(
  core: ViewerCore,
  mounted: MountedShell,
  settle: () => Promise<void> = async () => {},
): Promise<void> {
  const step = async (ms = 30) => {
    await wait(ms);
    await settle();
  };
  const api = core.api.as({ kind: 'person' });
  const scene = core.views.scene;
  const shell = mounted.shell;
  const root = mounted.root;

  shell.show('structure');
  await step(100);
  // Picks.
  const nodes = [...scene.roots.value, ...scene.roots.value.flatMap((r) => r.children.value)];
  for (const n of nodes.slice(0, 6)) {
    scene.pick(n);
    await step();
  }
  scene.pick(undefined);
  await step();

  // Hover across the structure: residues in turn, then nothing.
  const structure = core.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
  if (structure) {
    const all = StructureElement.Loci.all(structure);
    const first = StructureElement.Loci.firstResidue(all);
    const hover = (loci: StructureElement.Loci | typeof EmptyLoci) =>
      core.plugin.behaviors.interaction.hover.next({
        current: { loci },
        buttons: 0,
        button: 0,
        modifiers: { shift: false, control: false, alt: false, meta: false },
      } as never);
    for (let i = 0; i < 20; i++) {
      hover(i % 3 === 0 ? first : i % 3 === 1 ? all : EmptyLoci);
      await wait(5);
    }
    hover(EmptyLoci);
    await step();
  }

  // Select and clear.
  for (const target of ['A:1-20', 'A:30-40']) {
    expectOk(await api.select.set({ target, structure: '1cbs' }));
    await step();
    expectOk(await api.select.clear());
    await step();
  }

  // A representation's color.
  expectOk(await api.representation.color({ target: 'A:10-20', theme: 'hydrophobicity', structure: '1cbs' } as never));
  await step(100);

  // Tabs.
  for (const tab of ['state-tree.tab.tree', 'history', 'settings', 'structure']) {
    shell.show(tab);
    await step(80);
  }

  // The command menu: opened, typed in, closed.
  root.querySelector<HTMLElement>('button[title="Command menu"], [aria-label="Command menu"]')?.click();
  await step(80);
  const input = root.querySelector<HTMLInputElement>('[role=dialog] input, [role=combobox]');
  if (input) {
    for (const text of ['c', 'co', 'col', 'colo', 'color']) {
      input.value = text;
      input.dispatchEvent(new Event('input', { bubbles: true, composed: true }));
      await wait(20);
    }
    input.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowDown', bubbles: true, composed: true }));
    input.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowDown', bubbles: true, composed: true }));
    await step();
    input.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true, composed: true }));
    await step(80);
  }

  // The theme.
  core.managers.commands.open('view.theme');
  await step(80);

  // Trajectory frames.
  for (const index of [1, 2, 3, 4, 0, 1, 2]) {
    expectOk(await api.trajectory.frame({ index }));
    await step(20);
  }

  // Notifications.
  const notes = core.managers.notifications;
  for (let i = 0; i < 3; i++) {
    notes.post(`Saved ${i}`);
    await step(20);
  }
  notes.putTask({ id: 'audit', label: 'Working' });
  await step(20);
  notes.removeTask('audit');
  await step(100);
}
