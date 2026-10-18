/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RuntimeProbe, type Scope } from '@molstar/prism-ui/dom/scope';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { afterAll, beforeAll, describe, expect, it, vi } from 'vitest';
import { ensemble } from '../../api/_test/trajectory-fixtures.js';
import { cbs, expectOk } from '../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../default-spec.js';
import { createDefaultPlugin } from '../../sources/create-plugin.js';
import { createViewerCore, type ViewerCore } from '../../viewer.js';
import { type MountedShell, mountShell } from '../shell.js';

// Hidden UI does no UI updates (styling.md §5): with one tab shown, the others, the folded sections and the folded
// toolbar are inert and paused (`activity`). The model changes a person makes (a load, picks, hiding and showing,
// renames, selections, trajectory frames) write nothing inside them but the `inert` toggles; the shown tab updates.
// The runtime's probe counts binding writes and `watch` renders in paused scopes (none). A binding that does not
// pause, planted in a hidden panel, is caught: the check can fail.

class NoResize {
  observe() {}
  unobserve() {}
  disconnect() {}
}

const wait = (ms = 60) => new Promise((r) => setTimeout(r, ms));

/** Where a node is: its element's tag and classes, and the hidden part it is in (a panel, a section). */
function where(node: Node): string {
  const el = node.nodeType === Node.ELEMENT_NODE ? (node as Element) : node.parentElement;
  if (!el) return '?';
  const hidden = el.closest('[inert]');
  const name = (e: Element | null) =>
    e
      ? `${e.tagName.toLowerCase()}${e.getAttribute('role') ? `[${e.getAttribute('role')}]` : ''}${e.getAttribute('aria-label') ? `"${e.getAttribute('aria-label')}"` : ''}.${(e.getAttribute('class') ?? '').split(' ')[0]}`
      : '';
  return `${name(el)} in ${name(hidden)}`;
}

/** A change as seen when delivered: whether it was inside an inert (hidden) part, and where. */
interface Seen {
  readonly record: MutationRecord;
  readonly hidden: boolean;
  readonly where: string;
}

/** The changes inside inert (hidden) parts, but the `inert` toggles themselves. */
function hiddenMutations(seen: readonly Seen[]): string[] {
  return seen
    .filter((s) => s.hidden && !(s.record.type === 'attributes' && s.record.attributeName === 'inert'))
    .map((s) => `${s.record.type}${s.record.attributeName ? ` ${s.record.attributeName}` : ''}: ${s.where}`);
}

let core: ViewerCore;
let mounted: MountedShell;
let host: HTMLElement;
const counts = { write: 0, render: 0, hiddenWrite: 0, hiddenRender: 0 };

beforeAll(async () => {
  vi.stubGlobal('ResizeObserver', NoResize);
  RuntimeProbe.on = (kind, scope: Scope) => {
    counts[kind]++;
    if (!scope.active) counts[kind === 'write' ? 'hiddenWrite' : 'hiddenRender']++;
  };
  const base = DefaultPrismViewerSpec();
  core = await createViewerCore({
    // The toolbar folded: its tools are hidden UI too.
    spec: { ...base, ui: { ...base.ui, layout: { ...base.ui?.layout, toolbar: 'collapsed' } } },
    createPlugin: createDefaultPlugin,
    modules: ['files', 'page'],
  });
  vi.spyOn(core, 'mountCanvas').mockResolvedValue(false);
  host = document.createElement('div');
  document.body.append(host);
  mounted = await mountShell(core, host);
  expectOk(await core.api.load.files({ files: [cbs()] }));
  await wait(400);
});

afterAll(async () => {
  RuntimeProbe.on = undefined;
  mounted.dispose();
  await core.dispose();
  host.remove();
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
});

/** Records the shell's mutations while `drive` runs (and a moment after). */
async function record(drive: () => Promise<void>): Promise<Seen[]> {
  const records: Seen[] = [];
  // Classified when delivered (a microtask after the change), while the node is where it was written.
  const see = (list: MutationRecord[]) => {
    for (const r of list) {
      const el = r.target.nodeType === Node.ELEMENT_NODE ? (r.target as Element) : r.target.parentElement;
      records.push({ record: r, hidden: !!el?.closest('[inert]'), where: where(r.target) });
    }
  };
  const observer = new MutationObserver(see);
  observer.observe(mounted.root.getRootNode(), {
    subtree: true,
    childList: true,
    attributes: true,
    characterData: true,
  });
  await drive();
  await wait(300);
  see(observer.takeRecords());
  observer.disconnect();
  return records;
}

/** The model changes a person makes: a load, picks, a selection, hiding and showing, a rename, trajectory frames. */
async function drive(round: number) {
  const api = core.api.as({ kind: 'person' });
  const scene = core.views.scene;
  expectOk(await api.load.files({ files: [ensemble(5)] }));
  await wait(300);
  const roots = scene.roots.value;
  scene.pick(roots[0]);
  await wait(50);
  scene.pick(roots[roots.length - 1]);
  await wait(50);
  expectOk(await api.select.set({ target: 'A:1-20', structure: '1cbs' }));
  await wait(30);
  expectOk(await api.select.clear());
  expectOk(await api.component.visibility({ component: 'Polymer', structure: '1cbs', visible: false }));
  expectOk(await api.component.visibility({ component: 'Polymer', structure: '1cbs', visible: true }));
  expectOk(await api.component.rename({ component: 'Polymer', structure: '1cbs', label: `Chains ${round}` }));
  await wait(30);
  expectOk(await api.component.rename({ component: `Chains ${round}`, structure: '1cbs', label: 'Polymer' }));
  for (const index of [1, 2, 3, 0]) {
    expectOk(await api.trajectory.frame({ index }));
    await wait(30);
  }
  // Gone again: the next round's frames have one ensemble to step.
  expectOk(await api.structure.remove({ structure: 'ensemble.pdb' }));
  await wait(100);
}

const zero = () => {
  for (const k of Object.keys(counts) as (keyof typeof counts)[]) counts[k] = 0;
};

describe('hidden UI', () => {
  const setups: [string, () => Promise<void>, string][] = [
    [
      'the Structure tab shown, the toolbar folded',
      async () => {
        mounted.shell.show('structure');
      },
      'Structure',
    ],
    [
      'the Selection tab shown, its sections folded',
      async () => {
        mounted.shell.show('selection');
      },
      'Selection',
    ],
    [
      'the rail compact and closed: the docks behind it',
      async () => {
        mounted.shell.show('structure');
        core.managers.viewState.setCompactRail(true);
      },
      '',
    ],
  ];
  let round = 0;
  for (const [name, setup, shownTab] of setups) {
    it(`${name}: nothing hidden is written while the model changes; what shows updates`, async () => {
      await setup();
      await wait(400);
      const inert = [...mounted.root.querySelectorAll('[inert]')];
      expect(inert.length).toBeGreaterThan(0);
      zero();
      const records = await record(() => drive(++round));
      expect(hiddenMutations(records)).toEqual([]);
      expect({ writes: counts.hiddenWrite, renders: counts.hiddenRender }).toEqual({ writes: 0, renders: 0 });
      expect(counts.write + counts.render).toBeGreaterThan(0);
      if (shownTab) {
        const shown = [...mounted.root.querySelectorAll('[role=tabpanel]:not([inert])')].find(
          (p) => p.getAttribute('aria-label') === shownTab,
        );
        expect(shown).toBeDefined();
        expect(records.filter((r) => shown!.contains(r.record.target)).length).toBeGreaterThan(0);
      }
    });
  }

  it('catches a binding that does not pause (the check can fail)', async () => {
    const hidden = mounted.root.querySelector<HTMLElement>('[inert]')!;
    const leak = observableState('a');
    const span = document.createElement('span');
    hidden.append(span);
    // Not `bind`: a raw subscription, which no scope pauses.
    const sub = leak.asReadonly.subscribe((v) => {
      span.textContent = v;
    });
    const records = await record(async () => {
      leak.set('b');
    });
    sub.unsubscribe();
    span.remove();
    expect(hiddenMutations(records).length).toBeGreaterThan(0);
  });
});
