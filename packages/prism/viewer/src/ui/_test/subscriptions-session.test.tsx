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
import { bind } from '@molstar/prism-ui/dom/scope';
import { describe, expect, it, vi } from 'vitest';
import { cbs, expectOk } from '../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../default-spec.js';
import { createDefaultPlugin } from '../../sources/create-plugin.js';
import { createViewerCore } from '../../viewer.js';
import { mountShell } from '../shell.js';
import { installAudit } from './waste-audit.js';

// The app's subscriptions scale with what is shown (view-models.md §3), in the real shell: with each tab shown in turn,
// the live listeners by the view that added them (the subscription inventory, waste-audit.ts) are the same with 24
// structures loaded as with 48 (both fill the trees' screens). A view that listens to every structure (the positive
// control, planted after the tabs) is caught. `SUBS_REPORT=<file>` writes the inventory.

class NoResize {
  observe() {}
  unobserve() {}
  disconnect() {}
}

const Tabs = ['structure', 'state-tree.tab.tree', 'selection', 'history'] as const;
const Control = 'control: a view listening to every structure';
const wait = (ms = 30) => new Promise((r) => setTimeout(r, ms));

/** Live listeners by the file that added them, with each tab shown, and with the control planted. */
async function session(structures: number) {
  const audit = installAudit();
  const core = await createViewerCore({
    spec: DefaultPrismViewerSpec(),
    createPlugin: createDefaultPlugin,
    modules: ['files', 'page'],
  });
  vi.spyOn(core, 'mountCanvas').mockResolvedValue(false);
  const host = document.createElement('div');
  document.body.append(host);
  const mounted = await mountShell(core, host);
  for (let i = 0; i < structures; i++) expectOk(await core.api.load.files({ files: [cbs()] }));
  // The sequence band's chunk in (it comes with the first structure shown, at its own pace): both sessions measured
  // with it, then settled.
  await core.api.sequence.runtime();
  await audit.quiet(500, 10_000);
  const byFile = () => {
    const m = new Map<string, number>();
    for (const r of audit.subscriptions()) {
      const file = r.listener.replace(/:\d+$/, '');
      m.set(file, (m.get(file) ?? 0) + r.count);
    }
    return m;
  };
  const out = new Map<string, Map<string, number>>();
  for (const tab of Tabs) {
    mounted.shell.show(tab);
    await wait(150);
    out.set(tab, byFile());
  }
  // The control: a label per structure, each followed (not the narrowest source).
  const planted = document.createElement('div');
  document.body.append(planted);
  const unplant = render(() => {
    const el = document.createElement('div');
    for (const n of core.views.scene.roots.value) bind(n.label, (l) => (el.title = l));
    return el;
  }, planted);
  out.set(Control, byFile());
  unplant();
  planted.remove();
  if (process.env.SUBS_REPORT)
    (await import('node:fs')).writeFileSync(
      `${process.env.SUBS_REPORT}.${structures}`,
      audit
        .subscriptions()
        .map((r) => `${r.count} | ${r.state} <- ${r.listener}`)
        .join('\n'),
    );
  await audit.quiet(300, 5000);
  audit.dispose();
  mounted.dispose();
  await core.dispose();
  host.remove();
  return out;
}

/** The files whose listeners grew from `a` to `b`. */
function grown(a: Map<string, number> | undefined, b: Map<string, number> | undefined): string[] {
  const out: string[] = [];
  for (const [file, n] of b ?? []) {
    const was = a?.get(file) ?? 0;
    if (n > was) out.push(`${file}: ${was} -> ${n}`);
  }
  return out;
}

describe('the app’s subscriptions', () => {
  it('are the same for 24 structures as for 48, in each tab; a view listening to each structure is caught', async () => {
    vi.stubGlobal('ResizeObserver', NoResize);
    const a = await session(24);
    const b = await session(48);
    for (const tab of Tabs)
      expect(grown(a.get(tab), b.get(tab)), `${tab}: listeners that grew with the model`).toEqual([]);
    expect(grown(a.get(Control), b.get(Control)).length).toBeGreaterThan(0);
    vi.unstubAllGlobals();
  }, 120_000);
});
