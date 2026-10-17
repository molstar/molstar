/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { afterAll, beforeAll, describe, expect, it, vi } from 'vitest';
import { ensemble } from '../../api/_test/trajectory-fixtures.js';
import { cbs, expectOk } from '../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../default-spec.js';
import { createDefaultPlugin } from '../../sources/create-plugin.js';
import { createViewerCore, type ViewerCore } from '../../viewer.js';
import { type MountedShell, mountShell } from '../shell.js';
import { render } from '@molstar/prism-ui/dom/render';
import { bind } from '@molstar/prism-ui/dom/scope';
import { setText } from '@molstar/prism-ui/dom/writes';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { type Audit, installAudit } from './waste-audit.js';
import { personSession } from './waste-session.js';

// The waste audit over a person's session in the real shell (view-models.md §3: a view follows the narrowest source
// that decides its output). `WASTE_REPORT=1` prints the table. The budget (waste-budget.json): per source file, the
// writes and renders that changed no DOM, and the derived recomputes that came out the same (the cheap gates that keep
// the rest from waking: budgeted apart). A source over its budget fails; `WASTE_UPDATE=1` writes today's numbers (after
// an intended change; they should only go down).

class NoResize {
  observe() {}
  unobserve() {}
  disconnect() {}
}

let core: ViewerCore;
let mounted: MountedShell;
let host: HTMLElement;
let audit: ReturnType<typeof installAudit>;

beforeAll(async () => {
  vi.stubGlobal('ResizeObserver', NoResize);
  audit = installAudit();
  core = await createViewerCore({
    spec: DefaultPrismViewerSpec(),
    createPlugin: createDefaultPlugin,
    modules: ['files', 'page'],
  });
  vi.spyOn(core, 'mountCanvas').mockResolvedValue(false);
  host = document.createElement('div');
  document.body.append(host);
  mounted = await mountShell(core, host);
  audit.observe(mounted.root.getRootNode());
  audit.observe(document.body);
  expectOk(await core.api.load.files({ files: [cbs()] }));
  expectOk(await core.api.load.files({ files: [ensemble(5)] }));
  await new Promise((r) => setTimeout(r, 400));
});

afterAll(async () => {
  audit.dispose();
  mounted.dispose();
  await core.dispose();
  host.remove();
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
});

export function report(a: Audit, top = 25) {
  const pad = (s: string | number, n: number) => String(s).padStart(n);
  const lines = ['source state | changes | listener calls | derives (same) | writes (no DOM) | renders (no DOM)'];
  for (const r of a.sources().slice(0, top))
    lines.push(
      `${r.source} | ${pad(r.changes, 4)} | ${pad(r.listenerCalls, 5)} | ${pad(r.derives, 4)} (${r.derivesSame}) | ${pad(r.writes, 4)} (${r.writesNoDom}) | ${pad(r.renders, 4)} (${r.rendersNoDom})`,
    );
  lines.push('', 'view | calls | wasted');
  for (const v of a.views().slice(0, top)) lines.push(`${v.view} | ${v.calls} | ${v.wasted}`);
  lines.push('', 'source -> view | calls | wasted');
  for (const v of a.pairs().slice(0, top)) lines.push(`${v.view} | ${v.calls} | ${v.wasted}`);
  lines.push('', `wasted in all: ${a.wasted()}`);
  return lines.join('\n');
}

const BudgetFile = path.join(import.meta.dirname, 'waste-budget.json');

interface Budget {
  /** Writes and renders that changed no DOM. */
  noDom: number;
  /** Derived recomputes that came out the same. */
  sameDerives: number;
}

/** Per source file (a line moves as code is edited; a file's waste is what is budgeted). */
function bySourceFile(a: Audit): Record<string, Budget> {
  const out: Record<string, Budget> = {};
  for (const r of a.sources()) {
    const file = r.source.replace(/:\d+$/, '');
    const b = (out[file] ??= { noDom: 0, sameDerives: 0 });
    b.noDom += r.writesNoDom + r.rendersNoDom;
    b.sameDerives += r.derivesSame;
  }
  return out;
}

describe('the waste audit', () => {
  it("a person's session wastes no more than its budget, per source", async () => {
    audit.start();
    await personSession(core, mounted);
    audit.stop();
    if (process.env.WASTE_REPORT) console.log(report(audit, 200));
    expect(audit.sources().length).toBeGreaterThan(0);
    const now = bySourceFile(audit);
    if (process.env.WASTE_UPDATE) {
      const sorted = Object.fromEntries(Object.entries(now).sort(([a], [b]) => a.localeCompare(b)));
      fs.writeFileSync(BudgetFile, `${JSON.stringify(sorted, null, 2)}\n`);
      return;
    }
    const budget = JSON.parse(fs.readFileSync(BudgetFile, 'utf8')) as Record<string, Budget>;
    const over: string[] = [];
    for (const [file, b] of Object.entries(now)) {
      const allowed = budget[file] ?? { noDom: 0, sameDerives: 0 };
      if (b.noDom > allowed.noDom)
        over.push(`${file}: ${b.noDom} writes/renders with no DOM change (budget ${allowed.noDom})`);
      if (b.sameDerives > allowed.sameDerives)
        over.push(`${file}: ${b.sameDerives} derives that came out the same (budget ${allowed.sameDerives})`);
    }
    expect(over, 'over the waste budget: a view follows more than decides its output (view-models.md §3)').toEqual([]);
  });

  it('catches waste (the check can fail): a binding told of every change, for a value that never changes', async () => {
    const noisy = observableState(0);
    const box = document.createElement('div');
    document.body.append(box);
    const text = document.createTextNode('the same');
    box.append(text);
    // Follows the whole state for a text that stays: each change is a write that changes no DOM.
    const dispose = render(() => {
      bind(noisy.asReadonly, () => setText(text, 'the same'));
      return null;
    }, box);
    audit.start();
    for (let i = 1; i <= 5; i++) noisy.set(i);
    audit.stop();
    dispose();
    box.remove();
    expect(audit.wasted()).toBe(5);
  });
});
