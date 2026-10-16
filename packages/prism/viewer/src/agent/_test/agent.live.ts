/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { headlessViewer } from '../../_test/harness.js';
import type { PrismViewer } from '../../viewer.js';
import { createViewerAgent } from '../session.js';
import { env, LiveModels } from './live-models.js';

// The agent against real models (`pnpm test:live`, never in `test:unit`): a few first-scope tasks on a headless viewer,
// checked on the scene. It reaches the network (the models, and the PDB for the entry it loads).

const reprTypes = (v: PrismViewer) =>
  v.plugin.managers.structure.hierarchy.current.structures.flatMap((s) =>
    s.components.flatMap((c) =>
      c.representations.map((r) => (r.cell.transform.params as { type: { name: string } }).type.name),
    ),
  );

for (const m of LiveModels) {
  describe.skipIf(!env[m.key])(`the agent on ${m.model}`, () => {
    it('loads an entry, styles it, selects and focuses the ligand; one undoable turn each', async () => {
      const v = await headlessViewer();
      const provider = m.create();
      const agent = createViewerAgent({
        viewer: v,
        models: {
          providerFor: (p) => (p === m.provider ? provider : undefined),
          default: { provider: m.provider, model: m.model },
        },
      });
      const log: string[] = [];
      const hooks = { onTool: (e: { name: string; status: string }) => log.push(`${e.name}:${e.status}`) };

      const first = await agent.turn('Load PDB entry 1CBS and show it as a surface.', hooks);
      expect(first.status, log.join(' ')).toBe('completed');
      expect(v.plugin.managers.structure.hierarchy.current.structures).toHaveLength(1);
      expect(reprTypes(v).some((t) => t.includes('surface'))).toBe(true);

      const second = await agent.turn('Select the ligand and focus on it.', hooks);
      expect(second.status, log.join(' ')).toBe('completed');
      expect(v.plugin.managers.structure.selection.elementCount()).toBeGreaterThan(0);

      // The whole of the first turn is one step: the person can take it back at once.
      expect(v.managers.history.steps.value.some((s) => s.origin.kind === 'agent')).toBe(true);
      console.info(`${m.model}: ${log.join(' ')}; usage ${JSON.stringify(second.usage ?? {})}`);
      agent.dispose();
    });
  });
}
