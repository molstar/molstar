/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html';
import { setDebugMode } from '@molstar/core/util/debug';
import { AllCases, OnDemandCases } from './cases/index.js';
import { getBackendProfile, ReferenceProfileId } from './framework/backend.js';
import { encodePngDataUrl } from './framework/image-io.js';
import { DefaultSnapshotTolerance, summaryToJson, VisualTestRunner } from './framework/runner.js';
import { matchesAny } from './framework/test-case.js';
import { Gallery } from './ui/gallery.js';

declare global {
  interface Window {
    /** JSON summary for automation, updated after each test */
    __visualTestResults?: ReturnType<typeof summaryToJson>;
    __visualTestImage?: (id: string, kind: 'reference' | 'target' | 'diff') => Promise<string | undefined>;
    /** The case running, for a runner that reports a stalled page */
    __visualTestCurrent?: string;
  }
}

function splitList(value: string | null) {
  return (value ?? '')
    .split(',')
    .map((s) => s.trim())
    .filter((s) => s.length > 0);
}

async function main() {
  const params = new URLSearchParams(window.location.search);
  const targetId = params.get('target') || ReferenceProfileId;
  const referenceId = params.get('reference') || ReferenceProfileId;
  const filter = splitList(params.get('filter'));
  const tags = splitList(params.get('tags'));
  const snapshots = params.get('snapshots') === '1';
  const zoom = Number(params.get('zoom')) || 1;
  // enables GL error and program/uniform consistency checks, slows rendering down
  if (params.get('debug') === '1') setDebugMode(true);
  // the dev server serves the repository root
  const referencesUrl =
    params.get('references') || new URL('../../../src/tests/visual/references/', window.location.href).href;

  const target = getBackendProfile(targetId);
  const reference = getBackendProfile(referenceId);
  if (!target || !reference) throw new Error(`Unknown backend profile '${!target ? targetId : referenceId}'`);

  // `shard=i/n`: every n-th of the selected cases from the i-th (0-based), to run a suite in several pages
  const [shardIndex, shardCount] = (params.get('shard') ?? '0/1').split('/').map(Number);
  // the on-demand cases only where a filter or tags select cases
  const candidates = filter.length > 0 || tags.length > 0 ? [...AllCases, ...OnDemandCases] : AllCases;
  const tests = candidates
    .filter(
      (c) =>
        (filter.length === 0 || matchesAny(c.id, filter)) &&
        (tags.length === 0 || tags.some((t) => c.tags?.includes(t))),
    )
    .filter((_, i) => i % shardCount === shardIndex);

  const runner = new VisualTestRunner({
    target,
    reference,
    snapshots,
    snapshotTolerance: DefaultSnapshotTolerance,
    referencesUrl,
    sharedReferencesUrl: params.get('shared') ?? undefined,
    drawLists: params.get('drawLists') === '1',
  });

  const rerun = (id: string) => runner.runTest(tests.find((t) => t.id === id)!);
  const gallery = new Gallery(
    document.getElementById('app')!,
    {
      target: targetId,
      reference: referenceId,
      filter: filter.join(', '),
      tags: tags.join(', '),
      snapshots,
      zoom,
    },
    runner.summary,
    rerun,
  );

  window.__visualTestImage = async (id, kind) => {
    const r = runner.summary.results.find((r) => r.id === id);
    const image = kind === 'diff' ? r?.diff : r?.[kind]?.image;
    return image ? encodePngDataUrl(image) : undefined;
  };

  await runner.run(
    tests,
    (result) => {
      if (runner.summary.results.length === 1) gallery.setReports(runner.summary.reports);
      gallery.addResult(result);
      window.__visualTestResults = summaryToJson(runner.summary);
    },
    (test) => {
      window.__visualTestCurrent = test.id;
    },
  );
  gallery.setReports(runner.summary.reports);
  gallery.updateSummary();
  window.__visualTestResults = summaryToJson(runner.summary);
}

main().catch((e) => {
  console.error(e);
  document.body.append(Object.assign(document.createElement('pre'), { textContent: `${e}` }));
});
