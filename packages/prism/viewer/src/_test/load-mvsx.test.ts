/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { Zip } from '@molstar/core/util/zip/zip';
import { afterEach, describe, expect, it } from 'vitest';
import { DefaultPrismViewerSpec } from '../default-spec.js';
import { createPrismViewerHeadless } from '../headless.js';
import type { PrismViewer } from '../viewer.js';

// MVSX archives opened one after another (`load.files`). An MVS view replaces the scene, and its node refs are
// deterministic (`!mvs:<hash>:<n>`), so the second view's nodes can take the first one's refs.

const dataDir = path.resolve(import.meta.dirname, '../../../../../data/examples');

/** An MVSX archive of one structure file from the examples, shown as a cartoon colored `color`. */
async function mvsx(name: string, file: string, color: string) {
  const index = {
    root: {
      kind: 'root',
      children: [
        {
          kind: 'download',
          params: { url: file },
          children: [
            {
              kind: 'parse',
              params: { format: 'mmcif' },
              children: [
                {
                  kind: 'structure',
                  params: { type: 'model' },
                  children: [
                    {
                      kind: 'component',
                      params: { selector: 'polymer' },
                      children: [
                        {
                          kind: 'representation',
                          params: { type: 'cartoon' },
                          children: [{ kind: 'color', params: { color } }],
                        },
                      ],
                    },
                  ],
                },
              ],
            },
          ],
        },
      ],
    },
    metadata: { version: '1.0', timestamp: '2026-10-20T00:00:00Z' },
  };
  const bytes = await Zip({
    'index.mvsj': new TextEncoder().encode(JSON.stringify(index)),
    [file]: new Uint8Array(fs.readFileSync(path.join(dataDir, file))),
  }).run();
  return new File([bytes], name);
}

let viewer: PrismViewer | undefined;
afterEach(async () => {
  await viewer?.dispose();
  viewer = undefined;
});

/** The plugin's errors from now on. */
function errorsOf(v: PrismViewer) {
  const errors: string[] = [];
  v.plugin.events.log.subscribe((e) => {
    if (e.type === 'error') errors.push(e.message);
  });
  return errors;
}

describe('MVSX files opened in a row', () => {
  it('opens the second one, which replaces the first, and keeps it', async () => {
    viewer = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec() });
    const errors = errorsOf(viewer);
    const first = await viewer.api.load.files({ files: [await mvsx('a.mvsx', '1crn.cif', 'red')] as never });
    expect(first).toMatchObject({ ok: true, did: 'Opened a.mvsx (1 structure).' });
    // Its root nodes take the first view's refs: judged by ref alone it "made nothing", failed, and the step was
    // rolled back to the first view, whose archive's files the second had replaced (`fetch failed`).
    const second = await viewer.api.load.files({ files: [await mvsx('b.mvsx', '1cbs_updated.cif', 'blue')] as never });
    expect(second).toMatchObject({ ok: true, did: 'Opened b.mvsx (1 structure).' });
    const listed = await viewer.api.structure.list();
    expect(listed.ok && listed.did).toMatch(/1CBS/);
    expect(viewer.api.structure.loaded.value).toHaveLength(1);
    expect(errors).toEqual([]);
  });

  it('opens the same one again', async () => {
    viewer = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec() });
    const file = await mvsx('a.mvsx', '1crn.cif', 'red');
    expect(await viewer.api.load.files({ files: [file] as never })).toMatchObject({ ok: true });
    expect(await viewer.api.load.files({ files: [file] as never })).toMatchObject({ ok: true });
    expect(viewer.api.structure.loaded.value).toHaveLength(1);
  });

  it('still fails a file no reader takes', async () => {
    viewer = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec() });
    expect(await viewer.api.load.files({ files: [await mvsx('a.mvsx', '1crn.cif', 'red')] as never })).toMatchObject({
      ok: true,
    });
    const r = await viewer.api.load.files({ files: [new File(['x'], 'notes.unknownformat')] as never });
    expect(r).toMatchObject({ ok: false, code: 'failed' });
    expect(viewer.api.structure.loaded.value).toHaveLength(1);
  });
});
