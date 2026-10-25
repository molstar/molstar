/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { setCanvasModule } from '@molstar/graphics/geo/geometry/text/font-atlas';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterAll, afterEach, beforeAll, describe, expect, it } from 'vitest';
import { type ExampleScript, ExampleScripts } from '../../examples.js';
import { CommandLanguageExtension } from '../../prism-viewer.js';
import { ConsoleModel } from '../console-model.js';

// Every example script (src/examples.ts: the stories replay them) imported into the console with the structures it
// loads, as Import or a drop does: it runs whole, with no error and no warning but the ones it lists as known. The headless viewer has no canvas, so
// the lines that need one (the camera, the background, the display) are left out here; the stories run them.

const NeedsCanvas = {
  pymol: /^\s*(orient|bg_color|zoom|center|turn|move|reset)\b/,
  vmd: /^\s*(color Display|display|axes|rotate|translate|scale)\b/,
};

// Measurement texts are drawn into a font atlas, which Node has no canvas for: a stand-in that draws nothing (as
// host.test.ts).
beforeAll(() => {
  const context = {
    measureText: () => ({ width: 8 }),
    clearRect() {},
    fillText() {},
    getImageData: (_x: number, _y: number, w: number, h: number) => ({ data: new Uint8ClampedArray(w * h * 4) }),
  };
  setCanvasModule({ createCanvas: () => ({ getContext: () => context }) });
});
// The module state this file changed, put back (files share their worker's modules: AGENTS.md).
afterAll(() => setCanvasModule(undefined as never));

const viewers: PrismViewer[] = [];
afterEach(async () => {
  for (const v of viewers.splice(0)) await v.dispose();
});

async function replay(example: ExampleScript) {
  const base = DefaultPrismViewerSpec();
  const v = await createPrismViewerHeadless({
    spec: { ...base, extensions: [...(base.extensions ?? []), CommandLanguageExtension] },
  });
  viewers.push(v);
  const console = new ConsoleModel(v.api);
  console.setDialect(example.dialect);
  const files = [
    ...example.scripts.map(
      (s) =>
        new File(
          [
            s.text
              .split('\n')
              .filter((line) => !NeedsCanvas[example.dialect].test(line))
              .join('\n'),
          ],
          s.name,
        ),
    ),
    ...example.structures.map(
      (s) => new File([readFileSync(new URL(`../../../../../../data/examples/${s.example}`, import.meta.url))], s.name),
    ),
  ];
  await console.openFiles(files);
  return { v, entries: console.entries.value };
}

describe('the example scripts, imported', () => {
  for (const example of ExampleScripts) {
    it(`${example.title}: runs whole, no error, no warning but the known ones`, async () => {
      const { entries } = await replay(example);
      // One entry: the script the others are included by.
      expect(entries.map((e) => e.input)).toEqual([example.scripts[0].name]);
      const [entry] = entries;
      const said = entry.lines.map((l) => `${l.kind}: ${l.text}`);
      expect(entry.state, said.join('\n')).toBe('done');
      expect(
        entry.lines.filter((l) => l.kind === 'error').map((l) => l.text),
        said.join('\n'),
      ).toEqual([]);
      const known = example.warns ?? [];
      const warned = entry.lines.filter((l) => l.kind === 'warning').map((l) => l.text);
      expect(
        warned.filter((w) => !known.some((k) => w.includes(k))),
        said.join('\n'),
      ).toEqual([]);
    });
  }
});
