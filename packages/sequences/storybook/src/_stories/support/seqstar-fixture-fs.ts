/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `node:fs` for the Seq* test fixtures in the browser (main.ts sends their `node:fs` import here): the format's
// load-fixtures.ts and the model's alignment-fixtures.ts read the Seq* JSON fixtures (packages/sequences/core/_test/fixtures) with
// `readFileSync(new URL(…, import.meta.url))`; here those files come bundled, by name (or content). Nothing else is served.

import brs from '../../../../core/_test/fixtures/1brs.json?raw';
import iep from '../../../../core/_test/fixtures/1iep.json?raw';
import pae from '../../../../core/_test/fixtures/af-barnase-barstar-pae.json?raw';
import p53 from '../../../../core/_test/fixtures/p53.json?raw';
import pf00042 from '../../../../core/_test/fixtures/pf00042.json?raw';

const FILES: Readonly<Record<string, string>> = {
  '1iep.json': iep,
  '1brs.json': brs,
  'af-barnase-barstar-pae.json': pae,
  'p53.json': p53,
  'pf00042.json': pf00042,
};

/** The text of a `data:` URL, which a build inlines for a small asset. */
function dataUrlText(url: string): string | undefined {
  const m = /^data:[^,]*?(;base64)?,(.*)$/s.exec(url);
  if (!m) return undefined;
  if (!m[1]) return decodeURIComponent(m[2]);
  return new TextDecoder().decode(Uint8Array.from(atob(m[2]), (c) => c.charCodeAt(0)));
}

/**
 * The fixture a `new URL('…/name.json', import.meta.url)` names. In a Storybook build the URL is the asset's, whose name
 * carries a content hash (`name-AbC123_x.json`), or a `data:` URL for a small file (the file itself).
 */
export function readFileSync(path: string | URL): string {
  const url = String(path);
  const inlined = dataUrlText(url);
  if (inlined !== undefined) return inlined;
  const name = String(path instanceof URL ? path.pathname : path)
    .split('/')
    .pop() as string;
  const text = FILES[name] ?? FILES[name.replace(/-[\w-]{8}(\.[^.]+)$/, '$1')];
  if (text === undefined) throw new Error(`no bundled fixture ${name}`);
  return text;
}
