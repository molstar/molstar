/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { StructureElement } from '@molstar/model/model/structure';
import { afterEach, describe, expect, it, vi } from 'vitest';
import type { ResolvedTarget } from '../../api/targets.js';
import { DefaultPrismViewerSpec } from '../../default-spec.js';
import { createPrismViewerHeadless } from '../../headless.js';
import type { PrismViewer } from '../../viewer.js';
import { queryLanguage } from '../../targets/languages.js';
import { createTargetResolver } from '../targets.js';

// The resolver over a real plugin (targets.md §1–2).

const dataDir = path.resolve(import.meta.dirname, '../../../../../../data/examples');
const file = (name: string, as: string) => new File([fs.readFileSync(path.join(dataDir, name))], as);

let viewer: PrismViewer | undefined;
afterEach(async () => {
  await viewer?.dispose();
  viewer = undefined;
});

async function setup(...files: [string, string][]) {
  viewer = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec() });
  for (const [name, as] of files) await viewer.api.load.files({ files: [file(name, as)] as never });
  return { viewer, resolver: createTargetResolver(viewer.plugin) };
}

function ok(r: Awaited<ReturnType<ReturnType<typeof createTargetResolver>['resolve']>>): ResolvedTarget {
  if (!r.ok) throw new Error(`${r.code}: ${r.message}`);
  return r;
}

describe('the target resolver', () => {
  it('resolves keywords, addresses and query text on the only structure', async () => {
    const { resolver } = await setup(['1cbs_updated.cif', '1cbs.cif']);
    const ligand = ok(await resolver.resolve('ligand'));
    expect(ligand.scope).toBe('only');
    expect(ligand.parts).toHaveLength(1);
    expect(ligand.parts[0].selected.atoms).toBeGreaterThan(0);

    const range = ok(await resolver.resolve('A:10-20'));
    expect(range.parts[0].selected.residues).toBe(11);
    expect(range.summary).toMatch(/11 residues/);

    const query = ok(await resolver.resolve('resn REA'));
    expect(query.parts[0].selected.atoms).toBe(ligand.parts[0].selected.atoms);
    expect(query.stored).toHaveProperty('molql');
  });

  it('says what does not parse, and when nothing is loaded', async () => {
    const { viewer: v, resolver } = await setup(['1crn.cif', '1crn.cif']);
    const bad = await resolver.resolve('resn (((');
    expect(bad).toMatchObject({ ok: false, code: 'invalid' });
    await v.api.session.clear();
    expect(await resolver.resolve('ligand')).toMatchObject({ ok: false, code: 'not-found' });
  });

  it('applies to every loaded structure unless one is named', async () => {
    const { resolver } = await setup(['1crn.cif', '1crn.cif'], ['1cbs_updated.cif', '1cbs.cif']);
    const all = ok(await resolver.resolve('polymer'));
    expect(all.scope).toBe('all');
    expect(all.parts).toHaveLength(2);
    const named = ok(await resolver.resolve('polymer', { structure: '1CBS' }));
    expect(named.scope).toBe('named');
    expect(named.parts.map((p) => p.structure.label)).toHaveLength(1);
    expect(await resolver.resolve('polymer', { structure: 'nope' })).toMatchObject({ ok: false, code: 'not-found' });
  });

  it('reads and compiles a target once for every structure in scope', async () => {
    const { resolver } = await setup(['1crn.cif', '1crn.cif'], ['1cbs_updated.cif', '1cbs.cif']);
    const transpile = vi.spyOn(queryLanguage('pymol')!, 'transpile');
    const all = ok(await resolver.resolve('resn ALA and name CA'));
    expect(all.parts).toHaveLength(2);
    expect(all.parts.every((p) => p.selected.atoms > 0)).toBe(true);
    expect(transpile).toHaveBeenCalledTimes(1);
    transpile.mockRestore();
  });

  it("reads 'selection' as the person's selection on that structure", async () => {
    const { viewer: v, resolver } = await setup(['1cbs_updated.cif', '1cbs.cif']);
    expect(await resolver.resolve('selection')).toMatchObject({ ok: false });
    const range = ok(await resolver.resolve('A:30-34'));
    v.plugin.managers.structure.selection.fromLoci('set', range.parts[0].loci as StructureElement.Loci);
    const selection = ok(await resolver.resolve('selection'));
    expect(selection.scope).toBe('selected');
    expect(selection.parts[0].selected.residues).toBe(5);
  });
});
