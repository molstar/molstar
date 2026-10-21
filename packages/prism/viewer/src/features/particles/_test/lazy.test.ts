/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { describe, expect, it } from 'vitest';
import { StateTransformer } from '@molstar/core/state';
import { ParticleActions, ParticleFormats } from '@molstar/plugin/particles-registry';
import { expectOk, headlessViewer } from '../../../_test/harness.js';
import { LazyParticles } from '../../lazy/index.js';
import { particlesFile } from '../../../api/_test/volume-fixtures.js';
import type { PrismViewer } from '../../../viewer.js';

// Particles are opt-in (`ParticlesRegistry`) and the particles feature lazy (`LazyParticles`): nothing of them is in
// the viewer until a particle file, an mmCIF holding particles or a session with particles needs them.

const loaded = (v: PrismViewer) => v.extensions.loaded.value.includes('particles');
const bytes = (text: string) => async () => new TextEncoder().encode(text);
const crambin = fs.readFileSync(
  path.resolve(import.meta.dirname, '../../../../../../../data/examples/1crn.cif'),
  'utf8',
);

describe('the particles feature, lazy', () => {
  it('is not loaded, and the plugin knows no particles, until something needs them', async () => {
    const v = await headlessViewer();
    expect(loaded(v)).toBe(false);
    expect(v.plugin.dataFormats.has('relion_star_particles')).toBe(false);
    expect(v.plugin.representation.particles.registry.list).toEqual([]);

    expectOk(await v.api.load.files({ files: [particlesFile()] }));
    expect(loaded(v)).toBe(true);
    expect(v.plugin.dataFormats.has('relion_star_particles')).toBe(true);
    expect(v.plugin.managers.particles.hierarchy.current.lists).toHaveLength(1);
  });

  it('loads for an mmCIF holding particles, not for a structure', async () => {
    const v = await headlessViewer();
    await v.extensions.ensureFor({ fileExtensions: ['cif'], files: [{ name: '1crn.cif', content: bytes(crambin) }] });
    expect(loaded(v)).toBe(false);

    const cellpack = 'data_cellpack_model\n# cellPACK recipe\n_entry.id cellpack\n';
    await v.extensions.ensureFor({ fileExtensions: ['cif'], files: [{ name: 'hiv.cif', content: bytes(cellpack) }] });
    expect(loaded(v)).toBe(true);
  });

  it('loads for a session holding particle transformers', async () => {
    const v = await headlessViewer();
    await v.extensions.ensureFor({ transformers: ['ms-plugin.parse-cif', 'ms-plugin.trajectory-from-mmcif'] });
    expect(loaded(v)).toBe(false);
    await v.extensions.ensureFor({ transformers: ['ms-plugin.particles-representation-3d'] });
    expect(loaded(v)).toBe(true);
  });

  it('loads when one of its operations is called', async () => {
    const v = await headlessViewer();
    const r = await v.api.particles.list();
    expectOk(r);
    expect(r.did).toMatch(/No particles/);
    expect(loaded(v)).toBe(true);
  });

  it('declares what ParticlesRegistry brings', () => {
    const owns = LazyParticles.stub.owns!;
    for (const id of owns.transformers!) expect(StateTransformer.get(id), id).toBeDefined();
    for (const t of ParticleActions.actions!) {
      if ('id' in t && typeof t.id === 'string' && t.id.startsWith('ms-plugin.'))
        expect(owns.transformers).toContain(t.id);
    }
    const formats = ParticleFormats.formats!;
    expect([...owns.formats!].sort()).toEqual(formats.map((f) => f.name).sort());
    const extensions = formats.flatMap((f) => [...(f.stringExtensions ?? []), ...(f.binaryExtensions ?? [])]);
    expect([...owns.fileExtensions!, ...owns.sniff!.extensions].sort()).toEqual([...new Set(extensions)].sort());
  });
});
