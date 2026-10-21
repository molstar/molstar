/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import type { Model } from '@molstar/model/model/structure';
import { CustomPropertyDescriptor } from '@molstar/model/model/custom-property';
import { DefaultAssetManager, type Asset } from '@molstar/core/util/assets';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CustomModelProperty } from '../custom-model-property.js';
import type { CustomProperty } from '../custom-property.js';

/** Two models with the same topology: the second shares the first one's `customProperties`. */
async function models(): Promise<Model[]> {
  const atoms = [1, 2].map((m) => `ATOM 1 C CA ALA A 1 1 ${m} 0 0 1 10 1 ALA A CA ${m}`);
  const cif = `data_T
loop_
${'group_PDB id type_symbol label_atom_id label_comp_id label_asym_id label_entity_id label_seq_id Cartn_x Cartn_y Cartn_z occupancy B_iso_or_equiv auth_seq_id auth_comp_id auth_asym_id auth_atom_id pdbx_PDB_model_num'
  .split(' ')
  .map((c) => `_atom_site.${c}`)
  .join('\n')}
${atoms.join('\n')}
`;
  const parsed = await CIF.parseText(cif).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  return [0, 1].map((i) => trajectory.getFrameAtIndex(i) as Model);
}

function provider(type: 'static' | 'dynamic', disposed: string[]) {
  const defaultParams = { tag: PD.Text('') };
  return CustomModelProperty.createProvider({
    label: 'Test',
    descriptor: CustomPropertyDescriptor({ name: `test-assets-${type}` }),
    type,
    defaultParams,
    getParams: () => defaultParams,
    isApplicable: () => true,
    obtain: async (_ctx, model, props) => {
      const name = `${model.modelNum}:${props.tag}`;
      const asset = { dispose: () => disposed.push(name) } as unknown as Asset.Wrapper;
      return { value: name, assets: [asset] };
    },
  });
}

const ctx: CustomProperty.Context = { runtime: SyncRuntimeContext, assetManager: new DefaultAssetManager() };

describe('CustomModelProperty assets', () => {
  it('a dynamic property keeps each model its own assets', async () => {
    const [m1, m2] = await models();
    expect(m2.customProperties).toBe(m1.customProperties);

    const disposed: string[] = [];
    const p = provider('dynamic', disposed);
    await p.attach(ctx, m1);
    await p.attach(ctx, m2);
    expect(p.get(m1).value).toBe('1:');
    expect(p.get(m2).value).toBe('2:');
    expect(disposed).toEqual([]);

    // new props on one model replace that model's assets only
    await p.attach(ctx, m2, { tag: 'x' });
    expect(disposed).toEqual(['2:']);

    m1.customProperties.dispose();
    expect(disposed.sort()).toEqual(['1:', '2:', '2:x']);

    // disposed assets are forgotten: new props do not dispose them again
    await p.attach(ctx, m1, { tag: 'z' });
    expect(disposed.sort()).toEqual(['1:', '2:', '2:x']);
  });

  it('a static property replaces the assets it shares', async () => {
    const [m1, m2] = await models();
    const disposed: string[] = [];
    const p = provider('static', disposed);
    await p.attach(ctx, m1);
    await p.attach(ctx, m2, { tag: 'x' });
    expect(p.get(m2).value).toBe('2:x');
    expect(disposed).toEqual(['1:']);
  });

  it('set() with new props disposes that model’s assets only', async () => {
    const [m1, m2] = await models();
    const disposed: string[] = [];
    const p = provider('dynamic', disposed);
    await p.attach(ctx, m1);
    await p.attach(ctx, m2);
    p.set(m2, { tag: 'y' });
    expect(disposed).toEqual(['2:']);
  });
});
