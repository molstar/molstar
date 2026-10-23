/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { PolymerCenter } from '@molstar/plugin/state/helpers/periodic';
import { periodicModelProps } from '../load-helpers.js';
import { convertMvsToMolstar } from '../tree/molstar/conversion.js';
import type { MolstarNode } from '../tree/molstar/molstar-tree.js';

// A simulation's periodic box on a structure node, through `custom.molstar_periodic` (load.md §6).

afterEach(() => vi.restoreAllMocks());

describe('custom.molstar_periodic', () => {
  it('is absent: the default (molecules whole)', () => {
    expect(periodicModelProps(undefined)).toEqual({});
    expect(periodicModelProps({ other: 1 })).toEqual({});
  });

  it('takes a mode, or a mode and a center', () => {
    expect(periodicModelProps({ molstar_periodic: 'nojump' })).toEqual({
      periodic: 'nojump',
      periodicCenter: { name: 'box', params: {} },
    });
    expect(periodicModelProps({ molstar_periodic: { mode: 'box', center: 'polymers' } })).toEqual({
      periodic: 'box',
      periodicCenter: { name: 'selection', params: { script: PolymerCenter } },
    });
    expect(
      periodicModelProps({ molstar_periodic: { mode: 'box', center: { expression: 'resn ALA', language: 'pymol' } } }),
    ).toEqual({
      periodic: 'box',
      periodicCenter: { name: 'selection', params: { script: { expression: 'resn ALA', language: 'pymol' } } },
    });
  });

  it('skips an unknown mode, with a warning', () => {
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    expect(periodicModelProps({ molstar_periodic: 'sideways' })).toEqual({});
    expect(warn).toHaveBeenCalledTimes(1);
  });

  it("goes to the structure's shown model, not the topology's", () => {
    const tree = convertMvsToMolstar(
      {
        kind: 'root',
        children: [
          {
            kind: 'download',
            params: { url: 'https://files.test/a.gro' },
            children: [
              {
                kind: 'parse',
                params: { format: 'gro' },
                children: [
                  {
                    kind: 'structure',
                    params: { type: 'model', model_index: 3 },
                    custom: { molstar_periodic: 'box', other: 1 },
                  },
                ],
              },
            ],
          },
        ],
      } as never,
      undefined,
    );
    const models: MolstarNode[] = [];
    const walk = (n: MolstarNode & { children?: MolstarNode[] }) => {
      if (n.kind === 'model') models.push(n);
      for (const c of n.children ?? []) walk(c as never);
    };
    walk(tree as never);
    expect(models).toHaveLength(1);
    expect(models[0].custom).toEqual({ molstar_periodic: 'box' });
  });

  it("resolves custom.molstar_boxes's url relative to the view's, leaving the given tree as it is", () => {
    const coordinates = { kind: 'coordinates', params: {}, custom: { molstar_boxes: { url: 'protein.xsc' } } };
    const tree = {
      kind: 'root',
      children: [
        {
          kind: 'download',
          params: { url: 'protein.dcd' },
          children: [{ kind: 'parse', params: { format: 'dcd' }, children: [coordinates] }],
        },
      ],
    };
    const converted = convertMvsToMolstar(tree as never, 'https://files.test/views/view.mvsj');
    let found: unknown;
    const walk = (n: MolstarNode & { children?: MolstarNode[] }) => {
      if (n.kind === 'coordinates') found = n.custom?.molstar_boxes;
      for (const c of n.children ?? []) walk(c as never);
    };
    walk(converted as never);
    expect(found).toEqual({ url: 'https://files.test/views/protein.xsc' });
    expect(coordinates.custom.molstar_boxes.url).toBe('protein.xsc');
  });
});
