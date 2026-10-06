/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransformer } from '@molstar/core/state';
import { PluginState } from '@molstar/plugin/state';

const Valid = StateTransformer.ROOT.id;

function tree(...ids: string[]): any {
  return {
    tree: {
      transforms: ids.map((transformer, i) => ({
        parent: i === 0 ? `ref-${i}` : 'ref-0',
        transformer,
        params: void 0,
        ref: `ref-${i}`,
        version: '1',
      })),
    },
  };
}

function snapshot(parts: Partial<PluginState.Snapshot>): PluginState.Snapshot {
  return { id: 'snapshot-id' as any, ...parts };
}

const message = (ids: string) =>
  `Snapshot uses transformers that are not available in this plugin: ${ids}. Import the modules that define them.`;

describe('StateTransformer.has', () => {
  it('does not throw for unknown ids', () => {
    expect(StateTransformer.has(Valid)).toBe(true);
    expect(StateTransformer.has('test.not-registered')).toBe(false);
  });
});

describe('PluginState.validateSnapshotTransformers', () => {
  it('accepts a snapshot whose transformers are all registered', () => {
    const s = snapshot({
      behaviour: tree(Valid),
      data: tree(Valid, Valid),
      transition: { frames: [{ durationInMs: 1, data: tree(Valid) }] },
    });
    expect(() => PluginState.validateSnapshotTransformers(s)).not.toThrow();
    expect(() => PluginState.validateSnapshotTransformers(snapshot({}))).not.toThrow();
  });

  it('rejects an unknown id in the data tree', () => {
    expect(() =>
      PluginState.validateSnapshotTransformers(snapshot({ data: tree(Valid, 'test.missing-data') })),
    ).toThrow(message('test.missing-data'));
  });

  it('rejects an unknown id in the behavior tree', () => {
    expect(() =>
      PluginState.validateSnapshotTransformers(snapshot({ behaviour: tree('test.missing-behavior') })),
    ).toThrow(message('test.missing-behavior'));
  });

  it('rejects an unknown id in a transition frame', () => {
    const s = snapshot({
      data: tree(Valid),
      transition: {
        frames: [
          { durationInMs: 1, data: tree(Valid) },
          { durationInMs: 1, data: tree('test.missing-frame') },
        ],
      },
    });
    expect(() => PluginState.validateSnapshotTransformers(s)).toThrow(message('test.missing-frame'));
  });

  it('lists every missing id once in a single error', () => {
    const s = snapshot({
      behaviour: tree('test.a'),
      data: tree('test.b', 'test.a'),
      transition: { frames: [{ durationInMs: 1, data: tree('test.c') }] },
    });
    expect(() => PluginState.validateSnapshotTransformers(s)).toThrow(message('test.a, test.b, test.c'));
  });
});

describe('PluginState.setSnapshot', () => {
  it('fails before any side effect when a transformer is not registered', async () => {
    const stop = jest.fn();
    const sideEffect = jest.fn();
    const fake: any = {
      plugin: {
        managers: { animation: { stop }, structure: { component: { _setSnapshotState: sideEffect } } },
        runTask: sideEffect,
      },
      behaviors: { setSnapshot: sideEffect },
      data: { setSnapshot: sideEffect },
    };
    const s = snapshot({
      structureComponentManager: { options: {} as any },
      behaviour: tree(Valid),
      data: tree('test.missing-data'),
    });

    await expect(PluginState.prototype.setSnapshot.call(fake, s)).rejects.toThrow(message('test.missing-data'));
    expect(stop).not.toHaveBeenCalled();
    expect(sideEffect).not.toHaveBeenCalled();
  });
});
