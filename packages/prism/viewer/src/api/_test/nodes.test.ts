/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, headlessViewer } from '../../_test/harness.js';
import { ApiResult } from '../area.js';
import { node, nodeKind, parseNodeSelector, selectNodes, selectOneNode } from '../nodes.js';

// Node selectors (nodes.md): CSS-like text over the state tree, compiled to Mol*'s StateSelection.

describe('node selectors: parsing', () => {
  it('reads kinds, attributes, states, combinators and lists', () => {
    const p = parseNodeSelector('structure[label="1CBS"] > component[type=ligand] representation:hidden, #abc');
    expect('error' in p).toBe(false);
    if ('error' in p) return;
    expect(p.alternatives).toHaveLength(2);
    expect(p.alternatives[0].parts.map((x) => x.kind)).toEqual(['structure', 'component', 'representation']);
    expect(p.alternatives[0].combinators).toEqual(['child', 'descendant']);
    expect(p.alternatives[1].parts[0].attrs).toEqual([{ name: 'ref', op: '=', value: 'abc' }]);
  });

  it('says where text does not parse', () => {
    expect(parseNodeSelector('molecule')).toMatchObject({ error: expect.stringMatching(/Unknown kind 'molecule'/) });
    expect(parseNodeSelector('structure[color=red]')).toMatchObject({ error: expect.stringMatching(/attribute/) });
    expect(parseNodeSelector('structure[label="x]')).toMatchObject({ error: expect.stringMatching(/quote/) });
    expect(parseNodeSelector('structure:shiny')).toMatchObject({ error: expect.stringMatching(/:hidden/) });
    expect(parseNodeSelector('')).toMatchObject({ error: 'An empty selector.' });
  });

  it('writes the same text from the script builder', () => {
    const q = node('structure', { label: '1CBS' })
      .child('component', { label: { contains: 'lig' } })
      .find('representation');
    expect(q.text).toBe('structure[label=1CBS] > component[label*=lig] representation');
    expect(node('structure', { label: 'Ligand A' }).text).toBe('structure[label="Ligand A"]');
    expect(node.ref('abc').hidden().text).toBe('#abc:hidden');
    expect(node('structure').or('volume').child('representation').text).toBe(
      'structure > representation, volume > representation',
    );
    expect(JSON.stringify({ node: node('model') })).toBe('{"node":"model"}');
  });
});

describe('node selectors: over a loaded structure', () => {
  it('selects by kind, label, type and state, through decorators', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const state = v.plugin.state.data;
    const sel = (text: string) => {
      const r = selectNodes(state, text);
      if (ApiResult.isFailure(r)) throw new Error(r.message);
      return r;
    };
    const structures = sel('structure');
    expect(structures).toHaveLength(1);
    expect(nodeKind(state, structures[0])).toBe('structure');
    const reps = sel('representation');
    expect(reps.length).toBeGreaterThan(1);
    // Case does not matter for text; a component's static kind is its type.
    expect(sel('component[type=LIGAND]').length).toBe(1);
    expect(sel('component[type=ligand] > representation').length).toBeGreaterThan(0);
    expect(sel('structure > component').length).toBe(sel('component').length);
    expect(sel('representation[type=cartoon]').every((c) => (c.transform.params as any).type.name === 'cartoon')).toBe(
      true,
    );
    // A list is each once, in tree order; :not and :hidden narrow.
    expect(sel('representation, component representation')).toHaveLength(reps.length);
    expect(sel('representation:not([type=cartoon])').length).toBe(
      reps.length - sel('representation[type=cartoon]').length,
    );
    // The unit cell is a representation, hidden by the default preset.
    expect(sel('representation:hidden').every((c) => c.state.isHidden)).toBe(true);
    expect(sel('representation:hidden').length + sel('representation:visible').length).toBe(reps.length);
    // The builder compiles to Mol*'s StateSelection, which keeps composing.
    expect(node('structure').compile().subtree().select(state).length).toBeGreaterThan(reps.length);
    expect(sel(`#${structures[0].transform.ref}`)).toEqual(structures);
  });

  it('picks one, or fails naming the candidates', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const state = v.plugin.state.data;
    expect(ApiResult.isFailure(selectOneNode(state, 'structure', 'structure'))).toBe(false);
    expect(selectOneNode(state, 'representation')).toMatchObject({ ok: false, code: 'ambiguous' });
    expect(selectOneNode(state, 'volume')).toMatchObject({ ok: false, code: 'not-found' });
    expect(selectOneNode(state, 'structure[')).toMatchObject({ ok: false, code: 'invalid' });
  });
});
