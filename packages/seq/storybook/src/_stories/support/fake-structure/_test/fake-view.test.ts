// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import type { ChainTable } from '../chain-table.js';
import { FakeStructureViewer } from '../fake-viewer.js';
import { fakeStructureElement } from '../fake-view.js';

// The fake structure viewer's card (the story's view, in a DOM): it names the residue under the pointer, and the one the
// host hovered (no pointer here), and goes with the hover.

const res = (n: number, comp: string, observed = true) => ({ label: n, auth: String(n), ins: '', comp, observed });
const table: ChainTable = {
  entry: 'TEST',
  chains: [
    {
      label_asym_id: 'A',
      auth_asym_id: 'A',
      entity: '1',
      operators: ['1_555'],
      residues: [res(1, 'MET'), res(2, 'GLY'), res(3, 'SER', false)],
    },
  ],
};

let mounted: HTMLElement | undefined;
afterEach(() => {
  mounted?.remove();
  mounted = undefined;
});

describe('the fake structure viewer’s card', () => {
  it('names the hovered residue, whoever hovers it, and hides with the hover', () => {
    const fake = new FakeStructureViewer(table);
    mounted = fakeStructureElement(fake);
    document.body.append(mounted);
    const card = mounted.querySelector<HTMLElement>('.fv-card')!;
    expect(card).not.toBeNull();
    expect(card.style.display).not.toBe('block');
    const cells = mounted.querySelectorAll<HTMLElement>('.fv-cell');
    // the pointer: over the second residue
    cells[1].dispatchEvent(new Event('pointerover', { bubbles: true }));
    expect(card.style.display).toBe('block');
    expect(card.textContent).toBe('A 2 GLY');
    // the host's hover (no pointer here): a range, named by its ends
    fake.hover.set([{ row: 0, start: 0, end: 2 }]);
    expect(card.textContent).toBe('A 1 MET to 2 GLY');
    fake.pointerOver(null);
    expect(card.style.display).toBe('none');
  });
});

describe('the fake structure viewer’s "Within 5 Å" card', () => {
  it('reads how many residues of which chain’s neighbors a chain selection took, and goes with the next selection', () => {
    const two: ChainTable = { ...table, chains: [{ ...table.chains[0], operators: ['1_555', '2_555'] }] };
    const fake = new FakeStructureViewer(two);
    mounted = fakeStructureElement(fake);
    document.body.append(mounted);
    const note = mounted.querySelector<HTMLElement>('.fv-note')!;
    expect(note.style.display).toBe('none');
    fake.click({ row: 0, index: 0 });
    fake.selectChain();
    fake.within();
    expect(note.style.display).toBe('block');
    expect(note.textContent).toBe('2 residues within 5 Å of chain A');
    fake.clear();
    expect(note.style.display).toBe('none');
  });
});
