/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { type SequenceAlphabet, StandardComponents, standardComponent, standardLetter } from '../components.js';

const Alphabets: readonly SequenceAlphabet[] = ['protein', 'dna', 'rna'];

describe('standard components', () => {
  it('give each letter its component by alphabet: A is ALA, DA, A', () => {
    expect(Alphabets.map((a) => standardComponent(a, 'A'))).toEqual(['ALA', 'DA', 'A']);
    expect(standardComponent('protein', 'm')).toBe('MET');
    expect(standardComponent('protein', 'X')).toBe('UNK');
    expect(standardComponent('dna', 'N')).toBe('DN');
    expect(standardComponent('rna', 'U')).toBe('U');
  });

  it('have no component for an ambiguity code, the gap, or a letter not of the alphabet', () => {
    for (const l of ['B', 'Z', 'J', '-', '*']) expect(standardComponent('protein', l)).toBeUndefined();
    for (const l of ['R', 'Y', 'U', '-']) expect(standardComponent('dna', l)).toBeUndefined();
    expect(standardComponent('rna', 'T')).toBeUndefined();
    expect(standardComponent('protein', 'constructor')).toBeUndefined();
  });

  it('are one to one: a component’s letter gives the component back, in every alphabet', () => {
    for (const a of Alphabets) {
      const table = StandardComponents[a];
      expect(new Set(Object.values(table)).size).toBe(Object.keys(table).length);
      for (const [letter, component] of Object.entries(table)) {
        expect(standardLetter(a, component)).toBe(letter);
        expect(standardComponent(a, standardLetter(a, component)!)).toBe(component);
      }
    }
  });

  it('read a modified component as no standard one (its parent is the writer’s or the CCD’s)', () => {
    expect(standardLetter('protein', 'MSE')).toBeUndefined();
    expect(standardLetter('rna', 'PSU')).toBeUndefined();
    expect(standardLetter('dna', 'A')).toBeUndefined();
    expect(standardLetter('protein', 'ala')).toBe('A');
  });
});
