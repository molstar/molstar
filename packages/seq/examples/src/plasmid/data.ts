/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import { createSeqSpec } from '@molstar/seqx-core/lang';

// The Plasmid dataset, a stand-in until its dataset is recorded: a cyclic construct of pUC19's length (2,686 nt) whose
// letters are generated from a seed (never stored), and features laid out as a cloning vector's (bla on the minus
// strand, lacZα and its multiple cloning site, the origin), with restriction sites. Stand-in facts, not pUC19's: the
// recorded entry replaces the letters and derives the features from them. `rotation` moves the origin (the rotated
// variant: inside bla, so the CDS crosses it, with a plus-strand primer site across the same origin, EC-11).

/** The construct's length, pUC19's. */
export const PlasmidLength = 2686;

/** The stand-in's letters: generated from a seed. */
export function plasmidLetters(): string {
  const r = new PCG(19);
  let s = '';
  for (let i = 0; i < PlasmidLength; i++) s += 'ACGT'[r.int() & 3];
  return s;
}

/** A feature on the published origin: its start, end (exclusive), id, label and strand. */
interface Feature {
  readonly start: number;
  readonly end: number;
  readonly id: string;
  readonly label: string;
  readonly strand: '+' | '-';
  readonly category: string;
}

/** The features on the published origin (stand-in places, a cloning vector's layout). */
export const PlasmidFeatures: readonly Feature[] = [
  { start: 146, end: 470, id: 'lacZa', label: 'lacZα', strand: '-', category: 'CDS' },
  { start: 396, end: 454, id: 'mcs', label: 'MCS', strand: '+', category: 'misc' },
  { start: 1153, end: 1742, id: 'ori', label: 'ori', strand: '+', category: 'origin' },
  { start: 1625, end: 2486, id: 'bla', label: 'bla (AmpR)', strand: '-', category: 'CDS' },
  { start: 2486, end: 2591, id: 'bla-promoter', label: 'AmpR promoter', strand: '-', category: 'promoter' },
  { start: 514, end: 545, id: 'lac-promoter', label: 'lac promoter', strand: '+', category: 'promoter' },
  { start: 570, end: 587, id: 'cap', label: 'CAP site', strand: '+', category: 'misc' },
];

/** Restriction sites (stand-in places): the cloning site's enzymes and a few unique cutters. */
export const PlasmidSites: readonly { readonly at: number; readonly name: string }[] = [
  { at: 396, name: 'EcoRI' },
  { at: 402, name: 'SacI' },
  { at: 408, name: 'KpnI' },
  { at: 414, name: 'SmaI' },
  { at: 420, name: 'BamHI' },
  { at: 426, name: 'XbaI' },
  { at: 432, name: 'SalI' },
  { at: 438, name: 'PstI' },
  { at: 448, name: 'HindIII' },
  { at: 683, name: 'NdeI' },
  { at: 2617, name: 'AatII' },
  { at: 1924, name: 'ScaI' },
];

/** Where the rotated variant's origin is on the published one: inside bla. */
export const RotatedOrigin = 2000;

/** A position on the published origin, on a construct whose origin is at `origin`. */
export const rotate = (p: number, origin: number) => (((p - origin) % PlasmidLength) + PlasmidLength) % PlasmidLength;

/** The data and its handles, on a construct whose origin is at `origin` of the published one. */
export function plasmidData(b: ReturnType<typeof createSeqSpec>, origin = 0) {
  const letters = plasmidLetters();
  const turned = origin ? letters.slice(origin) + letters.slice(0, origin) : letters;
  const seq = b.data.sequence('plasmid', { alphabet: 'dna', letters: turned, cyclic: true, name: 'pUC19 (stand-in)' });
  const rows = (strand: '+' | '-') => {
    const f = PlasmidFeatures.filter((x) => x.strand === strand);
    // the rotated variant's primer site across the origin, on the plus strand (generated, EC-11)
    const extra =
      origin && strand === '+' ? [{ start: origin - 30, end: origin + 25, id: 'primer', label: 'primer site' }] : [];
    const all = [...f, ...extra];
    return {
      start: all.map((x) => rotate(x.start, origin)),
      end: all.map((x) => rotate(x.end, origin)),
      id: all.map((x) => x.id),
      label: all.map((x) => x.label),
      strand: all.map(() => strand),
    };
  };
  const plus = b.data.ranges('plus', { target: seq, rows: rows('+') });
  const minus = b.data.ranges('minus', { target: seq, rows: rows('-') });
  const sites = b.data.ranges('sites', {
    target: seq,
    rows: {
      start: PlasmidSites.map((s) => rotate(s.at, origin)),
      end: PlasmidSites.map((s) => rotate(s.at + 6, origin)),
      id: PlasmidSites.map((s) => s.name),
      label: PlasmidSites.map((s) => s.name),
    },
  });
  return { seq, plus, minus, sites };
}
