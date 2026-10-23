/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @module-tag heavy

import { describe, expect, it } from 'vitest';
import { parseGRO } from '@molstar/io/reader/gro/parser';
import { trajectoryFromGRO } from '@molstar/model/formats/structure/gro';
import { Model } from '../model.js';
import { ModelPeriodicity } from '../periodic.js';

// The periodic box's cost by size: the unwrap plan (with its distance links) and a
// frame made whole and put in the box, at n and 4n atoms up to a million, about linear, and a million atoms in a frame's
// budget. A solvated system as a simulation writes it: mostly water, 50-residue chains and two-residue ligands, every
// atom wrapped into a cubic box at water's density (0.1 atoms per Å³).

/** A deterministic sequence in [0, 1) */
function random(seed: number) {
  let s = seed;
  return () => {
    s = (s * 1103515245 + 12345) % 2147483648;
    return s / 2147483648;
  };
}

/** A GRO frame of about `n` atoms. */
function system(n: number): string {
  const side = Math.cbrt(n / 0.1) / 10; // nm
  const rnd = random(n);
  const lines: string[] = [];
  let residue = 0;
  const atom = (res: string, name: string, x: number, y: number, z: number) => {
    const w = (v: number) => (v - Math.floor(v / side) * side).toFixed(3).padStart(8);
    lines.push(
      `${String((residue % 99999) + 1).padStart(5)}${res.padEnd(5)}${name.padStart(5)}${String((lines.length % 99999) + 1).padStart(5)}${w(x)}${w(y)}${w(z)}`,
    );
  };
  // 4% of the atoms in chains, 1% in ligands, the rest water
  const chains = Math.max(1, Math.round((0.04 * n) / 250));
  const ligands = Math.max(1, Math.round((0.01 * n) / 9));
  for (let c = 0; c < chains; c++) {
    const [x0, y0, z0] = [rnd() * side, rnd() * side, rnd() * side];
    for (let r = 0; r < 50; r++, residue++) {
      const x = x0 + 0.375 * r;
      atom('ALA', 'N', x, y0, z0);
      atom('ALA', 'CA', x + 0.121, y0 + 0.07, z0);
      atom('ALA', 'C', x + 0.242, y0, z0);
      atom('ALA', 'O', x + 0.242, y0 - 0.123, z0);
      atom('ALA', 'CB', x + 0.121, y0 + 0.14, z0 + 0.12);
    }
  }
  for (let l = 0; l < ligands; l++) {
    const [x0, y0, z0] = [rnd() * side, rnd() * side, rnd() * side];
    for (let k = 0; k < 8; k++) {
      if (k === 0 || k === 4) residue++;
      atom(k < 4 ? 'LG1' : 'LG2', `C${(k % 4) + 1}`, x0 + 0.14 * k, y0 + (k % 2 ? 0.025 : -0.025), z0);
    }
    atom('LG2', 'O1', x0 + 0.98, y0, z0);
  }
  while (lines.length < n) {
    residue++;
    const [x, y, z] = [rnd() * side, rnd() * side, rnd() * side];
    atom('SOL', 'OW', x, y, z);
    atom('SOL', 'HW1', x + 0.076, y + 0.059, z);
    atom('SOL', 'HW2', x - 0.076, y + 0.059, z);
  }
  const s = side.toFixed(5).padStart(10);
  return `Scaling system\n${String(lines.length).padStart(5)}\n${lines.join('\n')}\n${s}${s}${s}\n`;
}

async function frameOf(n: number): Promise<Model> {
  const parsed = await parseGRO(system(n)).run();
  if (parsed.isError) throw new Error(parsed.message);
  return (await trajectoryFromGRO(parsed.result).run()).representative;
}

const boxOf = (m: Model) => Model.PeriodicBox.get(m)!;

/** The best of `runs` timings of `f`, ms. */
function best(runs: number, f: () => void) {
  let min = Infinity;
  for (let i = 0; i < runs; i++) {
    const t0 = performance.now();
    f();
    min = Math.min(min, performance.now() - t0);
  }
  return min;
}

describe('the periodic box at a million atoms', () => {
  it('builds the plan and fixes a frame in about linear time, a million atoms within a frame budget', {
    timeout: 300_000,
  }, async () => {
    const small = await frameOf(250_000);
    const large = await frameOf(1_000_000);
    const cost = (m: Model) => ({
      // the plan once per topology (its distance links), then a frame made whole and in the box
      plan: best(1, () => ModelPeriodicity.getPlan(m)),
      whole: best(5, () => ModelPeriodicity.applyWhole(ModelPeriodicity.getPlan(m), boxOf(m), m.atomicConformation)),
      box: best(5, () => ModelPeriodicity.applyBox(ModelPeriodicity.getPlan(m), boxOf(m), m.atomicConformation)),
    });
    const a = cost(small);
    const b = cost(large);
    console.info(
      `periodic box: plan ${a.plan.toFixed(0)} → ${b.plan.toFixed(0)} ms, whole ${a.whole.toFixed(1)} → ${b.whole.toFixed(1)} ms, box ${a.box.toFixed(1)} → ${b.box.toFixed(1)} ms (250k → 1M atoms)`,
    );
    // about linear: 4× the atoms, well under 16×
    expect(b.plan / a.plan).toBeLessThan(6);
    expect(b.whole / a.whole).toBeLessThan(6);
    expect(b.box / a.box).toBeLessThan(6);
    // a frame of a million atoms within a playing trajectory's budget (generous: CI machines are slow)
    expect(b.whole).toBeLessThan(250);
    expect(b.box).toBeLessThan(250);
  });
});
