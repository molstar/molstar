/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { createScriptSession, EmptyEnv, type ScriptEnv } from '../env.js';
import type { Command } from '../ir.js';
import { compileScript } from '../script.js';
import { VmdCommands } from '../vmd/commands.js';
import { literal, splitTcl } from '../vmd/tcl.js';
import { VmdUnsupported } from '../vmd/unsupported.js';

// VMD's dialect (_spec §3.2, vmd.md): Tcl's word rules, the closed subset evaluated, every table example, molecules,
// reps and atomselect handles.

const env: ScriptEnv = {
  ...EmptyEnv,
  objects: [
    { name: '1cbs', molecule: 0 },
    { name: '1tqn', molecule: 1 },
  ],
  top: '1tqn',
  nextMolecule: 2,
  files: ['traj.pdb'],
};

function ok(text: string, e: ScriptEnv = env, session = createScriptSession()): readonly Command[] {
  const r = compileScript(text, 'vmd', e, session);
  const errors = r.diagnostics.filter((d) => d.severity === 'error');
  if (errors.length) expect.fail(`${text}: ${errors.map((d) => `${d.code}: ${d.message}`).join('; ')}`);
  return r.commands;
}

function errorOf(text: string, e: ScriptEnv = env) {
  const d = compileScript(text, 'vmd', e).diagnostics.find((x) => x.severity === 'error');
  if (!d) expect.fail(`${text}: no error`);
  return d;
}

describe('Tcl words', () => {
  it('splits words, keeps braces literal, and finds substitutions', () => {
    const text = 'mol selection {protein and chain A}; set s [atomselect top "resname $r"]';
    const { commands, diagnostics } = splitTcl(text);
    expect(diagnostics).toEqual([]);
    expect(commands.map((c) => c.words.map((w) => literal(w) ?? '*'))).toEqual([
      ['mol', 'selection', 'protein and chain A'],
      ['set', 's', '*'],
    ]);
    const sub = commands[1].words[2].parts[0];
    expect(sub.kind).toBe('command');
  });

  it('skips comments at a command’s start, joins continued lines', () => {
    const { commands } = splitTcl('# a comment\nmol addrep \\\n top\n  # another');
    expect(commands).toHaveLength(1);
    expect(commands[0].words.map((w) => literal(w))).toEqual(['mol', 'addrep', 'top']);
  });

  it('says what is not closed', () => {
    expect(splitTcl('mol selection {protein').diagnostics[0].message).toContain('brace');
    expect(splitTcl('puts "hello').diagnostics[0].message).toContain('quote');
    expect(splitTcl('set s [atomselect top all').diagnostics[0].message).toContain('bracket');
  });
});

describe('the VMD table', () => {
  for (const c of VmdCommands) {
    for (const example of c.examples) {
      it(`compiles ${example}`, () => {
        // The examples' `$sel` reads a handle the session has from before.
        const { session } = compileScript('set sel [atomselect top "protein"]', 'vmd', env);
        const r = compileScript(example, 'vmd', env, session);
        const errors = r.diagnostics.filter((d) => d.severity === 'error');
        if (errors.length) expect.fail(`${example}: ${errors.map((d) => d.message).join('; ')}`);
      });
    }
  }

  it('fails every command Mol* does not do, saying so', () => {
    for (const u of VmdUnsupported) {
      expect(errorOf(`${u.name} x`).code, u.name).toBe('unsupported');
    }
  });
});

describe('VMD commands to IR', () => {
  it('loads with VMD’s look, the new molecule the next top', () => {
    expect(ok('mol pdbload pdb_00004hhb')[0]).toMatchObject({ from: { entry: { db: 'pdb', id: '4hhb' } } });
    const cmds = ok('mol new 4hhb type webpdb\nmol delrep 0 top\nmol representation NewCartoon\nmol addrep top');
    expect(cmds[0]).toMatchObject({
      kind: 'load',
      name: '4hhb',
      molecule: 2,
      look: 'vmd-lines',
      from: { entry: { id: '4hhb' } },
    });
    expect(cmds[1]).toMatchObject({ kind: 'rep.remove', object: '4hhb', index: 0 });
    expect(cmds[3]).toMatchObject({ kind: 'rep.add', object: '4hhb', rep: { style: 'cartoon', selection: 'all' } });
  });

  it('takes the rep defaults set before addrep', () => {
    const cmds = ok('mol color ColorID 1\nmol selection {resname HEM}\nmol material Transparent\nmol addrep 0');
    expect(cmds[3]).toMatchObject({
      kind: 'rep.add',
      object: '1cbs',
      rep: { style: 'lines', color: { rgb: 0xff0000 }, selection: 'resname HEM', material: 'transparent' },
    });
  });

  it('reads variables and atomselect handles, and their methods', () => {
    const text = 'set r HEM\nset sel [atomselect top "resname $r"]\n$sel num\n$sel text\n$sel delete';
    const r = compileScript(text, 'vmd', env);
    expect(r.commands.map((c) => c.kind)).toEqual(['read', 'print']);
    expect(r.commands[1]).toMatchObject({ text: 'resname HEM' });
    // The session the script leaves: its variables, and no handle after `$sel delete`.
    expect(r.session.vars.get('sel')).toBe('atomselect0');
    expect(r.session.handles.size).toBe(0);
  });

  it('labels atoms and measures by molid/index', () => {
    const cmds = ok('label add Bonds 0/12 0/15\nlabel add Atoms 1/3');
    expect(cmds[0]).toMatchObject({ kind: 'measure', measure: 'distance' });
    expect((cmds[0] as { atoms: { text: string }[] }).atoms.map((a) => a.text)).toEqual(['0/12', '0/15']);
    expect(cmds[1]).toMatchObject({ kind: 'label' });
  });

  it('says what is wrong', () => {
    expect(errorOf('mol addrep 7').code).toBe('no-name');
    expect(errorOf('mol representation HBonds').code).toBe('unsupported');
    expect(errorOf('mol representation Cartun').hint).toContain('Cartoon');
    expect(errorOf('mol addrp top').hint).toContain('addrep');
    expect(errorOf('mol selection {resname HEM and and}').code).toBe('selection');
    expect(errorOf('puts $nothing').code).toBe('no-name');
    expect(errorOf('foreach x {1 2} {puts $x}').code).toBe('unsupported');
    expect(errorOf('show sticks').hint).toContain('pymol show');
  });

  it('sets the background and the projection', () => {
    expect(ok('color Display Background white')[0]).toMatchObject({ kind: 'canvas', background: 0xffffff });
    expect(ok('display projection Orthographic')[0]).toMatchObject({ kind: 'canvas', projection: 'orthographic' });
    expect(errorOf('color Display Foreground white').code).toBe('unsupported');
  });
});

describe('scripts that run scripts', () => {
  it('runs a VMD script given with this one, its commands in its place', () => {
    const scripts = new Map([['reps.tcl', 'mol representation VDW\nmol addrep top']]);
    const r = compileScript('mol delrep 0 top\nsource reps.tcl', 'vmd', { ...env, scripts });
    expect(r.diagnostics).toEqual([]);
    expect(r.commands.map((c) => c.kind)).toEqual(['rep.remove', 'rep.defaults', 'rep.add']);
    expect(r.commands[2].source).toBe('reps.tcl: mol addrep top');
  });

  it('names the file and its line in an included script’s error, at the include', () => {
    const scripts = new Map([['bad.tcl', 'puts ok\nmol addrp top']]);
    const text = 'source bad.tcl';
    const d = compileScript(text, 'vmd', { ...env, scripts }).diagnostics[0];
    expect(d.message).toContain('bad.tcl line 2');
    expect(text.slice(d.at.from, d.at.to)).toBe('source bad.tcl');
  });

  it('stops a script that runs itself', () => {
    const scripts = new Map([['loop.tcl', 'source loop.tcl']]);
    expect(compileScript('source loop.tcl', 'vmd', { ...env, scripts }).diagnostics[0].message).toContain(
      'includes itself',
    );
  });
});
