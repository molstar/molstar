/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolScriptBuilder as B } from '@molstar/query-language/language/builder';
import { fail, note, type Span, warning } from '../diagnostics.js';
import { PdbIdHint, pdbId } from '../entries.js';
import type { ColorSpec, CommandBody, LabelPart, Sel } from '../ir.js';
import type { CommandHelp, Compiled, Compiler, ParamSpec } from '../table.js';
import { PymolArgs, pymolRep } from './args.js';
import { compileScript } from '@molstar/query-language/compile';
import { fixedTo } from '../selection.js';
import { PymolSettings, settingNamed } from './settings.js';
import type { PymolCall } from './tokenizer.js';
import { unquote } from './tokenizer.js';

// PyMOL's commands (pymol.md §1, _spec §2): each its help, its signature by PyMOL's argument names, the arguments
// PyMOL takes that change nothing here, and how it compiles to the IR. `PymolUnsupported` lists PyMOL's commands that
// Mol* does not do, so they fail saying so (never a silent no-op).

export interface PymolCommand extends CommandHelp {
  readonly params: readonly ParamSpec[];
  /** Arguments PyMOL takes that change nothing here, with what to say (a warning when given). */
  readonly ignored?: Readonly<Record<string, string>>;
  compile(a: PymolArgs, c: Compiler): Compiled;
}

const p = (name: string, kind: ParamSpec['kind'], optional = true, values?: readonly string[]): ParamSpec =>
  values ? { name, kind, optional, values } : { name, kind, optional };

const sel = (name = 'selection', optional = true) => p(name, 'selection', optional);

/** PyMOL's element colors, for `util.cnc` and `color atomic` (its own table's carbon, nitrogen, …). */
const ElementColors: Readonly<Record<string, number>> = {
  N: 0x3333ff,
  O: 0xff4d4d,
  S: 0xe6c640,
  H: 0xe6e6e6,
  P: 0xff8000,
  F: 0xb3ffff,
  CL: 0x1aff1a,
  BR: 0xa62929,
  I: 0x940094,
  SE: 0xffa100,
  FE: 0xe06633,
  ZN: 0x7d80b0,
  MG: 0x8aff00,
  CA: 0x3dff00,
  NA: 0xab5cf2,
  K: 0x8f40d4,
  MN: 0x9c7ac7,
  CU: 0xc88033,
  CO: 0xf090a0,
  NI: 0x50d050,
};

const CarbonColor = 0x33ff33;

/** The atoms of an element within a selection: each part intersected with `elem X`. */
function ofElement(s: Sel, element: string): Sel {
  const test = B.struct.generator.atomGroups({
    'atom-test': B.core.rel.eq([B.acp('elementSymbol'), B.es(element)]),
  });
  return {
    ...s,
    text: `(${s.text}) and elem ${element}`,
    parts: s.parts.map((x) => ({ ...x, molql: B.struct.modifier.intersectBy({ 0: x.molql, by: test }) })),
  };
}

/** The atoms with no covalent bond within a selection (PyMOL's `nb_spheres`, `nonbonded`). */
function nonbonded(s: Sel): Sel {
  const test = B.struct.generator.atomGroups({
    'atom-test': B.core.rel.eq([B.struct.atomProperty.core.bondCount(), 0]),
  });
  return {
    ...s,
    text: `(${s.text}) and not bonded`,
    parts: s.parts.map((x) => ({ ...x, molql: B.struct.modifier.intersectBy({ 0: x.molql, by: test }) })),
  };
}

/** `util.cnc`: the non-carbons by element, as PyMOL's color commands do. */
function colorNonCarbons(s: Sel): CommandBody[] {
  return Object.entries(ElementColors).map(([e, rgb]) => ({
    kind: 'color',
    sel: ofElement(s, e),
    color: { rgb },
    optional: true,
  }));
}

/** `util.cba*`: carbons in one color, the rest by element. */
const carbonsAnd =
  (rgb: number) =>
  (a: PymolArgs): Compiled => {
    const s = a.sel('selection');
    return [
      { kind: 'require', sel: s },
      { kind: 'color', sel: ofElement(s, 'C'), color: { rgb }, optional: true },
      ...colorNonCarbons(s),
    ];
  };

const UtilCarbons: Readonly<Record<string, readonly [number, string]>> = {
  'util.cbag': [0x33ff33, 'green'],
  'util.cbac': [0x00ffff, 'cyan'],
  'util.cbam': [0xff00ff, 'magenta'],
  'util.cbay': [0xffff00, 'yellow'],
  'util.cbas': [0xff9999, 'salmon'],
  'util.cbaw': [0xe6e6e6, 'white'],
  'util.cbab': [0x8080ff, 'slate'],
  'util.cbao': [0xff8c26, 'orange'],
  'util.cbap': [0xbf00bf, 'purple'],
  'util.cbak': [0xffa6d9, 'pink'],
};

/** What `spectrum` colors by, as a theme. */
const SpectrumBy: Readonly<Record<string, string>> = {
  count: 'sequence-id',
  resi: 'sequence-id',
  b: 'uncertainty',
  chain: 'chain-id',
  segi: 'chain-id',
};

const SpectrumPalettes = ['rainbow', 'rainbow_rev', 'blue_white_red', 'red_white_blue', 'blue_red', 'red_blue'];

/** A label's PyMOL expression as text parts: `resn+resi`, `"%s-%s" % (resn, resi)`, `name`, `"text"`. */
const LabelProps: Readonly<Record<string, LabelPart>> = {
  resn: { property: 'resn' },
  resi: { property: 'resi' },
  chain: { property: 'chain' },
  name: { property: 'name' },
  elem: { property: 'elem' },
  b: { property: 'b' },
  index: { property: 'index' },
  id: { property: 'index' },
};

export function labelParts(expression: string, at: Span): LabelPart[] {
  const e = expression.trim();
  const property = (t: string): LabelPart => {
    const w = t
      .trim()
      .replace(/^str\((.*)\)$/, '$1')
      .trim();
    const q = w.match(/^(["'])(.*)\1$/);
    if (q) return { text: q[2] };
    const prop = LabelProps[w];
    if (!prop) {
      fail(
        'unsupported',
        `The label expression '${e}' uses '${w}', which Mol* cannot write.`,
        at,
        `Labels take ${Object.keys(LabelProps).join(', ')}, quoted text, + and "%s" % (…).`,
      );
    }
    return prop;
  };
  // "%s-%s" % (resn, resi)
  const formatted = e.match(/^(["'])(.*?)\1\s*%\s*\(?(.*?)\)?$/);
  if (formatted) {
    const values = formatted[3].split(',').filter((x) => x.trim());
    const pieces = formatted[2].split(/%(?:\.\d+)?[sdif]/);
    if (pieces.length - 1 !== values.length) {
      fail('argument', `The label '${e}' has ${pieces.length - 1} %-fields for ${values.length} value(s).`, at);
    }
    const out: LabelPart[] = [];
    pieces.forEach((piece, i) => {
      if (piece) out.push({ text: piece });
      if (i < values.length) out.push(property(values[i]));
    });
    return out;
  }
  // resn + resi, "text" + chain: `+` outside quotes.
  const out: LabelPart[] = [];
  let quote: string | undefined;
  let start = 0;
  for (let i = 0; i <= e.length; i++) {
    const ch = e[i];
    if (quote) {
      if (ch === quote) quote = undefined;
      continue;
    }
    if (ch === '"' || ch === "'") quote = ch;
    else if (ch === '+' || i === e.length) {
      out.push(property(e.slice(start, i)));
      start = i + 1;
    }
  }
  return out;
}

const camera =
  (action: 'frame' | 'orient' | 'center') =>
  (a: PymolArgs): Compiled => {
    const animate = a.number('animate', -1) !== 0;
    if (!a.has('selection')) return [{ kind: 'camera', action: action === 'orient' ? 'orient' : 'frame', animate }];
    const buffer = a.has('buffer') ? a.number('buffer') : undefined;
    return [
      {
        kind: 'camera',
        action,
        sel: a.sel('selection'),
        animate,
        ...(buffer !== undefined ? { buffer } : {}),
      },
    ];
  };

/** `distance`, `angle`, `dihedral`: a name first when all selections are given, one atom each. */
function measure(kind: 'distance' | 'angle' | 'dihedral', n: number) {
  return (a: PymolArgs, c: Compiler): Compiled => {
    const given = [a.has('name'), ...Array.from({ length: n }, (_, i) => a.has(`selection${i + 1}`))];
    // `distance s1, s2`: PyMOL's first positional is the name; with only the selections given, they are read as such.
    let names = Array.from({ length: n }, (_, i) => `selection${i + 1}`);
    if (given.filter(Boolean).length === n && !given[n]) names = ['name', ...names.slice(0, n - 1)];
    const atoms = names.map((x) => {
      if (!a.has(x)) fail('argument', `${a.call.name} needs ${n} selections, one atom each.`, a.call.at);
      return a.sel(x);
    });
    if (kind === 'distance' && a.has('mode') && a.int('mode') !== 0) {
      fail(
        'unsupported-argument',
        `distance mode ${a.text('mode')} (contacts between selections) is not supported.`,
        a.at('mode'),
        'Give one atom on each side (mode 0); Mol* shows interactions from the Structure panel.',
      );
    }
    // With a cutoff, PyMOL measures every pair of the two sets within it.
    if (kind === 'distance' && a.has('cutoff')) {
      return [{ kind: 'measure', measure: kind, atoms, cutoff: a.number('cutoff') }];
    }
    return [{ kind: 'measure', measure: kind, atoms }];
  };
}

const superpose =
  (method: 'sequence' | 'structure') =>
  (a: PymolArgs): Compiled => [{ kind: 'superpose', method, mobile: a.sel('mobile'), target: a.sel('target') }];

const Ignored = {
  state: 'Mol* shows the frame the trajectory is at.',
  quiet: 'Mol* says what it did.',
};

export const PymolCommands: readonly PymolCommand[] = [
  {
    name: 'fetch',
    summary: 'Load an entry from the PDB (1cbs, pdb_00001cbs) or AlphaFold DB (AF-…).',
    usage: 'fetch code [, name [, type=cif]]',
    params: [p('code', 'text', false), p('name', 'name'), p('type', 'enum', true, ['cif', 'pdb', 'mmtf'])],
    ignored: {
      ...Ignored,
      async: 'Loads are whole.',
      zoom: 'The camera frames what loads.',
      path: 'Nothing is saved.',
    },
    examples: ['fetch 1cbs', 'fetch 1tqn, heme', 'fetch pdb_00001cbs', 'fetch AF-P69905-F1'],
    compile(a, c) {
      const codes = a
        .text('code')
        .split(/[\s,]+/)
        .filter(Boolean);
      return codes.map((code, i): CommandBody => {
        const af = /^AF-/i.test(code);
        const pdb = af ? undefined : pdbId(code);
        if (!af && !pdb) {
          fail('argument', `'${code}' is not a PDB id or an AlphaFold DB id.`, a.at('code'), PdbIdHint);
        }
        const wish = i === 0 && a.has('name') ? a.text('name') : code.toLowerCase();
        // An object of that name already: PyMOL fetches into it; here the loaded one stays (no second copy).
        const known = c.state.object(wish);
        return {
          kind: 'load',
          from: { entry: { db: af ? 'afdb' : 'pdb', id: af ? code.toUpperCase() : (pdb ?? code) } },
          name: known?.name ?? c.state.freeName(wish),
          molecule: known?.molecule ?? c.state.nextMolecule + i,
          copy: !known,
          look: 'default',
        };
      });
    },
  },
  {
    name: 'load',
    summary: 'Load a file given with the script, or a URL.',
    usage: 'load filename [, object [, format=…]]',
    params: [p('filename', 'file', false), p('object', 'name'), p('state', 'int'), p('format', 'text')],
    ignored: { ...Ignored, finish: 'Loads are whole.', discrete: 'Models load as a trajectory.' },
    examples: ['load https://files.rcsb.org/download/1cbs.cif', 'load model.pdb, mine'],
    compile(a, c) {
      const file = a.text('filename');
      const url = /^(https?|ftp):\/\//i.test(file);
      if (!url && !c.state.files.some((f) => f === file || f.split('/').pop() === file.split('/').pop())) {
        fail(
          'no-name',
          `There is no file '${file}' with this script.`,
          a.at('filename'),
          c.state.files.length
            ? `Files given: ${c.state.files.join(', ')}.`
            : 'Open the file with the script (drop both), or give a URL.',
        );
      }
      if (a.has('state'))
        c.report(warning('unsupported-argument', 'state is ignored: every model loads.', a.at('state')));
      const stem = (file.split(/[\\/]/).pop() ?? file).replace(/\.(gz|zip)$/i, '').replace(/\.[^.]+$/, '');
      const wish = a.has('object') ? a.text('object') : stem;
      const name = c.state.freeName(wish);
      return [
        {
          kind: 'load',
          from: url ? { url: file } : { file },
          name,
          molecule: c.state.nextMolecule,
          look: 'default',
          copy: true,
          ...(a.has('format') ? { format: a.text('format') } : {}),
        },
      ];
    },
  },
  {
    name: 'delete',
    summary: 'Remove an object (a structure), or forget a named selection; `all` removes everything.',
    usage: 'delete name',
    params: [p('name', 'object', false)],
    examples: ['delete 1cbs', 'delete all'],
    compile(a, c) {
      const name = a.text('name');
      if (name === 'all')
        return [
          { kind: 'object.remove', objects: 'all' },
          { kind: 'name.remove', names: 'all' },
        ];
      const named = c.state.name(name);
      if (named) return [{ kind: 'name.remove', names: [named.name] }];
      const object = c.state.object(name);
      if (object) return [{ kind: 'object.remove', objects: [object.name] }];
      return fail('no-name', `No object or selection is named '${name}'.`, a.at('name'), namesHint(c));
    },
  },
  ...(['enable', 'disable'] as const).map(
    (verb): PymolCommand => ({
      name: verb,
      summary: verb === 'enable' ? 'Show an object.' : 'Hide an object.',
      usage: `${verb} name`,
      params: [p('name', 'object')],
      examples: [`${verb} 1cbs`],
      compile(a, c) {
        const name = a.text('name', 'all');
        if (name === 'all') return [{ kind: 'object.visibility', objects: 'all', visible: verb === 'enable' }];
        const object = c.state.object(name);
        if (object) return [{ kind: 'object.visibility', objects: [object.name], visible: verb === 'enable' }];
        const named = c.state.name(name);
        if (named)
          return verb === 'enable' ? [{ kind: 'select', mode: 'set', sel: named.sel }] : [{ kind: 'select.clear' }];
        return fail('no-name', `No object or selection is named '${name}'.`, a.at('name'), namesHint(c));
      },
    }),
  ),
  {
    name: 'select',
    summary: 'Name a selection and select it (one argument: the name is `sele`).',
    usage: 'select [name,] selection [, enable=0]',
    params: [p('name', 'name', false), sel('selection'), p('enable', 'bool')],
    ignored: { merge: 'The selection replaces the name.', state: Ignored.state, quiet: Ignored.quiet },
    examples: ['select resn HEM', 'select lig, resn HEM', 'select near, lig around 5, enable=0'],
    compile(a, c) {
      // `select expr`: PyMOL names it `sele`.
      const one = !a.has('selection');
      const name = one ? 'sele' : a.text('name');
      if (!one && !/^[A-Za-z_][\w]*$/.test(name)) {
        fail('argument', `'${name}' is not a name.`, a.at('name'), 'Names are letters, digits and _ (lig, near_5).');
      }
      const s = one ? a.sel('name') : a.sel('selection');
      // The name is the atoms selected now: objects loaded later are not in it (PyMOL's meaning).
      const fixed = fixedTo(
        s,
        c.state.objects.map((o) => o.name),
      );
      // A name the selection language reads as its own word (`b` is the B-factor) cannot be used after: said now.
      if (!one && !usableName(name)) {
        fail(
          'argument',
          `'${name}' cannot name a selection: PyMOL's selection words take it.`,
          a.at('name'),
          'Pick another name (lig, site_1).',
        );
      }
      const out: CommandBody[] = [{ kind: 'name.set', name, sel: fixed }];
      if (a.bool('enable', true)) out.push({ kind: 'select', mode: 'set', sel: s });
      return out;
    },
  },
  {
    name: 'deselect',
    summary: 'Clear the selection.',
    usage: 'deselect',
    params: [],
    examples: ['deselect'],
    compile: () => [{ kind: 'select.clear' }],
  },
  ...(['show', 'hide'] as const).map(
    (verb): PymolCommand => ({
      name: verb,
      summary:
        verb === 'show'
          ? 'Show atoms as a representation (no arguments: lines).'
          : 'Hide a representation of atoms (no arguments: everything).',
      usage: `${verb} [representation [, selection]]`,
      params: [p('representation', 'rep'), sel()],
      examples: [
        `${verb} cartoon`,
        `${verb} sticks, resn HEM`,
        verb === 'show' ? 'show spheres, solvent' : 'hide everything, solvent',
      ],
      compile(a) {
        const rep = a.text('representation', verb === 'show' ? 'lines' : 'everything');
        const reps = pymolRep(rep, a.at('representation'));
        if (verb === 'show' && reps === 'all') {
          fail(
            'unsupported',
            'show everything draws every representation at once: not supported.',
            a.at('representation'),
            'Name one: show sticks, …',
          );
        }
        const restricted = rep === 'nb_spheres' || rep === 'nonbonded';
        const s = restricted ? nonbonded(a.sel('selection')) : a.sel('selection');
        return [{ kind: 'display', mode: verb, reps, sel: s }];
      },
    }),
  ),
  {
    name: 'as',
    summary: 'Show atoms as only this representation.',
    usage: 'as representation [, selection]',
    params: [p('representation', 'rep', false), sel()],
    examples: ['as cartoon', 'as sticks, resn HEM'],
    compile(a) {
      const reps = pymolRep(a.text('representation'), a.at('representation'));
      if (reps === 'all')
        fail('unsupported', 'as everything is not supported.', a.at('representation'), 'Name one: as cartoon.');
      return [{ kind: 'display', mode: 'only', reps, sel: a.sel('selection') }];
    },
  },
  {
    name: 'cartoon',
    summary: 'Set how a cartoon is drawn: automatic, putty or tube.',
    usage: 'cartoon type [, selection]',
    params: [p('type', 'enum', false, ['automatic', 'putty', 'tube']), sel()],
    examples: ['cartoon putty', 'cartoon automatic, chain A'],
    compile(a) {
      const type = a.text('type').toLowerCase();
      if (type !== 'automatic' && type !== 'putty' && type !== 'tube') {
        const known = ['loop', 'rect', 'oval', 'arrow', 'dumbbell', 'skip', 'cylinder'];
        if (known.includes(type)) {
          fail(
            'unsupported',
            `The cartoon type '${type}' is not supported.`,
            a.at('type'),
            'Use automatic, putty or tube.',
          );
        }
        fail('argument', `'${type}' is not a cartoon type.`, a.at('type'), 'Use automatic, putty or tube.');
      }
      return [{ kind: 'cartoon', style: type, sel: a.sel('selection') }];
    },
  },
  {
    name: 'color',
    summary: 'Color atoms: a PyMOL color name, 0xRRGGBB or [r, g, b]; `atomic` colors by element.',
    usage: 'color color [, selection]',
    params: [p('color', 'color', false), sel()],
    examples: ['color red', 'color slate, chain A', 'color 0xff8800, resn HEM', 'color atomic, resn HEM'],
    compile(a) {
      const name = a.text('color').toLowerCase();
      const s = a.sel('selection');
      if (name === 'atomic') {
        return [
          { kind: 'require', sel: s },
          { kind: 'color', sel: ofElement(s, 'C'), color: { rgb: CarbonColor }, optional: true },
          ...colorNonCarbons(s),
        ];
      }
      if (name === 'chain' || name === 'by_chain') return [{ kind: 'color', sel: s, color: { theme: 'chain-id' } }];
      return [{ kind: 'color', sel: s, color: { rgb: a.color('color') } }];
    },
  },
  ...Object.entries(UtilCarbons).map(
    ([name, [rgb, word]]): PymolCommand => ({
      name,
      summary: `Color carbons ${word} and the other atoms by element.`,
      usage: `${name} [selection]`,
      params: [sel()],
      examples: [`${name} resn HEM`],
      compile: carbonsAnd(rgb),
    }),
  ),
  {
    name: 'util.cnc',
    summary: 'Color the atoms other than carbon by element.',
    usage: 'util.cnc [selection]',
    params: [sel()],
    examples: ['util.cnc', 'util.cnc resn HEM'],
    compile: (a) => {
      const s = a.sel('selection');
      return [{ kind: 'require', sel: s }, ...colorNonCarbons(s)];
    },
  },
  ...(
    [
      ['util.cbc', 'Color by chain.', 'chain-id'],
      ['util.chainbow', 'Color each chain as a rainbow.', 'sequence-id'],
      ['util.rainbow', 'Color as a rainbow from the N to the C terminus.', 'sequence-id'],
      ['util.ss', 'Color by secondary structure (Mol*’s colors).', 'secondary-structure'],
    ] as const
  ).map(
    ([name, summary, theme]): PymolCommand => ({
      name,
      summary,
      usage: `${name} [selection]`,
      params: [sel()],
      examples: [`${name}`],
      compile(a, c) {
        if (name === 'util.ss') {
          c.report(note('Mol* colors helices, sheets and loops with its own colors, not PyMOL’s.', a.call.nameAt));
        }
        return [{ kind: 'color', sel: a.sel('selection'), color: { theme } }];
      },
    }),
  ),
  {
    name: 'spectrum',
    summary: 'Color along a property: count or resi (residue number), b (B-factor), chain.',
    usage: 'spectrum [expression [, palette [, selection]]]',
    params: [
      p('expression', 'enum', true, Object.keys(SpectrumBy)),
      p('palette', 'enum', true, SpectrumPalettes),
      sel(),
    ],
    ignored: {
      minimum: 'The range is the data’s.',
      maximum: 'The range is the data’s.',
      byres: 'Residues are colored whole.',
      quiet: Ignored.quiet,
    },
    examples: ['spectrum', 'spectrum b, blue_white_red', 'spectrum count, rainbow, chain A'],
    compile(a) {
      const by = a.text('expression', 'count').toLowerCase();
      const theme = SpectrumBy[by];
      if (!theme) {
        fail(
          'unsupported',
          `spectrum by '${by}' is not supported.`,
          a.at('expression'),
          `Use ${Object.keys(SpectrumBy).join(', ')}.`,
        );
      }
      const palette = a.text('palette', 'rainbow').toLowerCase();
      if (!SpectrumPalettes.includes(palette)) {
        fail(
          'unsupported',
          `The palette '${palette}' is not supported.`,
          a.at('palette'),
          `Palettes: ${SpectrumPalettes.join(', ')}.`,
        );
      }
      const color: ColorSpec =
        palette === 'rainbow' && theme !== 'uncertainty' ? { theme } : { theme, params: { palette } };
      return [{ kind: 'color', sel: a.sel('selection'), color }];
    },
  },
  {
    name: 'bg_color',
    summary: 'Set the background color.',
    usage: 'bg_color color',
    params: [p('color', 'color', false)],
    examples: ['bg_color white', 'bg_color black'],
    compile: (a) => [{ kind: 'canvas', background: a.color('color') }],
  },
  {
    name: 'set',
    summary: 'Change a setting (the ones Mol* has: see `help set`).',
    usage: 'set name, value [, selection]',
    params: [p('name', 'setting', false, Object.keys(PymolSettings)), p('value', 'value'), sel()],
    ignored: { state: Ignored.state, updates: 'Mol* draws again by itself.', quiet: Ignored.quiet },
    examples: [
      'set bg_rgb, white',
      'set cartoon_transparency, 0.5',
      'set sphere_scale, 0.3, resn HOH',
      'set orthoscopic, 1',
    ],
    compile(a, c) {
      const s = settingNamed(a.text('name'), a.at('name'));
      return s.set(a.has('value') ? a.text('value') : '1', a.at('value'), a, c);
    },
  },
  {
    name: 'unset',
    summary: 'Put a setting back to its default.',
    usage: 'unset name [, selection]',
    params: [p('name', 'setting', false, Object.keys(PymolSettings)), sel()],
    examples: ['unset cartoon_transparency', 'unset orthoscopic'],
    compile(a, c) {
      const s = settingNamed(a.text('name'), a.at('name'));
      return s.set(s.default, a.at('name'), a, c);
    },
  },
  {
    name: 'zoom',
    summary: 'Frame a selection (no selection: everything).',
    usage: 'zoom [selection [, buffer [, state [, complete [, animate]]]]]',
    params: [sel(), p('buffer', 'number'), p('state', 'int'), p('complete', 'bool'), p('animate', 'number')],
    examples: ['zoom', 'zoom resn HEM', 'zoom chain A, 5'],
    compile: camera('frame'),
  },
  {
    name: 'orient',
    summary: 'Orient along the principal axes (with a selection: its structure’s, then framed on it).',
    usage: 'orient [selection [, state [, animate]]]',
    params: [sel(), p('state', 'int'), p('animate', 'number')],
    examples: ['orient', 'orient resn HEM'],
    compile(a, c) {
      if (a.has('selection')) {
        c.report(note('Mol* orients along the structure’s axes, then frames the selection.', a.at('selection')));
      }
      return camera('orient')(a);
    },
  },
  {
    name: 'center',
    summary: 'Center the view on a selection, keeping the zoom.',
    usage: 'center [selection [, state [, origin [, animate]]]]',
    params: [sel(), p('state', 'int'), p('origin', 'bool'), p('animate', 'number')],
    examples: ['center', 'center resn HEM'],
    compile: camera('center'),
  },
  {
    name: 'reset',
    summary: 'Reset the camera.',
    usage: 'reset',
    params: [],
    examples: ['reset'],
    compile: () => [{ kind: 'camera', action: 'reset', animate: true }],
  },
  {
    name: 'turn',
    summary: 'Rotate the camera about an axis of the screen, in degrees.',
    usage: 'turn axis, angle',
    params: [p('axis', 'enum', false, ['x', 'y', 'z']), p('angle', 'number', false)],
    examples: ['turn y, 90', 'turn x, -45'],
    compile(a) {
      const axis = readAxis(a.text('axis'), a.at('axis'));
      return [{ kind: 'camera.turn', axis, angle: a.number('angle') }];
    },
  },
  {
    name: 'move',
    summary: 'Move the camera along an axis of the screen, in Ångströms.',
    usage: 'move axis, distance',
    params: [p('axis', 'enum', false, ['x', 'y', 'z']), p('distance', 'number', false)],
    examples: ['move z, -10', 'move x, 5'],
    compile(a) {
      const axis = readAxis(a.text('axis'), a.at('axis'));
      return [{ kind: 'camera.move', axis, distance: a.number('distance') }];
    },
  },
  {
    name: 'distance',
    aliases: ['dist'],
    summary: 'Measure the distance between two atoms, or every pair of two sets within a cutoff.',
    usage: 'distance [name,] selection1, selection2 [, cutoff [, mode]]',
    params: [p('name', 'name'), sel('selection1'), sel('selection2'), p('cutoff', 'number'), p('mode', 'int')],
    ignored: {
      label: 'Mol* labels every measurement.',
      reset: 'Measurements add up.',
      state: Ignored.state,
      quiet: Ignored.quiet,
    },
    examples: ['distance /1cbs//A/45/CA, /1cbs//A/60/CA', 'distance d1, A/45/OG, A/60/N'],
    compile: measure('distance', 2),
  },
  {
    name: 'angle',
    summary: 'Measure the angle of three atoms.',
    usage: 'angle [name,] selection1, selection2, selection3',
    params: [p('name', 'name'), sel('selection1'), sel('selection2'), sel('selection3')],
    ignored: {
      mode: 'Three atoms are measured.',
      label: 'Mol* labels every measurement.',
      state: Ignored.state,
      quiet: Ignored.quiet,
    },
    examples: ['angle A/45/N, A/45/CA, A/45/C'],
    compile: measure('angle', 3),
  },
  {
    name: 'dihedral',
    summary: 'Measure the dihedral angle of four atoms.',
    usage: 'dihedral [name,] selection1, selection2, selection3, selection4',
    params: [p('name', 'name'), sel('selection1'), sel('selection2'), sel('selection3'), sel('selection4')],
    ignored: {
      mode: 'Four atoms are measured.',
      label: 'Mol* labels every measurement.',
      state: Ignored.state,
      quiet: Ignored.quiet,
    },
    examples: ['dihedral A/45/N, A/45/CA, A/45/C, A/46/N'],
    compile: measure('dihedral', 4),
  },
  {
    name: 'label',
    summary: 'Label atoms: resn, resi, chain, name, elem, b, quoted text, + and "%s" % (…); empty clears.',
    usage: 'label [selection [, expression]]',
    params: [sel(), p('expression', 'text')],
    examples: ['label n. CA and resi 45, resn+resi', 'label resn HEM and name FE, "heme iron"', 'label all, ""'],
    compile(a) {
      const s = a.sel('selection');
      const e = a.has('expression') ? rawText(a, 'expression') : '';
      if (!e || /^(["'])\1$/.test(e.trim())) return [{ kind: 'label.clear', sel: s }];
      return [{ kind: 'label', sel: s, text: labelParts(e, a.at('expression')) }];
    },
  },
  {
    name: 'align',
    summary: 'Superpose by sequence alignment (Mol*’s chain superposition); prints the RMSD.',
    usage: 'align mobile, target',
    params: [sel('mobile', false), sel('target', false)],
    ignored: {
      cutoff: 'Mol* does not reject outliers.',
      cycles: 'Mol* does not reject outliers.',
      gap: 'Mol*’s alignment scores apply.',
      extend: 'Mol*’s alignment scores apply.',
      max_gap: 'Mol*’s alignment scores apply.',
      object: 'No alignment object is made.',
      matrix: 'Mol*’s alignment scores apply.',
      transform: 'The mobile structure is always moved.',
      quiet: Ignored.quiet,
    },
    examples: ['align 1tqn, 1cbs', 'align 1tqn and chain A, 1cbs and chain A'],
    compile: superpose('sequence'),
  },
  {
    name: 'super',
    summary: 'Superpose by structure (TM-align); prints the RMSD and TM-score.',
    usage: 'super mobile, target',
    params: [sel('mobile', false), sel('target', false)],
    ignored: {
      cutoff: 'TM-align applies.',
      cycles: 'TM-align applies.',
      object: 'No alignment object is made.',
      quiet: Ignored.quiet,
    },
    examples: ['super 1tqn, 1cbs'],
    compile: superpose('structure'),
  },
  {
    name: 'cealign',
    summary: 'Superpose by structure (TM-align, Mol*’s structure alignment); prints the RMSD.',
    usage: 'cealign target, mobile',
    params: [sel('target', false), sel('mobile', false)],
    ignored: {
      d0: 'TM-align applies.',
      d1: 'TM-align applies.',
      window: 'TM-align applies.',
      gap_max: 'TM-align applies.',
      quiet: Ignored.quiet,
    },
    examples: ['cealign 1cbs, 1tqn'],
    compile: (a) => [{ kind: 'superpose', method: 'structure', mobile: a.sel('mobile'), target: a.sel('target') }],
  },
  {
    name: 'pair_fit',
    summary: 'Superpose by pairs of atoms (mobile, target, mobile, target, …).',
    usage: 'pair_fit mobile1, target1 [, mobile2, target2, …]',
    params: Array.from({ length: 16 }, (_, i) =>
      sel(`${i % 2 === 0 ? 'mobile' : 'target'}${Math.floor(i / 2) + 1}`, i > 1),
    ),
    examples: ['pair_fit 1tqn and name CA and resi 10-50, 1cbs and name CA and resi 10-50'],
    compile(a) {
      const mobile: Sel[] = [];
      const target: Sel[] = [];
      for (let i = 1; i <= 8; i++) {
        if (!a.has(`mobile${i}`)) break;
        if (!a.has(`target${i}`)) fail('argument', `pair_fit needs a target for mobile${i}.`, a.at(`mobile${i}`));
        mobile.push(a.sel(`mobile${i}`));
        target.push(a.sel(`target${i}`));
      }
      return [
        {
          kind: 'superpose',
          method: 'atoms',
          mobile: joinSels(mobile),
          target: joinSels(target),
          pairs: mobile.map((m, i) => ({ mobile: m, target: target[i] })),
        },
      ];
    },
  },
  {
    name: 'remove',
    summary: 'Leave whole entities or chains out of the scene (hydrogens: `remove hydro`).',
    usage: 'remove selection',
    params: [sel('selection', false)],
    examples: ['remove solvent', 'remove hydro', 'remove chain B'],
    compile(a) {
      const t = a.text('selection').trim().toLowerCase();
      if (t === 'hydro' || t === 'h.' || t === 'hydrogens') return [{ kind: 'hydrogens', visible: false }];
      return [{ kind: 'hide-contents', sel: a.sel('selection') }];
    },
  },
  {
    name: 'png',
    summary: 'Save an image (PNG).',
    usage: 'png filename [, width [, height [, dpi [, ray]]]]',
    params: [p('filename', 'file'), p('width', 'text'), p('height', 'text'), p('dpi', 'number'), p('ray', 'bool')],
    ignored: { quiet: Ignored.quiet, prior: 'The image is drawn now.', format: 'Images are PNG.' },
    examples: ['png image.png', 'png image.png, 1200, 900'],
    compile(a, c) {
      if (a.has('dpi')) c.report(note('dpi is ignored: the image has the pixels asked for.', a.at('dpi')));
      if (a.has('ray')) c.report(note('ray is ignored: Mol* draws the same image without ray tracing.', a.at('ray')));
      const size = (name: string) => {
        if (!a.has(name)) return undefined;
        const t = a.text(name);
        const v = Number(t.replace(/px$/, ''));
        if (!Number.isFinite(v) || v <= 0) fail('argument', `${name} is a size in pixels, not '${t}'.`, a.at(name));
        return Math.round(v);
      };
      const width = size('width');
      const height = size('height');
      const name = a.text('filename').replace(/\.png$/i, '');
      return [
        {
          kind: 'export.image',
          ...(name ? { name } : {}),
          ...(width ? { width } : {}),
          ...(height ? { height } : {}),
          ...(c.state.session.vars.get('ray_opaque_background') === '0' ? { transparent: true } : {}),
        },
      ];
    },
  },
  ...(['ray', 'draw'] as const).map(
    (name): PymolCommand => ({
      name,
      summary: 'Nothing to do: Mol* renders live (`png` saves the image).',
      usage: `${name} [width [, height]]`,
      params: [p('width', 'number'), p('height', 'number')],
      examples: [name],
      compile(a, c) {
        c.report(note(`${name} is not needed: Mol* renders live; png saves the image.`, a.call.nameAt));
        return [];
      },
    }),
  ),
  {
    name: 'count_atoms',
    summary: 'Count the atoms of a selection.',
    usage: 'count_atoms [selection]',
    params: [sel()],
    ignored: { quiet: Ignored.quiet, state: Ignored.state },
    examples: ['count_atoms', 'count_atoms resn HEM'],
    compile: (a) => [{ kind: 'read', read: { kind: 'count', sel: a.sel('selection') } }],
  },
  {
    name: 'get_names',
    summary: 'List the objects and the named selections.',
    usage: 'get_names',
    params: [p('type', 'enum', true, ['all', 'objects', 'selections'])],
    examples: ['get_names'],
    compile: () => [{ kind: 'read', read: { kind: 'names' } }],
  },
  {
    name: 'get_object_list',
    summary: 'List the objects.',
    usage: 'get_object_list',
    params: [sel()],
    examples: ['get_object_list'],
    compile: () => [{ kind: 'read', read: { kind: 'objects' } }],
  },
  {
    name: 'get_view',
    summary: 'Describe the camera.',
    usage: 'get_view',
    params: [],
    examples: ['get_view'],
    compile: () => [{ kind: 'read', read: { kind: 'view' } }],
  },
  {
    name: 'frame',
    summary: 'Show a frame of the trajectory (1-based, as PyMOL’s states).',
    usage: 'frame number',
    params: [p('number', 'int', false)],
    examples: ['frame 1', 'frame 10'],
    compile(a) {
      const n = a.int('number');
      if (n < 1) fail('argument', 'Frames are counted from 1.', a.at('number'));
      return [{ kind: 'trajectory', frame: n - 1 }];
    },
  },
  {
    name: 'mplay',
    summary: 'Play the trajectory.',
    usage: 'mplay',
    params: [],
    examples: ['mplay'],
    compile: () => [{ kind: 'trajectory', play: 'forward' }],
  },
  {
    name: 'mstop',
    summary: 'Stop playing the trajectory.',
    usage: 'mstop',
    params: [],
    examples: ['mstop'],
    compile: () => [{ kind: 'trajectory', play: 'stop' }],
  },
  ...(['undo', 'redo'] as const).map(
    (action): PymolCommand => ({
      name: action,
      summary: action === 'undo' ? 'Undo the last step (alone on its line).' : 'Redo what Undo took back.',
      usage: action,
      params: [],
      examples: [action],
      compile: () => [{ kind: 'history', action }],
    }),
  ),
];

function readAxis(text: string, at: Span): 'x' | 'y' | 'z' {
  const t = text.trim().toLowerCase();
  if (t === 'x' || t === 'y' || t === 'z') return t;
  return fail('argument', `'${text}' is not an axis.`, at, 'Give x, y or z.');
}

/** The selections as one: their parts merged per scope. */
function joinSels(xs: readonly Sel[]): Sel {
  if (xs.length === 1) return xs[0];
  const scopes = new Map<string, { scope: Sel['parts'][number]['scope']; molqls: Sel['parts'][number]['molql'][] }>();
  for (const s of xs) {
    for (const part of s.parts) {
      const key = part.scope === 'others' ? '' : part.scope.object.toLowerCase();
      const e = scopes.get(key) ?? { scope: part.scope, molqls: [] };
      e.molqls.push(part.molql);
      scopes.set(key, e);
    }
  }
  return {
    text: xs.map((x) => `(${x.text})`).join(' or '),
    dialect: xs[0].dialect,
    parts: [...scopes.values()].map((e) => ({ scope: e.scope, molql: B.struct.combinator.merge(e.molqls) })),
  };
}

/** Whether a name reads as itself in a selection (not as a keyword or a property of PyMOL's). */
function usableName(name: string): boolean {
  // A marker no keyword compiles to: the name reads as itself when the text is the marker.
  const marker = B.struct.generator.atomGroups({ 'atom-test': B.core.rel.eq([B.acp('sourceIndex'), -1]) });
  try {
    const probe = compileScript('pymol', name, { resolveName: (n) => (n === name ? marker : undefined) });
    return JSON.stringify(probe) === JSON.stringify(marker);
  } catch {
    return false;
  }
}

/** The raw text of an argument (quotes kept: a label's expression reads them). */
function rawText(a: PymolArgs, name: string): string {
  const at = a.at(name);
  const arg = a.call.args.find((x) => x.at.from === at.from && x.at.to === at.to);
  return arg?.text ?? '';
}

function namesHint(c: Compiler): string {
  const names = [...c.state.objects.map((o) => o.name), ...c.state.names.keys()];
  return names.length ? `Names: ${names.join(', ')}.` : 'Nothing is loaded or named yet.';
}

/** Every PyMOL command by name and alias. */
export const PymolCommandsByName: ReadonlyMap<string, PymolCommand> = new Map(
  PymolCommands.flatMap((c) => [c.name, ...(c.aliases ?? [])].map((n) => [n, c] as const)),
);

/** Compiles one PyMOL command call. */
export function compilePymolCall(call: PymolCall, command: PymolCommand, c: Compiler): Compiled {
  const args = new PymolArgs(call, command.params, c, command.ignored);
  return command.compile(args, c);
}

/** `unquote`, for the table's users. */
export { unquote };
