/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Color schemes docs page's live parts: every registered scheme from the kernel's table itself (label, what it
// fits, colorblind-safe, provenance) with its light and dark swatches through `schemeTheme`, so the page never drifts
// from the code; and the page's code snippet (kept here so its indentation survives MDX).

import {
  Color,
  type ColorSchemeEntry,
  createRegistries,
  KernelColorSchemes,
  KernelRegistry,
  register,
  schemeTheme,
} from '@molstar/seq-core/format';
import f from '../support/seqstar-format.module.css';
import cs from '../support/seqstar-schemes.module.css';

const reg = createRegistries();
register(reg, KernelRegistry);

/** A few colors that show a scheme: its letters' (one each), its ramp sampled, its categories, a table's rules or legend. */
function swatches(e: ColorSchemeEntry, dark: boolean): string[] {
  const named = schemeTheme(reg, { name: e.name }, { dark });
  // A named rule table (Clustal X, Percent identity) as its rules, as resolution expands it.
  const expand = named ? reg.get('themes', named.theme)?.expand : undefined;
  const t = named && expand ? expand(named.params ?? {}, 'color') : named;
  const p = (t?.params ?? {}) as Record<string, unknown>;
  if (t?.theme === 'letter') return [...new Set(Object.values(p.colors as Record<string, string>))].slice(0, 12);
  if (t?.theme === 'field_palette') {
    const colors = (p.palette as { colors: string[] }).colors;
    const step = Math.max(1, Math.floor(colors.length / 8));
    return colors.filter((_, i) => i % step === 0).slice(0, 10);
  }
  return [...new Set(((p.rules ?? []) as { color?: string }[]).flatMap((r) => (r.color ? [r.color] : [])))];
}
const fits = (e: ColorSchemeEntry) =>
  e.scheme.kind === 'letters' || e.scheme.kind === 'theme'
    ? `${e.scheme.alphabets.join(', ')} letters${e.scheme.kind === 'theme' && e.scheme.columns ? ' (alignment)' : ''}`
    : e.scheme.kind === 'categorical'
      ? 'categories'
      : `values (${e.scheme.kind})`;
function Swatches({ colors }: { colors: readonly string[] }) {
  return (
    <span className={cs.swatches}>
      {colors.map((c) => (
        <span key={c} className={cs.swatch} style={{ background: Color.toHex(Color.decodePrimary(c)) }} />
      ))}
    </span>
  );
}

/** Every kernel scheme: its name, label, what it fits, colorblind-safe, light and dark swatches, provenance. */
export function SchemeTable() {
  return (
    <table className={f.table}>
      <thead>
        <tr>
          <th>Scheme</th>
          <th>Fits</th>
          <th>Colorblind-safe</th>
          <th>Light</th>
          <th>Dark</th>
          <th>Provenance</th>
        </tr>
      </thead>
      <tbody>
        {KernelColorSchemes.map((e) => (
          <tr key={e.name}>
            <td>
              <code>{e.name}</code> {e.label}: {e.description}
            </td>
            <td>{fits(e)}</td>
            <td>{e.colorblind ? 'yes' : ''}</td>
            <td>
              <Swatches colors={swatches(e, false)} />
            </td>
            <td>
              <Swatches colors={swatches(e, true)} />
            </td>
            <td>{e.provenance}</td>
          </tr>
        ))}
      </tbody>
    </table>
  );
}

export const IN_DOCUMENT = `{ "kind": "layer", "ref": "family.rows",
  "params": { "type": "rows", "color": { "theme": "scheme", "params": { "name": "clustal_x" } } } }

{ "kind": "layer", "ref": "plddt",
  "params": { "type": "color", "data": "plddt-values",
    "color": { "theme": "scheme", "params": { "name": "viridis", "domain": [0, 100] } } } }`;
