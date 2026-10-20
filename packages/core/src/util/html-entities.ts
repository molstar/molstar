/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// HTML character references in text (Markdown's, a label's): `&#160;`, `&#xA0;`, and the named ones descriptions and
// labels use. DOM-free (headless, workers). Not the whole HTML table (2000+ names): a name not here stays as written.

/** A Map, not an object literal: `&constructor;` must not find `Object.prototype.constructor`. */
const Named = new Map<string, string>([
  // Markup and spaces
  ['amp', '&'],
  ['lt', '<'],
  ['gt', '>'],
  ['quot', '"'],
  ['apos', "'"],
  ['nbsp', ' '],
  ['ensp', ' '],
  ['emsp', ' '],
  ['thinsp', ' '],
  ['shy', '­'],
  ['zwj', '‍'],
  ['zwnj', '‌'],
  // Punctuation
  ['ndash', '–'],
  ['mdash', '—'],
  ['hellip', '…'],
  ['lsquo', '‘'],
  ['rsquo', '’'],
  ['ldquo', '“'],
  ['rdquo', '”'],
  ['laquo', '«'],
  ['raquo', '»'],
  ['bull', '•'],
  ['middot', '·'],
  ['para', '¶'],
  ['sect', '§'],
  ['dagger', '†'],
  ['Dagger', '‡'],
  ['prime', '′'],
  ['Prime', '″'],
  ['iexcl', '¡'],
  ['iquest', '¿'],
  // Symbols
  ['copy', '©'],
  ['reg', '®'],
  ['trade', '™'],
  ['deg', '°'],
  ['plusmn', '±'],
  ['times', '×'],
  ['divide', '÷'],
  ['minus', '−'],
  ['le', '≤'],
  ['ge', '≥'],
  ['ne', '≠'],
  ['asymp', '≈'],
  ['infin', '∞'],
  ['micro', 'µ'],
  ['permil', '‰'],
  ['sup1', '¹'],
  ['sup2', '²'],
  ['sup3', '³'],
  ['frac12', '½'],
  ['frac14', '¼'],
  ['frac34', '¾'],
  ['euro', '€'],
  ['pound', '£'],
  ['yen', '¥'],
  ['cent', '¢'],
  ['check', '✓'],
  ['star', '☆'],
  // Arrows
  ['larr', '←'],
  ['uarr', '↑'],
  ['rarr', '→'],
  ['darr', '↓'],
  ['harr', '↔'],
  ['lArr', '⇐'],
  ['rArr', '⇒'],
  ['hArr', '⇔'],
  // Letters
  ['Aring', 'Å'],
  ['aring', 'å'],
  ['angst', 'Å'],
  ['auml', 'ä'],
  ['ouml', 'ö'],
  ['uuml', 'ü'],
  ['Auml', 'Ä'],
  ['Ouml', 'Ö'],
  ['Uuml', 'Ü'],
  ['szlig', 'ß'],
  ['eacute', 'é'],
  ['egrave', 'è'],
  ['aacute', 'á'],
  ['iacute', 'í'],
  ['oacute', 'ó'],
  ['uacute', 'ú'],
  ['ntilde', 'ñ'],
  ['ccedil', 'ç'],
  // Greek
  ['alpha', 'α'],
  ['beta', 'β'],
  ['gamma', 'γ'],
  ['delta', 'δ'],
  ['epsilon', 'ε'],
  ['zeta', 'ζ'],
  ['eta', 'η'],
  ['theta', 'θ'],
  ['kappa', 'κ'],
  ['lambda', 'λ'],
  ['mu', 'μ'],
  ['nu', 'ν'],
  ['xi', 'ξ'],
  ['pi', 'π'],
  ['rho', 'ρ'],
  ['sigma', 'σ'],
  ['tau', 'τ'],
  ['phi', 'φ'],
  ['chi', 'χ'],
  ['psi', 'ψ'],
  ['omega', 'ω'],
  ['Gamma', 'Γ'],
  ['Delta', 'Δ'],
  ['Theta', 'Θ'],
  ['Lambda', 'Λ'],
  ['Pi', 'Π'],
  ['Sigma', 'Σ'],
  ['Phi', 'Φ'],
  ['Psi', 'Ψ'],
  ['Omega', 'Ω'],
]);

/** A character reference: decimal (`&#160;`), hexadecimal (`&#xA0;`) or named (`&nbsp;`). */
const Reference = /&(?:#(\d{1,7})|#[xX]([0-9a-fA-F]{1,6})|([a-zA-Z][a-zA-Z0-9]{1,31}));/g;

/** One reference's character: a numeric one outside Unicode, a surrogate or zero is U+FFFD (CommonMark). */
function decodeReference(whole: string, dec?: string, hex?: string, name?: string): string {
  if (name !== undefined) return Named.get(name) ?? whole;
  const code = dec !== undefined ? Number.parseInt(dec, 10) : Number.parseInt(hex!, 16);
  return code > 0 && code <= 0x10ffff && (code < 0xd800 || code > 0xdfff) ? String.fromCodePoint(code) : '�';
}

/** The text with its character references as the characters they name; a name not known stays as written. */
export function decodeHtmlEntities(text: string): string {
  return text.includes('&') ? text.replace(Reference, decodeReference) : text;
}
