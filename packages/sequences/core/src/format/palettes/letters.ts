/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * The letter dictionaries of the `letter` theme, keyed by one-letter code, upper case (seqspec:format#themes).
 *
 * Jalview's residue colour schemes (https://www.jalview.org/help/html/colourSchemes/), copied from Nightingale
 * (MIT, (c) 2018 Plotly, Inc, (c) 2023 EMBL-European Bioinformatics Institute; nightingale-msa 5.11.0,
 * `packages/nightingale-msa/src/colorschemes/{zappo,taylor,hydrophobicity,helix,strand,turn,buried,nucleotide,
 * purine}.ts`). Nightingale's white placeholders (`#fff` for letters a scheme does not color) are left out: those
 * letters are not colored.
 */

import { Color, colorDict } from "./color.ts";
import { ColorDictionaries, type ColorDictionaryName, ElementSymbolColors } from "./dictionaries.ts";

const dict = colorDict;

export const LetterDictionaries = {
  Zappo: dict("#ffafaf A I L M V #6464ff R H K #00ff00 N Q S T #ff0000 D E #ffff00 C #ff00ff G P #ffc800 F W Y"),
  Taylor: dict(`
    #ccff00 A #0000ff R #cc00ff N #ff0000 D #ffff00 C #ff00cc Q #ff0066 E #ff9900 G #0066ff H #66ff00 I #33ff00 L
    #6600ff K #00ff00 M #00ff66 F #ffcc00 P #ff3300 S #ff6600 T #00ccff W #00ffcc Y #99ff00 V
  `),
  Hydrophobicity: dict(`
    #ad0052 A #0c00f3 B D E N Q Z #c2003d C #cb0034 F #6a0095 G #1500ea H #ff0000 I #0000ff K R #ea0015 L #b0004f M
    #4600b9 P #5e00a1 S #61009e T #f60009 V #5b00a4 W #680097 X #4f00b0 Y
  `),
  HelixPropensity: dict(`
    #e718e7 A #6f906f R #1be41b N #778877 D #23dc23 C #926d92 Q #ff00ff E #00ff00 G P #758a75 H X #8a758a I W
    #ae51ae L #a05fa0 K #ef10ef M #986798 F #36c936 S #47b847 T #21de21 Y #857a85 V #49b649 B #c936c9 Z
  `),
  StrandPropensity: dict(`
    #5858a7 A #6b6b94 R #64649b N #2121de D #9d9d62 C T #8c8c73 Q #0000ff E #4949b6 G S #60609f H #ecec13 I
    #b2b24d L #4747b8 K Z #82827d M #c2c23d F #2323dc P #c0c03f W #d3d32c Y #ffff00 V #4343bc B #797986 X
  `),
  TurnPropensity: dict(`
    #2cd3d3 A #708f8f R H #ff0000 N G #e81717 D #a85757 C #3fc0c0 Q #778888 E #00ffff I #1ce3e3 L #7e8181 K
    #1ee1e1 M F #f60909 P #e11e1e S #738c8c T W #9d6262 Y #07f8f8 V #f30c0c B #7c8383 X #5ba4a4 Z
  `),
  BuriedIndex: dict(`
    #00a35c A #00fc03 R #00eb14 N D B #0000ff C #00f10e Q E Z #009d62 G #00d52a H S #0054ab I #007b84 L #00ff00 K
    #009768 M #008778 F #00e01f P #00db24 T #00a857 W #00e619 Y #005fa0 V #00b649 X
  `),
  Nucleotide: dict("#64f73f A #ffb340 C #eb413c G #3c88ee T U"),
  PurinePyrimidine: dict("#ff83fa A G R #40e0d0 C T U Y"),
};

export type LetterDictionaryName = keyof typeof LetterDictionaries;
export const LetterDictionaryNames = Object.keys(LetterDictionaries) as readonly LetterDictionaryName[];

/**
 * The colors of a `letter` theme: `colors` a dictionary name (a letter dictionary or MVS's; `ElementSymbol` with
 * carbon, as `letter` colors atoms with it), a `{letter: Color}` object, or `null` for the default of `alphabet`
 * (`null`: a ligand's atoms): `Zappo` for protein, `Nucleotide` for DNA and RNA, `ElementSymbol` with carbon for
 * atoms. `undefined` for an unknown name and for the `generic` alphabet's default (it has none).
 */
export function letterColors(
  colors: string | Readonly<Record<string, string>> | null,
  alphabet: "protein" | "dna" | "rna" | "generic" | null,
): Readonly<Record<string, Color>> | undefined {
  if (colors === null) {
    if (alphabet === "protein") return LetterDictionaries.Zappo;
    if (alphabet === "dna" || alphabet === "rna") return LetterDictionaries.Nucleotide;
    return alphabet === null ? ElementSymbolColors : undefined;
  }
  if (typeof colors === "string") {
    if (colors === "ElementSymbol") return ElementSymbolColors;
    if (Object.hasOwn(LetterDictionaries, colors)) return LetterDictionaries[colors as LetterDictionaryName];
    if (Object.hasOwn(ColorDictionaries, colors)) return ColorDictionaries[colors as ColorDictionaryName];
    return undefined;
  }
  const out: Record<string, Color> = {};
  for (const [key, color] of Object.entries(colors)) out[key] = Color.decodePrimary(color);
  return out;
}
