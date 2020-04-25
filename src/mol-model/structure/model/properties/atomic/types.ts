/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementSymbol } from '../../types';
import { AtomNumber } from './measures';

/**
 * Enum of element symbols
 */
export const enum Elements {
    H = 'H', D = 'D', T = 'T', HE = 'HE', LI = 'LI', BE = 'BE', B = 'B', C = 'C', N = 'N', O = 'O', F = 'F', NE = 'NE', NA = 'NA', MG = 'MG', AL = 'AL', SI = 'SI', P = 'P', S = 'S', CL = 'CL', AR = 'AR', K = 'K', CA = 'CA', SC = 'SC', TI = 'TI', V = 'V', CR = 'CR', MN = 'MN', FE = 'FE', CO = 'CO', NI = 'NI', CU = 'CU', ZN = 'ZN', GA = 'GA', GE = 'GE', AS = 'AS', SE = 'SE', BR = 'BR', KR = 'KR', RB = 'RB', SR = 'SR', Y = 'Y', ZR = 'ZR', NB = 'NB', MO = 'MO', TC = 'TC', RU = 'RU', RH = 'RH', PD = 'PD', AG = 'AG', CD = 'CD', IN = 'IN', SN = 'SN', SB = 'SB', TE = 'TE', I = 'I', XE = 'XE', CS = 'CS', BA = 'BA', LA = 'LA', CE = 'CE', PR = 'PR', ND = 'ND', PM = 'PM', SM = 'SM', EU = 'EU', GD = 'GD', TB = 'TB', DY = 'DY', HO = 'HO', ER = 'ER', TM = 'TM', YB = 'YB', LU = 'LU', HF = 'HF', TA = 'TA', W = 'W', RE = 'RE', OS = 'OS', IR = 'IR', PT = 'PT', AU = 'AU', HG = 'HG', TL = 'TL', PB = 'PB', BI = 'BI', PO = 'PO', AT = 'AT', RN = 'RN', FR = 'FR', RA = 'RA', AC = 'AC', TH = 'TH', PA = 'PA', U = 'U', NP = 'NP', PU = 'PU', AM = 'AM', CM = 'CM', BK = 'BK', CF = 'CF', ES = 'ES', FM = 'FM', MD = 'MD', NO = 'NO', LR = 'LR', RF = 'RF', DB = 'DB', SG = 'SG', BH = 'BH', HS = 'HS', MT = 'MT', DS = 'DS', RG = 'RG', CN = 'CN', NH = 'NH', FL = 'FL', MC = 'MC', LV = 'LV', TS = 'TS', OG = 'OG'
}

export const ElementNames: { [k: string]: string } = {
    H: 'Hydrogen', HE: 'Helium', LI: 'Lithium', BE: 'Beryllium', B: 'Boron', C: 'Carbon', N: 'Nitrogen', O: 'Oxygen', F: 'Fluorine', NE: 'Neon', NA: 'Sodium', MG: 'Magnesium', AL: 'Aluminum', SI: 'Silicon', P: 'Phosphorus', S: 'Sulfur', CL: 'Chlorine', AR: 'Argon', K: 'Potassium', CA: 'Calcium', SC: 'Scandium', TI: 'Titanium', V: 'Vanadium', CR: 'Chromium', MN: 'Manganese', FE: 'Iron', CO: 'Cobalt', NI: 'Nickel', CU: 'Copper', ZN: 'Zinc', GA: 'Gallium', GE: 'Germanium', AS: 'Arsenic', SE: 'Selenium', BR: 'Bromine', KR: 'Krypton', RB: 'Rubidium', SR: 'Strontium', Y: 'Yttrium', ZR: 'Zirconium', NB: 'Niobium', MO: 'Molybdenum', TC: 'Technetium', RU: 'Ruthenium', RH: 'Rhodium', PD: 'Palladium', AG: 'Silver', CD: 'Cadmium', IN: 'Indium', SN: 'Tin', SB: 'Antimony', TE: 'Tellurium', I: 'Iodine', XE: 'Xenon', CS: 'Cesium', BA: 'Barium', LA: 'Lanthanum', CE: 'Cerium', PR: 'Praseodymium', ND: 'Neodymium', PM: 'Promethium', SM: 'Samarium', EU: 'Europium', GD: 'Gadolinium', TB: 'Terbium', DY: 'Dysprosium', HO: 'Holmium', ER: 'Erbium', TM: 'Thulium', YB: 'Ytterbium', LU: 'Lutetium', HF: 'Hafnium', TA: 'Tantalum', W: 'Wolfram', RE: 'Rhenium', OS: 'Osmium', IR: 'Iridium', PT: 'Platinum', AU: 'Gold', HG: 'Mercury', TL: 'Thallium', PB: 'Lead', BI: 'Bismuth', PO: 'Polonium', AT: 'Astatine', RN: 'Radon', FR: 'Francium', RA: 'Radium', AC: 'Actinium', TH: 'Thorium', PA: 'Protactinium', U: 'Uranium', NP: 'Neptunium', PU: 'Plutonium', AM: 'Americium', CM: 'Curium', BK: 'Berkelium', CF: 'Californium', ES: 'Einsteinium', FM: 'Fermium', MD: 'Mendelevium', NO: 'Nobelium', LR: 'Lawrencium', RF: 'Rutherfordium', DB: 'Dubnium', SG: 'Seaborgium', BH: 'Bohrium', HS: 'Hassium', MT: 'Meitnerium', DS: 'Darmstadtium', RG: 'Roentgenium', CN: 'Copernicium', NH: 'Nihonium', FL: 'Flerovium', MC: 'Moscovium', LV: 'Livermorium', TS: 'Tennessine', OG: 'Oganesson'
};

export const AlkaliMetals = new Set<ElementSymbol>(['LI', 'NA', 'K', 'RB', 'CS', 'FR'] as ElementSymbol[]);
export function isAlkaliMetal(element: ElementSymbol) { return AlkaliMetals.has(element); }

export const AlkalineEarthMetals = new Set<ElementSymbol>(['BE', 'MG', 'CA', 'SR', 'BA', 'RA'] as ElementSymbol[]);
export function isAlkalineEarthMetal(element: ElementSymbol) { return AlkalineEarthMetals.has(element); }

export const PolyatomicNonmetals = new Set<ElementSymbol>(['C', 'P', 'S', 'SE'] as ElementSymbol[]);
export function isPolyatomicNonmetal(element: ElementSymbol) { return PolyatomicNonmetals.has(element); }

export const DiatomicNonmetals = new Set<ElementSymbol>(['H', 'N', 'O', 'F', 'CL', 'BR', 'I'] as ElementSymbol[]);
export function isDiatomicNonmetal(element: ElementSymbol) { return DiatomicNonmetals.has(element); }

export const NobleGases = new Set<ElementSymbol>(['HE', 'NE', 'AR', 'KR', 'XE', 'RN'] as ElementSymbol[]);
export function isNobleGas(element: ElementSymbol) { return NobleGases.has(element); }

export const PostTransitionMetals = new Set<ElementSymbol>(['ZN', 'GA', 'CD', 'IN', 'SN', 'HG', 'TI', 'PB', 'BI', 'PO', 'CN'] as ElementSymbol[]);
export function isPostTransitionMetal(element: ElementSymbol) { return PostTransitionMetals.has(element); }

export const Metalloids = new Set<ElementSymbol>(['B', 'SI', 'GE', 'AS', 'SB', 'TE', 'AT'] as ElementSymbol[]);
export function isMetalloid(element: ElementSymbol) { return Metalloids.has(element); }

export const Halogens = new Set<ElementSymbol>(['F', 'CL', 'BR', 'I', 'AT'] as ElementSymbol[]);
export function isHalogen(element: ElementSymbol) { return Halogens.has(element); }

export function isTransitionMetal(element: ElementSymbol) {
    const no = AtomNumber(element);
    return (
        (no >= 21 && no <= 29) ||
        (no >= 39 && no <= 47) ||
        (no >= 72 && no <= 79) ||
        (no >= 104 && no <= 108)
    );
}

export function isLanthanide (element: ElementSymbol) {
    const no = AtomNumber(element);
    return no >= 57 && no <= 71;
}

export function isActinide (element: ElementSymbol) {
    const no = AtomNumber(element);
    return no >= 89 && no <= 103;
}

export function isMetal(element: ElementSymbol) {
    return (
        isAlkaliMetal(element) ||
        isAlkalineEarthMetal(element) ||
        isLanthanide(element) ||
        isActinide(element) ||
        isTransitionMetal(element) ||
        isPostTransitionMetal(element)
    );
}

export function isNonmetal(element: ElementSymbol) {
    return (
        isDiatomicNonmetal(element) ||
        isPolyatomicNonmetal(element) ||
        isNobleGas(element)
    );
}