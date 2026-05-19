/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Prihoda <david.prihoda@gmail.com>
 */

import { Structure } from '../../../mol-model/structure';
import { Representation, RepresentationContext, RepresentationParamsGetter } from '../../../mol-repr/representation';
import { ThemeRegistryContext, Theme } from '../../../mol-theme/theme';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation';
import { UnitsRepresentation } from '../units-representation';
import { PolymerTerminaParams, PolymerNtermConeVisual, PolymerCtermConeVisual } from '../visual/polymer-terminus-cone';
import { Color } from '../../../mol-util/color';
import { UniformColorTheme } from '../../../mol-theme/color/uniform';

const DefaultNtermColor = Color(0x0000AA); // dark blue
const DefaultCtermColor = Color(0xCC0000); // dark red

export const PolymerTerminusRepresentationParams = {
    ...PolymerTerminaParams,
    visuals: PD.MultiSelect(['nterm-cone', 'cterm-cone'] as TerminusVisualKey[], PD.arrayToOptions(['nterm-cone', 'cterm-cone'] as TerminusVisualKey[])),
    ntermColor: PD.Color(DefaultNtermColor),
    ctermColor: PD.Color(DefaultCtermColor),
};
export type PolymerTerminusRepresentationParams = typeof PolymerTerminusRepresentationParams
export type TerminusVisualKey = 'nterm-cone' | 'cterm-cone'

export function getPolymerTerminusParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(PolymerTerminusRepresentationParams);
}

export type PolymerTerminusRepresentation = StructureRepresentation<PolymerTerminusRepresentationParams>

function makeUniformTheme(base: Theme, color: Color): Theme {
    return {
        color: UniformColorTheme({}, { value: color, saturation: 0, lightness: 0 }),
        size: base.size,
    };
}

export function PolymerTerminusRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, PolymerTerminusRepresentationParams>): PolymerTerminusRepresentation {
    const ntermRepr = UnitsRepresentation('N-terminus cone', ctx, getParams, PolymerNtermConeVisual);
    const ctermRepr = UnitsRepresentation('C-terminus cone', ctx, getParams, PolymerCtermConeVisual);

    const base = Representation.createMulti('Polymer Terminus', ctx, getParams, StructureRepresentationStateBuilder, {
        'nterm-cone': () => ntermRepr,
        'cterm-cone': () => ctermRepr,
    } as any);

    const origSetTheme = base.setTheme.bind(base);
    base.setTheme = (theme: Theme) => {
        origSetTheme(theme);
        // override each visual's color after the base sets the shared theme
        const props = base.props as PD.Values<PolymerTerminusRepresentationParams>;
        const ntermColor = props?.ntermColor ?? DefaultNtermColor;
        const ctermColor = props?.ctermColor ?? DefaultCtermColor;
        ntermRepr.setTheme(makeUniformTheme(theme, ntermColor));
        ctermRepr.setTheme(makeUniformTheme(theme, ctermColor));
    };

    return base;
}

export const PolymerTerminusRepresentationProvider = StructureRepresentationProvider({
    name: 'polymer-terminus',
    label: 'Polymer Terminus',
    description: 'Displays cones at N/C termini and chain break endpoints indicating chain direction.',
    factory: PolymerTerminusRepresentation,
    getParams: getPolymerTerminusParams,
    defaultValues: PD.getDefaultValues(PolymerTerminusRepresentationParams),
    defaultColorTheme: { name: 'uniform', props: { value: DefaultNtermColor } },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (structure: Structure) => structure.polymerResidueCount > 0,
});
