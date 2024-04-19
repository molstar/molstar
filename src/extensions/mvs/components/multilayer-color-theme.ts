/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Location } from '../../../mol-model/location';
import { Bond, Structure, StructureElement } from '../../../mol-model/structure';
import { ColorTheme, LocationColor } from '../../../mol-theme/color';
import { ThemeDataContext } from '../../../mol-theme/theme';
import { Color } from '../../../mol-util/color';
import { ColorNames } from '../../../mol-util/color/names';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { stringToWords } from '../../../mol-util/string';
import { isMVSStructure } from './is-mvs-model-prop';
import { ElementSet, SelectorParams, isSelectorAll } from './selector';


/** Special value that can be used as color with null-like semantic (i.e. "no color provided").
 * By some lucky coincidence, Mol* treats -1 as white. */
export const NoColor = Color(-1);

/** Return true if `color` is a real color, false if it is `NoColor`. */
function isValidColor(color: Color): boolean {
    return color >= 0;
}

const DefaultBackgroundColor = ColorNames.white;

/** Parameter definition for color theme "Multilayer" */
export function makeMultilayerColorThemeParams(colorThemeRegistry: ColorTheme.Registry, ctx: ThemeDataContext) {
    const colorThemeInfo = {
        help: (value: { name: string, params: {} }) => {
            const { name, params } = value;
            const p = colorThemeRegistry.get(name);
            const ct = p.factory({}, params);
            return { description: ct.description, legend: ct.legend };
        }
    };
    const nestedThemeTypes = colorThemeRegistry.types.filter(([name, label, category]) => name !== MultilayerColorThemeName && colorThemeRegistry.get(name).isApplicable(ctx)); // Adding 'multilayer' theme itself would cause infinite recursion
    return {
        layers: PD.ObjectList(
            {
                theme: PD.Mapped<any>(
                    'uniform',
                    nestedThemeTypes,
                    name => PD.Group<any>(colorThemeRegistry.get(name).getParams({ structure: Structure.Empty })),
                    colorThemeInfo),
                selection: SelectorParams,
            },
            obj => stringToWords(obj.theme.name),
            { description: 'A list of layers, each defining a color theme. The last listed layer is the top layer (applies first). If the top layer does not provide color for a location or its selection does not cover the location, the underneath layers will apply.' }),
        background: PD.Color(DefaultBackgroundColor, { description: 'Color for elements where no layer applies' }),
    };
}
/** Parameter definition for color theme "Multilayer" */
export type MultilayerColorThemeParams = ReturnType<typeof makeMultilayerColorThemeParams>

/** Parameter values for color theme "Multilayer" */
export type MultilayerColorThemeProps = PD.Values<MultilayerColorThemeParams>

/** Default values for `MultilayerColorThemeProps` */
export const DefaultMultilayerColorThemeProps: MultilayerColorThemeProps = { layers: [], background: DefaultBackgroundColor };


/** Return color theme that assigns colors based on a list of nested color themes (layers).
 * The last layer in the list whose selection covers the given location
 * and which provides a valid (non-negative) color value will be used.
 * If a nested theme provider has `ensureCustomProperties` methods, these will not be called automatically
 * (the caller must ensure that any required custom properties be attached). */
function makeMultilayerColorTheme(ctx: ThemeDataContext, props: MultilayerColorThemeProps, colorThemeRegistry: ColorTheme.Registry): ColorTheme<MultilayerColorThemeParams> {
    const colorLayers: { color: LocationColor, elementSet: ElementSet | undefined }[] = []; // undefined elementSet means 'all'
    for (let i = props.layers.length - 1; i >= 0; i--) { // iterate from end to get top layer first, bottom layer last
        const layer = props.layers[i];
        const themeProvider = colorThemeRegistry.get(layer.theme.name);
        if (!themeProvider) {
            console.warn(`Skipping color theme '${layer.theme.name}', cannot find it in registry.`);
            continue;
        }
        if (themeProvider.ensureCustomProperties?.attach) {
            console.warn(`Multilayer color theme: layer "${themeProvider.name}" has ensureCustomProperties.attach method, but Multilayer color theme does not call it. If the layer does not work, make sure you call ensureCustomProperties.attach somewhere.`);
        }
        const theme = themeProvider.factory(ctx, layer.theme.params);
        switch (theme.granularity) {
            case 'uniform':
            case 'instance':
            case 'group':
            case 'groupInstance':
            case 'vertex':
            case 'vertexInstance':
                const elementSet = isSelectorAll(layer.selection) ? undefined : ElementSet.fromSelector(ctx.structure, layer.selection); // treating 'all' specially for performance reasons (it's expected to be used most often)
                colorLayers.push({ color: theme.color, elementSet });
                break;
            default:
                console.warn(`Skipping color theme '${layer.theme.name}', cannot process granularity '${theme.granularity}'`);
        }
    };

    function structureElementColor(loc: StructureElement.Location, isSecondary: boolean): Color {
        for (const layer of colorLayers) {
            const matches = !layer.elementSet || ElementSet.has(layer.elementSet, loc);
            if (!matches) continue;
            const color = layer.color(loc, isSecondary);
            if (!isValidColor(color)) continue;
            return color;
        }
        return props.background;
    }
    const auxLocation = StructureElement.Location.create(ctx.structure);

    const color: LocationColor = (location: Location, isSecondary: boolean) => {
        if (StructureElement.Location.is(location)) {
            return structureElementColor(location, isSecondary);
        } else if (Bond.isLocation(location)) {
            // this will be applied for each bond twice, to get color of each half (a* refers to the adjacent atom, b* to the opposite atom)
            auxLocation.unit = location.aUnit;
            auxLocation.element = location.aUnit.elements[location.aIndex];
            return structureElementColor(auxLocation, isSecondary);
        }
        return props.background;
    };

    return {
        factory: (ctx_, props_) => makeMultilayerColorTheme(ctx_, props_, colorThemeRegistry),
        granularity: 'group',
        preferSmoothing: true,
        color: color,
        props: props,
        description: 'Combines colors from multiple color themes.',
    };
}


/** Unique name for "Multilayer" color theme */
export const MultilayerColorThemeName = 'mvs-multilayer';

/** A thingy that is needed to register color theme "Multilayer" */
export function makeMultilayerColorThemeProvider(colorThemeRegistry: ColorTheme.Registry): ColorTheme.Provider<MultilayerColorThemeParams, typeof MultilayerColorThemeName> {
    return {
        name: MultilayerColorThemeName,
        label: 'MVS Multi-layer',
        category: ColorTheme.Category.Misc,
        factory: (ctx, props) => makeMultilayerColorTheme(ctx, props, colorThemeRegistry),
        getParams: (ctx: ThemeDataContext) => makeMultilayerColorThemeParams(colorThemeRegistry, ctx),
        defaultValues: DefaultMultilayerColorThemeProps,
        isApplicable: (ctx: ThemeDataContext) => !!ctx.structure && isMVSStructure(ctx.structure),
    };
}
