/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html'
import { Canvas3D } from 'mol-canvas3d/canvas3d';
import CIF, { CifFrame } from 'mol-io/reader/cif'
import { Model, Structure, StructureElement, Unit } from 'mol-model/structure';
import { ColorTheme, LocationColor } from 'mol-theme/color';
import { SizeTheme } from 'mol-theme/size';
import { CartoonRepresentationProvider } from 'mol-repr/structure/representation/cartoon';
import { trajectoryFromMmCIF } from 'mol-model-formats/structure/mmcif';
import { AccessibleSurfaceArea } from 'mol-model/structure/structure/accessible-surface-area';
import { Color, ColorScale } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { ThemeDataContext } from 'mol-theme/theme';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { ColorListName, ColorListOptions } from 'mol-util/color/scale';
import { resizeCanvas } from 'mol-canvas3d/util';

const parent = document.getElementById('app')!
parent.style.width = '100%'
parent.style.height = '100%'

const canvas = document.createElement('canvas')
parent.appendChild(canvas)
resizeCanvas(canvas, parent)

const canvas3d = Canvas3D.fromCanvas(canvas)
canvas3d.animate()


async function parseCif(data: string|Uint8Array) {
    const comp = CIF.parse(data);
    const parsed = await comp.run();
    if (parsed.isError) throw parsed;
    return parsed.result;
}

async function downloadCif(url: string, isBinary: boolean) {
    const data = await fetch(url);
    return parseCif(isBinary ? new Uint8Array(await data.arrayBuffer()) : await data.text());
}

async function downloadFromPdb(pdb: string) {
    // const parsed = await downloadCif(`https://files.rcsb.org/download/${pdb}.cif`, false);
    const parsed = await downloadCif(`https://webchem.ncbr.muni.cz/ModelServer/static/bcif/${pdb}`, true);
    return parsed.blocks[0];
}

async function getModels(frame: CifFrame) {
    return await trajectoryFromMmCIF(frame).run();
}

async function getStructure(model: Model) {
    return Structure.ofModel(model);
}

const reprCtx = {
    colorThemeRegistry: ColorTheme.createRegistry(),
    sizeThemeRegistry: SizeTheme.createRegistry()
}
function getCartoonRepr() {
    return CartoonRepresentationProvider.factory(reprCtx, CartoonRepresentationProvider.getParams)
}

let accessibleSurfaceArea: AccessibleSurfaceArea;
async function init(props = {}) {
    const cif = await downloadFromPdb(
        // '3j3q'
        // '1aon'
        // '1acj'
        // '1pga'
        '1brr'
        // '1hrc'
    )
    const models = await getModels(cif)
    const structure = await getStructure(models[0])

    // async compute ASA
    accessibleSurfaceArea = await AccessibleSurfaceArea.compute(structure)

    const cartoonRepr = getCartoonRepr()

    // create color theme
    cartoonRepr.setTheme({
        color: AccessibleSurfaceAreaColorTheme(reprCtx, { ...PD.getDefaultValues(AccessibleSurfaceAreaColorThemeParams), ...props }),
        size: reprCtx.sizeThemeRegistry.create('uniform', { structure })
    })
    await cartoonRepr.createOrUpdate({ ...CartoonRepresentationProvider.defaultValues, quality: 'auto' }, structure).run()

    canvas3d.add(cartoonRepr)
    canvas3d.resetCamera()
}

init()

const DefaultColor = Color(0xFFFFFF)
const Description = 'Assigns a color based on the relative accessible surface area of a residue.'

export const AccessibleSurfaceAreaColorThemeParams = {
    list: PD.ColorScale<ColorListName>('Rainbow', ColorListOptions)
}
export type AccessibleSurfaceAreaColorThemeParams = typeof AccessibleSurfaceAreaColorThemeParams
export function getAccessibleSurfaceAreaColorThemeParams(ctx: ThemeDataContext) {
    return AccessibleSurfaceAreaColorThemeParams // TODO return copy
}

export function AccessibleSurfaceAreaColorTheme(ctx: ThemeDataContext, props: PD.Values<AccessibleSurfaceAreaColorThemeParams>): ColorTheme<AccessibleSurfaceAreaColorThemeParams> {
    let color: LocationColor = () => DefaultColor
    const scale = ColorScale.create({
        listOrName: props.list,
        minLabel: '0.0 (buried)',
        maxLabel: '1.0 (exposed)',
        domain: [0.0, 1.0]
    })
    color = (location: Location): Color => {
        if (StructureElement.isLocation(location)) {
            if (Unit.isAtomic(location.unit)) {
                const value = accessibleSurfaceArea.relativeAccessibleSurfaceArea![location.unit.residueIndex[location.element]];
                return value !== AccessibleSurfaceArea.VdWLookup[0] /* signals missing value */ ? scale.color(value) : DefaultColor;
            }
        }

        return DefaultColor
    }

    return {
        factory: AccessibleSurfaceAreaColorTheme,
        granularity: 'group',
        color,
        props,
        description: Description,
        legend: scale ? scale.legend : undefined
    }
}

export const AccessibleSurfaceAreaColorThemeProvider: ColorTheme.Provider<AccessibleSurfaceAreaColorThemeParams> = {
    label: 'Accessible Surface Area',
    factory: AccessibleSurfaceAreaColorTheme,
    getParams: getAccessibleSurfaceAreaColorThemeParams,
    defaultValues: PD.getDefaultValues(AccessibleSurfaceAreaColorThemeParams),
    isApplicable: (ctx: ThemeDataContext) => !!ctx.structure
}