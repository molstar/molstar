/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html'
import { Canvas3D } from '../../mol-canvas3d/canvas3d';
import { CIF, CifFrame } from '../../mol-io/reader/cif'
import { Model, Structure } from '../../mol-model/structure';
import { ColorTheme } from '../../mol-theme/color';
import { SizeTheme } from '../../mol-theme/size';
import { CartoonRepresentationProvider } from '../../mol-repr/structure/representation/cartoon';
import { trajectoryFromMmCIF } from '../../mol-model-formats/structure/mmcif';
import { MolecularSurfaceRepresentationProvider } from '../../mol-repr/structure/representation/molecular-surface';
import { BallAndStickRepresentationProvider } from '../../mol-repr/structure/representation/ball-and-stick';
import { GaussianSurfaceRepresentationProvider } from '../../mol-repr/structure/representation/gaussian-surface';
import { ComputedSecondaryStructure } from '../../mol-model-props/computed/secondary-structure';
import { resizeCanvas } from '../../mol-canvas3d/util';

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

function getBallAndStickRepr() {
    return BallAndStickRepresentationProvider.factory(reprCtx, BallAndStickRepresentationProvider.getParams)
}

function getMolecularSurfaceRepr() {
    return MolecularSurfaceRepresentationProvider.factory(reprCtx, MolecularSurfaceRepresentationProvider.getParams)
}

function getGaussianSurfaceRepr() {
    return GaussianSurfaceRepresentationProvider.factory(reprCtx, GaussianSurfaceRepresentationProvider.getParams)
}

async function init() {
    const cif = await downloadFromPdb('1crn')
    const models = await getModels(cif)
    const structure = await getStructure(models[0])
    console.time('computeDSSP')
    await ComputedSecondaryStructure.attachFromCifOrCompute(structure)
    console.timeEnd('computeDSSP');

    const show = {
        cartoon: false,
        ballAndStick: true,
        molecularSurface: true,
        gaussianSurface: false,
    }

    const cartoonRepr = getCartoonRepr()
    const ballAndStickRepr = getBallAndStickRepr()
    const molecularSurfaceRepr = getMolecularSurfaceRepr()
    const gaussianSurfaceRepr = getGaussianSurfaceRepr()

    if (show.cartoon) {
        cartoonRepr.setTheme({
            color: reprCtx.colorThemeRegistry.create('secondary-structure', { structure }),
            size: reprCtx.sizeThemeRegistry.create('uniform', { structure })
        })
        await cartoonRepr.createOrUpdate({ ...CartoonRepresentationProvider.defaultValues, quality: 'auto' }, structure).run()
    }

    if (show.ballAndStick) {
        ballAndStickRepr.setTheme({
            color: reprCtx.colorThemeRegistry.create('secondary-structure', { structure }),
            size: reprCtx.sizeThemeRegistry.create('uniform', { structure })
        })
        await ballAndStickRepr.createOrUpdate({ ...BallAndStickRepresentationProvider.defaultValues, quality: 'auto' }, structure).run()
    }

    if (show.molecularSurface) {
        molecularSurfaceRepr.setTheme({
            color: reprCtx.colorThemeRegistry.create('secondary-structure', { structure }),
            size: reprCtx.sizeThemeRegistry.create('physical', { structure })
        })
        console.time('molecular surface')
        await molecularSurfaceRepr.createOrUpdate({ ...MolecularSurfaceRepresentationProvider.defaultValues, quality: 'custom', alpha: 0.5, flatShaded: true, doubleSided: true, resolution: 0.3 }, structure).run()
        console.timeEnd('molecular surface')
    }

    if (show.gaussianSurface) {
        gaussianSurfaceRepr.setTheme({
            color: reprCtx.colorThemeRegistry.create('secondary-structure', { structure }),
            size: reprCtx.sizeThemeRegistry.create('physical', { structure })
        })
        console.time('gaussian surface')
        await gaussianSurfaceRepr.createOrUpdate({ ...GaussianSurfaceRepresentationProvider.defaultValues, quality: 'custom', alpha: 1.0, flatShaded: true, doubleSided: true, resolution: 0.3 }, structure).run()
        console.timeEnd('gaussian surface')
    }

    if (show.cartoon) canvas3d.add(cartoonRepr)
    if (show.ballAndStick) canvas3d.add(ballAndStickRepr)
    if (show.molecularSurface) canvas3d.add(molecularSurfaceRepr)
    if (show.gaussianSurface) canvas3d.add(gaussianSurfaceRepr)
    canvas3d.resetCamera()
    // canvas3d.setProps({ trackball: { ...canvas3d.props.trackball, spin: true } })
}

init()