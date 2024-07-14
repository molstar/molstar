/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Herman Bergwerf <post@hbergwerf.nl>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AssetManager } from '../../mol-util/assets';
import { Canvas3D, Canvas3DProps, Canvas3DContext } from '../../mol-canvas3d/canvas3d';
import { resizeCanvas } from '../../mol-canvas3d/util';
import { parseSdf as _parseSdf, SdfFile } from '../../mol-io/reader/sdf/parser';
import { Structure } from '../../mol-model/structure';
import { trajectoryFromSdf } from '../../mol-model-formats/structure/sdf';
import { ColorTheme } from '../../mol-theme/color';
import { SizeTheme } from '../../mol-theme/size';
import { RepresentationContext } from '../../mol-repr/representation';
import { BallAndStickRepresentationProvider } from '../../mol-repr/structure/representation/ball-and-stick';
import { StructureRepresentationProvider } from '../../mol-repr/structure/representation';
import { ParamDefinition as PD } from '../../mol-util/param-definition';

async function downloadPubChemSdf(cid: number) {
    const root = 'https://pubchem.ncbi.nlm.nih.gov/rest';
    const r = await fetch(`${root}/pug/compound/cid/${cid}/sdf?record_type=3d`);
    return await r.text();
}

async function parseSdf(sdf: string): Promise<SdfFile> {
    const parsed = await _parseSdf(sdf).run();
    if (parsed.isError) {
        throw parsed;
    } else {
        return parsed.result;
    }
}

function getRepr<P extends PD.Params>(provider: StructureRepresentationProvider<P>, reprCtx: RepresentationContext) {
    return provider.factory(reprCtx, provider.getParams);
}

function main() {
    const parent = document.body;
    const canvas = document.createElement('canvas');
    parent.appendChild(canvas);
    resizeCanvas(canvas, parent);

    const ctx = Canvas3DContext.fromCanvas(canvas, new AssetManager());
    const cids = [
        6440397,
        2244,
        2519,
        241
    ];

    const reprCtx: RepresentationContext = {
        webgl: ctx.webgl,
        colorThemeRegistry: ColorTheme.createRegistry(),
        sizeThemeRegistry: SizeTheme.createRegistry()
    };

    const rect = parent.getBoundingClientRect();
    const w = Math.floor(rect.width / 2);
    const h = Math.floor(rect.height / 2);
    for (let i = 0; i < cids.length; i++) {
        const ix = i % 2;
        const iy = Math.floor(i / 2);
        addViewer(ctx, reprCtx, cids[i], {
            name: 'static-frame',
            params: {
                x: ix * w,
                y: iy * h,
                width: ix === 1 ? rect.width - w : w,
                height: iy === 1 ? rect.height - h : h,
            }
        }, {
            adjustCylinderLength: ix === 0
        });
    }
}

async function addViewer(
    ctx: Canvas3DContext,
    reprCtx: RepresentationContext,
    cid: number,
    viewport: Canvas3DProps['viewport'],
    reprValues: Partial<typeof BallAndStickRepresentationProvider.defaultValues>
) {
    const file = await parseSdf(await downloadPubChemSdf(cid));
    const models = await trajectoryFromSdf(file.compounds[0]).run();
    const structure = Structure.ofModel(models.representative);
    const canvas3d = Canvas3D.create(ctx, { viewport });

    const repr = getRepr(BallAndStickRepresentationProvider, reprCtx);
    repr.setTheme({
        color: reprCtx.colorThemeRegistry.create('element-symbol', { structure }, {
            carbonColor: { name: 'element-symbol' }
        } as ColorTheme.BuiltInParams<'element-symbol'>),
        size: reprCtx.sizeThemeRegistry.create('physical', { structure })
    });

    await repr.createOrUpdate({
        ...BallAndStickRepresentationProvider.defaultValues,
        ...reprValues
    }, structure).run();

    canvas3d.add(repr);
    canvas3d.animate();
}

main();
