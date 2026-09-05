/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Plugin integration for the HFF (EMDB-SFF / HDF5 segmentation) format.
 *
 * Mirrors the g3d extension's structure (see issue #886): a vendored HDF5
 * reader provides binary I/O, a typed parser produces SffData, and this
 * module wires it into mol*'s state-tree:
 *   Binary -> ParseHff -> ShapeFromSff -> Shape.Provider -> ShapeRepresentation3D
 *
 * This is the same shape of pipeline the built-in shape formats use
 * (`mol-plugin-state/formats/shape.ts` for ply/obj/vtp), so HFF gets the generic
 * `ShapeRepresentation3D` and, via `parseRaw`, out-of-state-tree parsing.
 *
 * The HffFormat PluginBehavior registers a DataFormatProvider so HFF files
 * are recognized by drag-and-drop, plus a LoadHff state action for URL load.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { applyTransformerRaw, DataFormatProvider, rawDataObject } from '../../mol-plugin-state/formats/provider';
import { ShapeFormatCategory } from '../../mol-plugin-state/formats/shape';
import { PluginStateObject as SO, PluginStateTransform } from '../../mol-plugin-state/objects';
import { Download } from '../../mol-plugin-state/transforms/data';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { PluginContext } from '../../mol-plugin/context';
import { StateAction, StateObjectRef } from '../../mol-state';
import { Task } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { parseHff } from './parser';
import { SffData } from './schema';
import { Mat4 } from '../../mol-math/linear-algebra';
import { shapeFromSff } from './shape';


/** PluginStateObject wrapping parsed SFF data. */
export class SffDataObject extends SO.Create<SffData>({ name: 'SFF Data', typeClass: 'Data' }) { }


export type ParseHff = typeof ParseHff
export const ParseHff = PluginStateTransform.BuiltIn({
    name: 'parse-hff',
    display: { name: 'Parse HFF', description: 'Parse an EMDB-SFF segmentation file in its HDF5 (.hff) serialization.' },
    from: SO.Data.Binary,
    to: SffDataObject,
})({
    apply({ a }) {
        return Task.create('Parse HFF', async ctx => {
            const parsed = await parseHff(a.data).runInContext(ctx);
            if (parsed.isError) throw new Error(parsed.message);
            const sff = parsed.result;
            const segCount = sff.segments.length;
            const meshCount = sff.segments.reduce((acc, s) => acc + s.meshes.length, 0);
            return new SffDataObject(sff, {
                label: sff.name?.trim() || 'EMDB-SFF',
                description: `${segCount} segment${segCount === 1 ? '' : 's'}, ${meshCount} mesh${meshCount === 1 ? '' : 'es'}`,
            });
        });
    },
});


export type ShapeFromSff = typeof ShapeFromSff
export const ShapeFromSff = PluginStateTransform.BuiltIn({
    name: 'shape-from-sff',
    display: { name: 'Shape from SFF', description: 'Create a Shape from parsed EMDB-SFF data.' },
    from: SffDataObject,
    to: SO.Shape.Provider,
    params: {
        transforms: PD.Optional(PD.Value<Mat4[]>([Mat4.identity()], { isHidden: true })),
        label: PD.Optional(PD.Text('', { isHidden: true })),
    },
})({
    apply({ a, params }) {
        return Task.create('Shape from SFF', async ctx => {
            const provider = await shapeFromSff(a.data, params).runInContext(ctx);
            return new SO.Shape.Provider(provider, { label: provider.label, description: a.description });
        });
    },
});


export const HffProvider = DataFormatProvider({
    label: 'HFF',
    description: 'EMDB-SFF segmentation (HDF5)',
    category: ShapeFormatCategory,
    binaryExtensions: ['hff', 'h5', 'hdf5'],
    isApplicable: (info, data) => {
        // HDF5 magic bytes: \x89HDF\r\n\x1a\n at offset 0 (or in the user block,
        // but mol*'s built-in apps drop us here only by file extension anyway).
        if (data instanceof Uint8Array && data.length >= 8) {
            return data[0] === 0x89 && data[1] === 0x48 && data[2] === 0x44 && data[3] === 0x46
                && data[4] === 0x0D && data[5] === 0x0A && data[6] === 0x1A && data[7] === 0x0A;
        }
        return false;
    },
    parse: async (plugin, data) => {
        const format = plugin.state.data.build()
            .to(data)
            .apply(ParseHff, {}, { state: { isGhost: true } });

        const shape = format.apply(ShapeFromSff);

        await format.commit();

        return { format: format.selector, shape: shape.selector };
    },
    parseRaw: async (plugin, ctx, data) => {
        const format = await applyTransformerRaw(plugin, ctx, ParseHff, rawDataObject(data));
        const shape = await applyTransformerRaw(plugin, ctx, ShapeFromSff, format);
        return { shape: shape.data };
    },
    visuals(plugin: PluginContext, data: { shape: StateObjectRef<SO.Shape.Provider> }) {
        const repr = plugin.state.data.build()
            .to(data.shape)
            .apply(StateTransforms.Representation.ShapeRepresentation3D);
        return repr.commit();
    },
});


export const LoadHff = StateAction.build({
    display: { name: 'Load EMDB-SFF (HFF)', description: 'Load an EMDB-SFF segmentation HDF5 file from the specified URL.' },
    from: SO.Root,
    params: { url: PD.Text('') },
})(({ params, state }, ctx: PluginContext) => Task.create('Load HFF', taskCtx => {
    return state.transaction(async () => {
        const url = params.url.trim();
        if (url.length === 0) throw new Error('Specify URL');

        ctx.behaviors.layout.leftPanelTabName.next('data');

        const shapeNode = await state.build().toRoot()
            .apply(Download, { url, isBinary: true, label: `HFF: ${url}` })
            .apply(ParseHff, {}, { state: { isGhost: true } })
            .apply(ShapeFromSff)
            .commit();

        await state.build().to(shapeNode)
            .apply(StateTransforms.Representation.ShapeRepresentation3D)
            .commit();
    }).runInContext(taskCtx);
}));


export const HffFormat = PluginBehavior.create({
    name: 'hff-format',
    category: 'misc',
    display: {
        name: 'HFF (EMDB-SFF)',
        description: 'Read EMDB-SFF segmentation data in HDF5 format (.hff/.h5/.hdf5)',
    },
    ctor: class extends PluginBehavior.Handler {
        register(): void {
            this.ctx.state.data.actions.add(LoadHff);
            this.ctx.dataFormats.add('hff', HffProvider);
        }
        unregister(): void {
            this.ctx.state.data.actions.remove(LoadHff);
            this.ctx.dataFormats.remove('hff');
        }
    },
});
