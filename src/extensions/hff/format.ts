/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Plugin integration for the HFF (EMDB-SFF / HDF5 segmentation) format.
 *
 * Mirrors the g3d extension's structure (see issue #886): a vendored HDF5
 * reader provides binary I/O, a typed parser produces SffData, and this
 * module wires it into mol*'s state-tree:
 *   Binary  -> ParseHff           -> SffDataObject
 *   SffData -> SffRepresentation3D -> Shape.Representation3D
 *
 * The Shape representation is a multi-visual one (see `representation.ts`)
 * with a `Visuals` MultiSelect to toggle Mesh and/or Wireframe — same UX as
 * the gaussian-surface representation.
 *
 * The HffFormat PluginBehavior registers a DataFormatProvider so HFF files
 * are recognized by drag-and-drop, plus a LoadHff state action for URL load.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { DataFormatProvider } from '../../mol-plugin-state/formats/provider';
import { ShapeFormatCategory } from '../../mol-plugin-state/formats/shape';
import { PluginStateObject as SO, PluginStateTransform } from '../../mol-plugin-state/objects';
import { Download } from '../../mol-plugin-state/transforms/data';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { PluginContext } from '../../mol-plugin/context';
import { StateAction, StateObjectRef, StateTransformer } from '../../mol-state';
import { Task } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { parseHff } from '../../mol-io/reader/hff/parser';
import { SffData } from '../../mol-io/reader/hff/schema';
import { SffRepresentation, SffRepresentationParams } from './representation';


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


export type SffRepresentation3D = typeof SffRepresentation3D
export const SffRepresentation3D = PluginStateTransform.BuiltIn({
    name: 'sff-representation-3d',
    display: { name: 'SFF Representation', description: 'Multi-visual representation of an EMDB-SFF segmentation (Mesh and/or Wireframe).' },
    from: SffDataObject,
    to: SO.Shape.Representation3D,
    params: () => SffRepresentationParams,
})({
    canAutoUpdate() {
        return true;
    },
    apply({ a, params }, plugin: PluginContext) {
        return Task.create('SFF Representation', async ctx => {
            const repr = SffRepresentation(
                {
                    webgl: plugin.canvas3d?.webgl,
                    colorThemeRegistry: plugin.representation.structure.themes.colorThemeRegistry,
                    sizeThemeRegistry: plugin.representation.structure.themes.sizeThemeRegistry,
                },
                () => SffRepresentationParams,
            );
            const props = { ...PD.getDefaultValues(SffRepresentationParams), ...params };
            await repr.createOrUpdate(props, a.data).runInContext(ctx);
            return new SO.Shape.Representation3D({ repr, sourceData: a.data }, { label: a.label, description: a.description });
        });
    },
    update({ a, b, newParams }) {
        return Task.create('SFF Representation', async ctx => {
            const props = { ...b.data.repr.props, ...newParams };
            await b.data.repr.createOrUpdate(props, a.data).runInContext(ctx);
            b.data.sourceData = a.data;
            return StateTransformer.UpdateResult.Updated;
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
        const sff = plugin.state.data.build()
            .to(data)
            .apply(ParseHff, {}, { state: { isGhost: true } });
        await sff.commit();
        return { format: sff.selector, shape: sff.selector };
    },
    visuals(plugin: PluginContext, data: { shape: StateObjectRef<SffDataObject> }) {
        return plugin.state.data.build()
            .to(data.shape)
            .apply(SffRepresentation3D)
            .commit();
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

        const sffNode = await state.build().toRoot()
            .apply(Download, { url, isBinary: true, label: `HFF: ${url}` })
            .apply(ParseHff, {}, { state: { isGhost: true } })
            .commit();

        await state.build().to(sffNode)
            .apply(SffRepresentation3D)
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
