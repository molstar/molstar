/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SimulariumGeometryResolver } from '../../mol-model-formats/particles/simularium';
import { ParticleList, ParticleTarget } from '../../mol-model/particles/particle-list';
import { Structure } from '../../mol-model/structure';
import { PluginConfig } from '../../mol-plugin/config';
import { PluginContext } from '../../mol-plugin/context';
import { RuntimeContext, Task } from '../../mol-task';
import { Asset } from '../../mol-util/assets';
import { getFileNameInfo } from '../../mol-util/file-info';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { rawDataObject } from '../formats/provider';

/** The registered formats that can be loaded as a reference object for a particle target. */
function targetFormatProviders(plugin: PluginContext) {
    return plugin.dataFormats.list.filter(e => !!e.provider.parseRaw);
}

/** Options for a `PD.Select` of the supported target formats, including auto-detection. */
export function particleTargetFormatOptions(plugin: PluginContext): [string, string][] {
    return [
        ['auto', 'Auto'],
        ...targetFormatProviders(plugin).map(e => [e.name, e.provider.label] as [string, string]),
    ];
}

/** Comma separated extension list for a `PD.File` accept attribute. */
export function particleTargetFileExtensions(plugin: PluginContext): string {
    const extensions = new Set<string>();
    for (const { provider } of targetFormatProviders(plugin)) {
        provider.stringExtensions?.forEach(e => extensions.add(e));
        provider.binaryExtensions?.forEach(e => extensions.add(e));
    }
    return Array.from(extensions).map(e => `.${e}`).join(',');
}

export type ParticleTargetSource =
    | { readonly kind: 'url', readonly url: string | Asset.Url, readonly format?: string, readonly isBinary?: boolean }
    | { readonly kind: 'file', readonly file: Asset.File, readonly format?: string }

export interface ParticleTargetEntry {
    readonly targetId: number
    readonly target: ParticleTarget
}

export interface MatchedParticleTargetFile {
    readonly file: Asset.File
    readonly targetIds: ReadonlyArray<number>
}

export interface ParticleTargetFileMatches {
    readonly matches: ReadonlyArray<MatchedParticleTargetFile>
    readonly warnings: ReadonlyArray<string>
}

/** Match files by basename to the entity name associated with each particle target id. */
export function matchParticleTargetFiles(particles: ParticleList, files: ReadonlyArray<Asset.File>, excludedTargetIds?: ReadonlySet<number>): ParticleTargetFileMatches {
    const warnings: string[] = [];
    const { entities, entityInfo, targets } = particles;
    if (!entities || !entityInfo) {
        if (files.length > 0) warnings.push('Cannot match particle target files because the particle list has no entity metadata.');
        return { matches: [], warnings };
    }

    const namesByTarget = new Map<number, Set<string>>();
    const invalidTargets = new Set<number>();
    for (let i = 0; i < particles.count; ++i) {
        const targetId = targets[i];
        const info = entityInfo.get(entities[i]);
        if (!info) {
            invalidTargets.add(targetId);
            continue;
        }
        let names = namesByTarget.get(targetId);
        if (!names) {
            names = new Set<string>();
            namesByTarget.set(targetId, names);
        }
        names.add(info.name);
    }

    const targetIdsByName = new Map<string, number[]>();
    const excludedNames = new Set<string>();
    for (const [targetId, names] of namesByTarget) {
        if (invalidTargets.has(targetId) || names.size !== 1) {
            if (!excludedTargetIds?.has(targetId)) warnings.push(`Cannot match particle target ${targetId} because it does not map to exactly one entity name.`);
            continue;
        }
        const name = names.values().next().value!;
        if (excludedTargetIds?.has(targetId)) {
            excludedNames.add(name);
            continue;
        }
        let targetIds = targetIdsByName.get(name);
        if (!targetIds) {
            targetIds = [];
            targetIdsByName.set(name, targetIds);
        }
        targetIds.push(targetId);
    }
    for (const targetId of invalidTargets) {
        if (!excludedTargetIds?.has(targetId) && !namesByTarget.has(targetId)) warnings.push(`Cannot match particle target ${targetId} because it has no entity name.`);
    }

    const filesByBase = new Map<string, Asset.File[]>();
    for (const file of files) {
        const base = getFileNameInfo(file.file?.name ?? file.name).base;
        let entries = filesByBase.get(base);
        if (!entries) {
            entries = [];
            filesByBase.set(base, entries);
        }
        entries.push(file);
    }

    const matches: MatchedParticleTargetFile[] = [];
    for (const [base, entries] of filesByBase) {
        if (entries.length !== 1) {
            warnings.push(`Cannot match particle target files named '${base}' because the filename is not unique.`);
            continue;
        }
        const targetIds = targetIdsByName.get(base);
        if (!targetIds) {
            if (excludedNames.has(base)) continue;
            warnings.push(`Could not match particle target file '${entries[0].name}' to an entity name.`);
            continue;
        }
        matches.push({ file: entries[0], targetIds });
    }
    return { matches, warnings };
}

/** Turn a parsed data object into the reference object of a particle target. */
async function particleTargetFromParsed(plugin: PluginContext, ctx: RuntimeContext, parsed: any, name: string): Promise<ParticleTarget> {
    const volume = parsed.volume ?? parsed.volumes?.[0];
    if (volume) return { kind: 'volume', volume };

    if (parsed.shape) {
        const provider = parsed.shape;
        const shape = await provider.getShape(ctx, provider.data, PD.getDefaultValues(provider.params));
        return { kind: 'shape', shape };
    }

    if (parsed.structure) {
        return { kind: 'structure', structure: parsed.structure };
    }

    if (parsed.trajectory) {
        const model = await Task.resolveInContext(parsed.trajectory.getFrameAtIndex(0), ctx);
        return { kind: 'structure', structure: Structure.ofModel(model) };
    }

    throw new Error(`'${name}' does not provide a structure, shape, or volume.`);
}

/**
 * Load a reference object (structure, shape, or volume) for a particle target from a URL or
 * a file using the registered data format providers. The object is meant to be attached to a
 * `ParticleList` via its `targetMapping`.
 *
 * The resolved asset is kept in the asset manager (so it is part of session snapshots) and
 * added to `assets`, whose owner must dispose the wrappers. Without `assets` it is released
 * again once parsed.
 */
export async function loadParticleTarget(plugin: PluginContext, ctx: RuntimeContext, source: ParticleTargetSource, assets?: Asset.Wrapper[]): Promise<ParticleTarget> {
    const asset = source.kind === 'url'
        ? Asset.getUrlAsset(plugin.managers.asset, source.url)
        : source.file;
    const info = getFileNameInfo(source.kind === 'url'
        ? Asset.getUrl(source.url)
        : source.file.file?.name ?? source.file.name);
    const isBinary = (source.kind === 'url' ? source.isBinary : undefined) ?? plugin.dataFormats.binaryExtensions.has(info.ext);

    const wrapper = await plugin.managers.asset.resolve(asset, isBinary ? 'binary' : 'string').runInContext(ctx);
    let keepAsset = false;
    try {
        const data = wrapper.data;
        const provider = !source.format || source.format === 'auto'
            ? plugin.dataFormats.auto(info, rawDataObject(data))
            : plugin.dataFormats.get(source.format);
        if (!provider) throw new Error(`Could not find a data format provider for '${info.name}'.`);
        if (!provider.parseRaw) throw new Error(`The '${provider.label}' format cannot be used as a particle target.`);

        const target = await particleTargetFromParsed(plugin, ctx, await provider.parseRaw(plugin, ctx, data), info.name);
        if (assets) {
            assets.push(wrapper);
            keepAsset = true;
        }
        return target;
    } finally {
        if (!keepAsset) wrapper.dispose();
    }
}

/** Load the reference objects of several targets, skipping the ones that fail with a warning. */
export async function loadParticleTargets(plugin: PluginContext, ctx: RuntimeContext, sources: ReadonlyArray<{ targetId: number, source: ParticleTargetSource }>, assets?: Asset.Wrapper[]): Promise<ParticleTargetEntry[]> {
    const entries: ParticleTargetEntry[] = [];
    for (const { targetId, source } of sources) {
        try {
            entries.push({ targetId, target: await loadParticleTarget(plugin, ctx, source, assets) });
        } catch (e) {
            console.error(e);
            plugin.log.warn(`Could not load the reference object for particle target ${targetId}.`);
        }
    }
    return entries;
}

/**
 * Resolve a value that is either a 4-character PDB id or a direct URL into a target source,
 * mirroring the `DownloadStructure` action: PDB ids are fetched from the default provider.
 */
export function particleTargetSourceFromPdbIdOrUrl(plugin: PluginContext, idOrUrl: string): ParticleTargetSource {
    const id = idOrUrl.trim();
    if (/^[1-9][a-z0-9]{3}$/i.test(id)) {
        const provider = plugin.config.get(PluginConfig.Download.DefaultPdbProvider) || 'pdbe';
        const url = provider === 'rcsb'
            ? `https://models.rcsb.org/${id.toUpperCase()}.bcif`
            : `https://www.ebi.ac.uk/pdbe/entry-files/download/${id.toLowerCase()}.bcif`;
        return { kind: 'url', url, format: 'mmcif', isBinary: true };
    }
    return { kind: 'url', url: id };
}

/** Resolves the geometries referenced by Simularium agent types into reference objects. */
export function createSimulariumGeometryResolver(plugin: PluginContext, ctx: RuntimeContext, assets?: Asset.Wrapper[]): SimulariumGeometryResolver {
    return async geometry => {
        const source: ParticleTargetSource = geometry.displayType === 'PDB'
            ? particleTargetSourceFromPdbIdOrUrl(plugin, geometry.url)
            : { kind: 'url', url: geometry.url };
        try {
            return await loadParticleTarget(plugin, ctx, source, assets);
        } catch (e) {
            console.error(e);
            plugin.log.warn(`Could not load the geometry of the Simularium agent type '${geometry.name}'.`);
            return undefined;
        }
    };
}
