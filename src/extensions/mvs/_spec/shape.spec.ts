/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { parseVtp } from '../../../mol-io/reader/vtp/parser';
import { shapeFromVtp } from '../../../mol-model-formats/shape/vtp';
import { Task } from '../../../mol-task';
import { Color } from '../../../mol-util/color';
import { ColorNames } from '../../../mol-util/color/names';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { shapeRepresentationProps } from '../helpers/load-shape';
import { MVSData } from '../mvs-data';
import { convertMvsToMolstar } from '../tree/molstar/conversion';
import { MVSTree } from '../tree/mvs/mvs-tree';

/** Octahedron of radius 12 with a per-vertex scalar, in ASCII VTP.
 * Deliberately at molecular scale so that "renders but the camera is wrong" is
 * distinguishable from "does not render". */
const OCTAHEDRON_VTP = `<?xml version="1.0"?>
<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian">
  <PolyData>
    <Piece NumberOfPoints="6" NumberOfPolys="8">
      <PointData Scalars="height">
        <DataArray type="Float32" Name="height" NumberOfComponents="1" format="ascii">
          0 0 0 0 12 -12
        </DataArray>
      </PointData>
      <Points>
        <DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii">
          12 0 0  -12 0 0  0 12 0  0 -12 0  0 0 12  0 0 -12
        </DataArray>
      </Points>
      <Polys>
        <DataArray type="Int32" Name="connectivity" format="ascii">
          0 2 4  2 1 4  1 3 4  3 0 4  2 0 5  1 2 5  3 1 5  0 3 5
        </DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">
          3 6 9 12 15 18 21 24
        </DataArray>
      </Polys>
    </Piece>
  </PolyData>
</VTKFile>
`;

describe('MVS shape node', () => {
    it('builder produces a download -> parse -> shape tree', () => {
        const builder = MVSData.createBuilder();
        builder.download({ url: 'octahedron.vtp' })
            .parse({ format: 'vtp' })
            .shape().color({ color: '#3b82f6' });
        const state = builder.getState();

        expect(MVSData.validationIssues(state)).toEqual(undefined);

        const parse = state.root.children![0].children![0];
        expect(parse.kind).toEqual('parse');
        expect(parse.params).toEqual({ format: 'vtp' });

        const shape = parse.children![0];
        expect(shape.kind).toEqual('shape');
        expect(shape.params).toEqual({});
        expect(shape.children![0].kind).toEqual('color');
        expect(shape.children![0].params).toEqual({ color: '#3b82f6' });
    });

    it('shape accepts vtp coloring options as custom props', () => {
        const builder = MVSData.createBuilder();
        builder.download({ url: 'surface.vtp' })
            .parse({ format: 'vtp' })
            .shape({}, { vtp_attribute: 'height', vtp_attribute_source: 'point', vtp_palette: { kind: 'continuous', colors: 'Viridis' } });
        expect(MVSData.validationIssues(builder.getState())).toEqual(undefined);
    });

    it('all three shape formats are valid parse formats', () => {
        for (const format of ['vtp', 'ply', 'obj'] as const) {
            const builder = MVSData.createBuilder();
            builder.download({ url: `mesh.${format}` }).parse({ format }).shape({});
            expect(MVSData.validationIssues(builder.getState())).toEqual(undefined);
        }
    });

    it('converts to a Molstar tree, marking vtp/ply as binary and obj as text', () => {
        for (const [format, isBinary] of [['vtp', true], ['ply', true], ['obj', false]] as const) {
            const builder = MVSData.createBuilder();
            builder.download({ url: `mesh.${format}` }).parse({ format }).shape({});
            const converted = convertMvsToMolstar(builder.getState().root as MVSTree, undefined);
            // `download` and `parse` are condensed into one node carrying both params
            const download = converted.children![0];
            expect(download.params).toMatchObject({ url: `mesh.${format}`, is_binary: isBinary });
        }
    });

    it('rejects a shape node that is not under parse', () => {
        const tree: MVSTree = {
            kind: 'root',
            children: [{ kind: 'shape', params: {} } as any],
        };
        expect(MVSData.validationIssues({ kind: 'single', root: tree, metadata: { version: '1', timestamp: '' } })).not.toEqual(undefined);
    });

    it('feeds the VTP shape provider explicit params, never its file-derived defaults', async () => {
        const parsed = await parseVtp(new TextEncoder().encode(OCTAHEDRON_VTP)).run();
        if (parsed.isError) throw new Error(parsed.message);
        const file = parsed.result;

        // The provider is molstar's own, i.e. identical to opening the file through the UI.
        const provider = await shapeFromVtp(file).run();
        expect(Object.keys(provider.params)).toContain('colorTheme');
        expect(Object.keys(provider.params)).toContain('scale');

        // This file has a point array, so the provider's own default is attribute coloring.
        // MVS must not inherit that -- a node with no `attribute` has to render uniform.
        expect((provider.params.colorTheme as any).defaultValue.name).toEqual('attribute');

        const node = { kind: 'shape', params: {}, children: [{ kind: 'color', params: { color: '#3b82f6' } }] } as any;
        const props = shapeRepresentationProps(node, 'vtp') as any;
        expect(props.colorTheme.name).toEqual('uniform');
        expect(props.colorTheme.params.color).toEqual(Color.fromHexStyle('#3b82f6'));

        // And with an attribute set, MVS selects it explicitly using the provider's own key format.
        const attrNode = { kind: 'shape', params: {}, custom: { vtp_attribute: 'height', vtp_attribute_source: 'point' } } as any;
        const attrProps = shapeRepresentationProps(attrNode, 'vtp') as any;
        expect(attrProps.colorTheme.name).toEqual('attribute');
        expect(attrProps.colorTheme.params.name).toEqual('point:height');
        expect(attrProps.colorTheme.params.colors).toEqual('viridis');
    });

    it('renders the same mesh and coloring as the shape provider does directly', async () => {
        const parsed = await parseVtp(new TextEncoder().encode(OCTAHEDRON_VTP)).run();
        if (parsed.isError) throw new Error(parsed.message);
        const provider = await shapeFromVtp(parsed.result).run();

        const node = { kind: 'shape', params: {}, custom: { vtp_attribute: 'height', vtp_attribute_source: 'point' } } as any;
        const props = { ...PD.getDefaultValues(provider.params), ...shapeRepresentationProps(node, 'vtp') } as any;
        const shape = await Task.create('t', ctx => provider.getShape(ctx, provider.data, props)).run();

        expect(shape.geometry.vertexCount).toEqual(6);
        expect(shape.geometry.triangleCount).toEqual(8);
        // `height` is [0,0,0,0,12,-12]: min and max vertices must not share a color.
        expect(shape.getColor(5, 0)).not.toEqual(shape.getColor(4, 0));
    });

    it('passes a palette name and domain straight through to the provider', () => {
        const props = (custom: any) => shapeRepresentationProps({ kind: 'shape', params: {}, custom } as any, 'vtp') as any;

        // Setting colorTheme replaces the whole mapped value, so `colors` is always stated --
        // omitting it would leave the param undefined rather than defaulted.
        expect(props({ vtp_attribute: 'height' }).colorTheme.params.colors).toEqual('viridis');

        // A molstar color list name is legitimate in a custom prop, which is a vendor extension.
        expect(props({ vtp_attribute: 'height', vtp_palette: 'viridis' }).colorTheme.params.colors).toEqual('viridis');

        // Domain: absent -> auto, given -> a custom interval.
        expect(props({ vtp_attribute: 'height' }).colorTheme.params.domain).toEqual({ name: 'auto', params: { symmetric: false } });
        expect(props({ vtp_attribute: 'height', vtp_domain: [0, 100] }).colorTheme.params.domain).toEqual({ name: 'custom', params: [0, 100] });
    });

    it('ply_coloring / obj_coloring select where colors come from', () => {
        const props = (custom: any, format: any) => shapeRepresentationProps({ kind: 'shape', params: {}, custom } as any, format) as any;

        // Default is uniform for both, i.e. the file's own colors are not used unless asked for.
        expect(props({}, 'ply').coloring.name).toEqual('uniform');
        expect(props({}, 'obj').coloring.name).toEqual('uniform');

        // PLY vertex/material coloring pins the conventional property names rather than relying
        // on whatever the provider auto-detects.
        const plyVertex = props({ ply_coloring: 'vertex' }, 'ply');
        expect(plyVertex.coloring.name).toEqual('vertex');
        expect(plyVertex.coloring.params).toEqual({ red: 'red', green: 'green', blue: 'blue' });
        expect(props({ ply_coloring: 'material' }, 'ply').coloring.name).toEqual('material');

        // OBJ vertex colors are embedded, so no params are needed.
        expect(props({ obj_coloring: 'vertex' }, 'obj').coloring).toEqual({ name: 'vertex', params: {} });

        // OBJ 'custom' colors each material group by name -- no MTL file involved, since the names
        // come from the OBJ's own usemtl directives. Colors are stated so the result is reproducible.
        const objCustom = props({ obj_coloring: 'custom', obj_material_colors: { shell: '#ff0000', core: 'blue' } }, 'obj');
        expect(objCustom.coloring.name).toEqual('custom');
        expect(objCustom.coloring.params).toEqual({ shell: Color.fromHexStyle('#ff0000'), core: ColorNames.blue });
    });

    it('uses each format own coloring param shape', () => {
        const node = { kind: 'shape', params: {}, children: [{ kind: 'color', params: { color: 'red' } }] } as any;
        expect((shapeRepresentationProps(node, 'vtp') as any).colorTheme.name).toEqual('uniform');
        expect((shapeRepresentationProps(node, 'ply') as any).coloring.name).toEqual('uniform');
        expect((shapeRepresentationProps(node, 'obj') as any).coloring.name).toEqual('uniform');
    });

    it('parses the ASCII octahedron fixture and reads its point attribute', async () => {
        const bytes = new TextEncoder().encode(OCTAHEDRON_VTP);
        const parsed = await parseVtp(bytes).run();
        if (parsed.isError) throw new Error(parsed.message);
        const file = parsed.result;

        expect(file.numberOfPoints).toEqual(6);
        expect(file.numberOfTriangles).toEqual(8);

        const height = file.pointData.get('height');
        expect(height?.numberOfComponents).toEqual(1);
        expect(Array.from(height!.values.toArray())).toEqual([0, 0, 0, 0, 12, -12]);
    });
});
