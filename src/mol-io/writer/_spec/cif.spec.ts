import * as Data from '../../reader/cif/data-model';
import { CifWriter } from '../cif';
import decodeMsgPack from '../../common/msgpack/decode';
import { EncodedFile, EncodedCategory } from '../../common/binary-cif';
import Field from '../../reader/cif/binary/field';
import * as C from '../cif/encoder';

const cartn_x = Data.CifField.ofNumbers([1.001, 1.002, 1.003, 1.004, 1.005, 1.006, 1.007, 1.008, 1.009]);
const cartn_y = Data.CifField.ofNumbers([-3.0, -2.666, -2.3333, -2.0, -1.666, -1.333, -1.0, -0.666, -0.333]);
const cartn_z = Data.CifField.ofNumbers([1, 2, 3, 4, 5, 6, 7, 8, 9].map(i => Math.sqrt(i)));
const label_seq_id = Data.CifField.ofNumbers([1, 2, 3, 6, 11, 23, 47, 106, 235]);
const atom_site = Data.CifCategory.ofFields('atom_site', { 'Cartn_x': cartn_x, 'Cartn_y': cartn_y, 'Cartn_z': cartn_z, 'label_seq_id': label_seq_id });
const field1 = Data.CifField.ofNumbers([1, 2, 3, 6, 11, 23, 47, 106, 235]);
const field2 = Data.CifField.ofNumbers([-1, -2, -3, -6, -11, -23, -47, -106, -235]);
const other_fields = Data.CifCategory.ofFields('other_fields', { 'field1': field1, 'field2': field2 });

const encoding_aware_encoder = CifWriter.createEncoder({
    binary: true,
    binaryAutoClassifyEncoding: true,
    binaryEncodingPovider: CifWriter.createEncodingProviderFromJsonConfig([
        {
            'categoryName': 'atom_site',
            'columnName': 'Cartn_y',
            'encoding': 'rle',
            'precision': 0
        },
        {
            'categoryName': 'atom_site',
            'columnName': 'Cartn_z',
            'encoding': 'delta',
            'precision': 1
        },
        {
            'categoryName': 'atom_site',
            'columnName': 'label_seq_id',
            'encoding': 'delta-rle'
        }
    ])
});

describe('encoding-config', () => {
    const decoded = process(encoding_aware_encoder);

    const decoded_atom_site = decoded.blocks[0].categories['atom_site'];
    const decoded_cartn_x = decoded_atom_site.getField('Cartn_x')!;
    const decoded_cartn_y = decoded_atom_site.getField('Cartn_y')!;
    const decoded_cartn_z = decoded_atom_site.getField('Cartn_z')!;
    const decoded_label_seq_id = decoded_atom_site.getField('label_seq_id')!;

    const delta = 0.001;
    function assert(e: ArrayLike<number>, a: ArrayLike<number>) {
        expect(e.length).toBe(a.length);
        for (let i = 0; i < e.length; i++) {
            expect(Math.abs(e[i] - a[i])).toBeLessThan(delta);
        }
    }

    function join(field: Data.CifField) {
        return field.binaryEncoding!.map(e => e.kind).join();
    }

    it('strategy', () => {
        expect(join(decoded_cartn_x)).toBe('FixedPoint,Delta,IntegerPacking,ByteArray');
        expect(join(decoded_cartn_y)).toBe('FixedPoint,RunLength,IntegerPacking,ByteArray');
        expect(join(decoded_cartn_z)).toBe('FixedPoint,Delta,IntegerPacking,ByteArray');
        expect(join(decoded_label_seq_id)).toBe('Delta,RunLength,IntegerPacking,ByteArray');
    });

    it('precision', () => {
        assert(decoded_cartn_x.toFloatArray(), cartn_x.toFloatArray());
        assert(decoded_cartn_y.toFloatArray(), cartn_y.toFloatArray().map(d => Math.round(d)));
        assert(decoded_cartn_z.toFloatArray(), cartn_z.toFloatArray().map(d => Math.round(d * 10) / 10));
        assert(decoded_label_seq_id.toIntArray(), label_seq_id.toIntArray());
    });
});

const filter_aware_encoder1 = CifWriter.createEncoder({
    binary: true,
    binaryAutoClassifyEncoding: true
});
filter_aware_encoder1.setFilter(C.Category.filterOf('atom_site\n' +
'\n' +
'atom_site.Cartn_x\n' +
'atom_site.Cartn_y\n'));

const filter_aware_encoder2 = CifWriter.createEncoder({
    binary: true
});
filter_aware_encoder2.setFilter(C.Category.filterOf('!atom_site\n' +
'\n' +
'!other_fields.field2\n'));

describe('filtering-config', () => {
    const decoded1 = process(filter_aware_encoder1);

    const atom_site1 = decoded1.blocks[0].categories['atom_site'];
    const cartn_x1 = atom_site1.getField('Cartn_x');
    const cartn_y1 = atom_site1.getField('Cartn_y');
    const cartn_z1 = atom_site1.getField('Cartn_z');
    const label_seq_id1 = atom_site1.getField('label_seq_id');
    const fields1 = decoded1.blocks[0].categories['other_fields'];

    it('whitelist-filtering', () => {
        expect(atom_site1).toBeDefined();
        expect(cartn_x1).toBeDefined();
        expect(cartn_y1).toBeDefined();
        expect(cartn_z1).toBeUndefined();
        expect(label_seq_id1).toBeUndefined();
        expect(fields1).toBeUndefined();
    });

    const decoded2 = process(filter_aware_encoder2);

    const atom_site2 = decoded2.blocks[0].categories['atom_site'];
    const fields2 = decoded2.blocks[0].categories['other_fields'];
    const field12 = fields2.getField('field1');
    const field22 = fields2.getField('field2');

    it('blacklist-filtering', () => {
        expect(atom_site2).toBeUndefined();
        expect(fields2).toBeDefined();
        expect(field12).toBeDefined();
        expect(field22).toBeUndefined();
    });
});

function process(encoder: C.Encoder) {
    encoder.startDataBlock('test');

    for (const cat of [atom_site, other_fields]) {
        const fields: CifWriter.Field[] = [];
        for (const f of cat.fieldNames) {
            fields.push(wrap(f, cat.getField(f)!));
        }
        encoder.writeCategory(getCategoryInstanceProvider(cat, fields));
    }

    const encoded = encoder.getData() as Uint8Array;

    const unpacked = decodeMsgPack(encoded) as EncodedFile;
    return Data.CifFile(unpacked.dataBlocks.map(block => {
        const cats = Object.create(null);
        for (const cat of block.categories) cats[cat.name.substr(1)] = Category(cat);
        return Data.CifBlock(block.categories.map(c => c.name.substr(1)), cats, block.header);
    }));
}


function getCategoryInstanceProvider(cat: Data.CifCategory, fields: CifWriter.Field[]): CifWriter.Category {
    return {
        name: cat.name,
        instance: () => CifWriter.categoryInstance(fields, { data: cat, rowCount: cat.rowCount })
    };
}

function wrap(name: string, field: Data.CifField): CifWriter.Field {
    const type = Data.getCifFieldType(field);
    if (type['@type'] === 'str') {
        return { name, type: CifWriter.Field.Type.Str, value: field.str, valueKind: field.valueKind };
    } else if (type['@type'] === 'float') {
        return { name, type: CifWriter.Field.Type.Float, value: field.float, valueKind: field.valueKind };
    } else {
        return { name, type: CifWriter.Field.Type.Int, value: field.int, valueKind: field.valueKind };
    }
}

function Category(data: EncodedCategory): Data.CifCategory {
    const map = Object.create(null);
    const cache = Object.create(null);
    for (const col of data.columns) map[col.name] = col;
    return {
        rowCount: data.rowCount,
        name: data.name.substr(1),
        fieldNames: data.columns.map(c => c.name),
        getField(name) {
            const col = map[name];
            if (!col) return void 0;
            if (!!cache[name]) return cache[name];
            cache[name] = Field(col);
            return cache[name];
        }
    };
}