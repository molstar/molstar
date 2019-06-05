import './index.html'
import { CIF, CifCategory, CifField, getCifFieldType } from '../../mol-io/reader/cif';
import { CifWriter } from '../../mol-io/writer/cif';

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
    const parsed = await downloadCif(`https://webchem.ncbr.muni.cz/ModelServer/static/bcif/${pdb}`, true);
    return parsed.blocks[0];
}

async function init(props = {}) {
    const cif = await downloadFromPdb('1brr')
    const encoder = CifWriter.createEncoder({
        binary: true,
        encoderName: 'mol*',
        binaryAutoClassifyEncoding: true,
        binaryEncodingPovider: CifWriter.createEncodingProviderFromJsonConfig([
            // {
            //     'categoryName': 'atom_site',
            //     'columnName': 'Cartn_x',
            //     'encoding': 'delta',
            //     'precision': 3
            // },
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

    encoder.startDataBlock(cif.header);
    for (const c of cif.categoryNames) {
        const cat = cif.categories[c];
        const fields: CifWriter.Field[] = [];
        for (const f of cat.fieldNames) {
            fields.push(wrap(f, cat.getField(f)!))
        }

        encoder.writeCategory(getCategoryInstanceProvider(cif.categories[c], fields));
    }
    const ret = encoder.getData() as Uint8Array;

    const cif2 = (await parseCif(ret)).blocks[0];
    // should be untouched: delta encoding
    console.log(cif2.categories['atom_site'].getField('Cartn_x'));
    // should have rle encoding, 0 decimal places
    console.log(cif2.categories['atom_site'].getField('Cartn_y'));
    // should have delta encoding, 1 decimal place
    console.log(cif2.categories['atom_site'].getField('Cartn_z'));
    // should use delta-rle encoding
    console.log(cif2.categories['atom_site'].getField('label_seq_id'));
}

init()

function getCategoryInstanceProvider(cat: CifCategory, fields: CifWriter.Field[]): CifWriter.Category {
    return {
        name: cat.name,
        instance: () => CifWriter.categoryInstance(fields, { data: cat, rowCount: cat.rowCount })
    };
}

function wrap(name: string, field: CifField): CifWriter.Field {
    const type = getCifFieldType(field);
    if (type['@type'] === 'str') {
        return { name, type: CifWriter.Field.Type.Str, value: field.str, valueKind: field.valueKind };
    } else if (type['@type'] === 'float') {
        return { name, type: CifWriter.Field.Type.Float, value: field.float, valueKind: field.valueKind };
    } else {
        return { name, type: CifWriter.Field.Type.Int, value: field.int, valueKind: field.valueKind };
    }
}