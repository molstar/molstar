import { CifWriter } from '../cif';

const cif = await downloadFromPdb('1brr')
const encoder = CifWriter.createEncoder({
    binary: true,
    encoderName: 'mol*',
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