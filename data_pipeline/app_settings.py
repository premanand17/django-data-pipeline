
# elastic indices used only by pipeline
ELASTIC = {
    'BAND': {
        'name': 'bands_hg38',
        'idx_type': {
            'BAND': {'type': 'bands'},
            'CHROM': {'type': 'chromosome'}
        }
    },
    'HAPMAP': {
        'name': 'hapmap_phase2_b38',
        'idx_type': {
            'HAPMAP': {'type': 'hapmap', 'auth_public': True}
        },
        'auth_public': True
    }
}
