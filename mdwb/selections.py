"""Shared residue-range shorthand for wizard, JSON and analysis planning."""
import re


def normalize_mask(value):
    """Interpret bare numbers/ranges as residues; preserve named/explicit masks."""
    if not isinstance(value, str) or not value.strip():
        raise ValueError('A selection must be nonempty text')
    value = value.strip()
    if re.fullmatch(r'[0-9,\-\s]+', value):
        compact = re.sub(r'\s*([,-])\s*', r'\1', value)
        for part in compact.split(','):
            match = re.fullmatch(r'([1-9]\d*)(?:-([1-9]\d*))?', part)
            if not match or (match[2] and int(match[2]) < int(match[1])):
                raise ValueError('Use positive residue numbers or ascending ranges, e.g. 6-8,12')
        return ':' + compact
    return value


def normalize_config_masks(config):
    """Normalize only selection fields, never names, filenames or frame ranges."""
    if 'fit_mask' in config:
        config['fit_mask'] = normalize_mask(config['fit_mask'])
    advanced = config.get('advanced', {})
    if 'matrix_mask' in advanced:
        advanced['matrix_mask'] = normalize_mask(advanced['matrix_mask'])
    imaging = config.get('imaging', {})
    if imaging.get('anchor'):
        imaging['anchor'] = normalize_mask(imaging['anchor'])
    for section in config.get('sections', []):
        section['mask'] = normalize_mask(section['mask'])
    for collection in ('interactions', 'solvation'):
        for item in config.get(collection, []):
            for key in ('mask1', 'mask2'):
                item[key] = normalize_mask(item[key])
    for monitor in config.get('monitors', []):
        monitor['masks'] = [normalize_mask(mask) for mask in monitor['masks']]
