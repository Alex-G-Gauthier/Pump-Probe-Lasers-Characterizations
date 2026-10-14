from pathlib import Path
import importlib
import json
import pandas as pd

ROOT = Path(__file__).resolve().parent
required = [
    'laser_characterization.py', 'laser_characterization_core.py', 'experiment_metadata.json',
    'outputs/data_inventory.csv', 'outputs/spatial_frame_metrics.csv',
    'outputs/temporal_record_metrics.csv', 'outputs/spectral_metrics.csv',
    'outputs/integrated_results.json', 'outputs/power_energy_summary.csv',
]
missing = [x for x in required if not (ROOT / x).exists()]
if missing:
    raise SystemExit('Missing required files:\n  ' + '\n  '.join(missing))
importlib.import_module('laser_characterization_core')
importlib.import_module('laser_characterization')
for csv in ['spatial_frame_metrics.csv','temporal_record_metrics.csv','spectral_metrics.csv','power_energy_summary.csv']:
    df = pd.read_csv(ROOT/'outputs'/csv)
    if df.empty:
        raise SystemExit(f'{csv} is unexpectedly empty')
print('PASS: helper modules import and all required precomputed outputs are present.')
