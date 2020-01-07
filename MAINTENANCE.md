# HIVFacts Data Maintenance

## Dependencies

1. Read access to HIVDB2 database (default: `127.0.0.1:3308/HIVDB2`). Can be
   changed by specifying environment parameter `DATABASE_URI`. Check
   `hivdbql/config.py` (in private repo) for more information.
2. Read/write access to HIVDB\_Rules database (default:
   `10.77.6.244:3306/HIVDB_Rules`). Can be changed by specifying environment
   parameter `DATABASE_URI_HIVRULES`.
3. Python 3.7 runtime + pipenv. Run this command at the root directory:
   `pipenv install`.

## `data/aapcnt`

To update this folder, run this command:

```bash
scripts/update-aapcnt.sh
```

### `data/codonpcnt`

To update this folder, run this command:

```bash
scripts/update-codonpcnt.sh
```

### `data/apobecs`

To update this folder, run following commands:

```bash
scripts/prepare-apobec-data.sh
scripts/update-apobec.sh
```

### `data/drm_hiv1.json`
List of HIV-1 drug resistance mutation based on latest HIVDB algorithm.

```bash
scripts/update-drmlist.sh
```

### `data/sdrm_hiv1.json`
List of HIV-1 drug resistance mutation based on latest HIVDB algorithm.

```bash
scripts/update-sdrmlist.sh
```

### HIVDB algorithm in `data/algorithms`

TODO
