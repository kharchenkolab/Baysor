# cellAdmix audit on xenium_pancreas_377_full (n_pool=60, fixed typing, fixed pairs)

Lower total admixture = cleaner segmentation.

| variant | typing | seed | total_admixture_rate | admixed molecules | detected/evaluated pairs | cells | molecules used | runtime (s) |
|---|---|---|---|---|---|---|---|---|
| vendor | celltypes | 1 | 0.155925 | 253371 | 20/20 | 39544 | 1624955 | 15.0 |
| border10 | celltypes | 1 | 0.157050 | 255199 | 19/20 | 39544 | 1624955 | 13.1 |
| border30 | celltypes | 1 | 0.181670 | 295205 | 19/20 | 39543 | 1624955 | 13.4 |
| dilate2 | celltypes | 1 | 0.216188 | 398072 | 19/20 | 39534 | 1841326 | 13.9 |

## Checks

- **border_degradation_monotone**: PASS — `{"chain": ["vendor", "border10", "border30"], "checked": true, "values": {"vendor": 0.1559249333058454, "border10": 0.15704988753534713, "border30": 0.1816696462363573}, "steps": [{"from": "vendor", "to": "border10", "rate_from": 0.1559249333058454, "rate_to": 0.15704988753534713, "increased": true}, {"from": "border10", "to": "border30", "rate_from": 0.15704988753534713, "rate_to": 0.1816696462363573, "increased": true}], "passed": true}`
- **dilation_worse_than_vendor**: PASS — `{"chain": ["vendor", "dilate2"], "checked": true, "values": {"vendor": 0.1559249333058454, "dilate2": 0.21618768213776376}, "steps": [{"from": "vendor", "to": "dilate2", "rate_from": 0.1559249333058454, "rate_to": 0.21618768213776376, "increased": true}], "passed": true}`
