# cellAdmix admixture audit — validation on the pancreas Xenium crop

Lower total admixture = cleaner segmentation.

| variant | typing | seed | total_admixture_rate | admixed molecules | detected/evaluated pairs | cells | molecules used | runtime (s) |
|---|---|---|---|---|---|---|---|---|
| vendor | quick_cluster | 1 | 0.075851 | 9103 | 9/38 | 2863 | 120011 | 4.7 |
| vendor_seed2 | quick_cluster | 2 | 0.075851 | 9103 | 9/38 | 2863 | 120011 | 4.7 |
| vendor_repeat | quick_cluster | 1 | 0.075851 | 9103 | 9/38 | 2863 | 120011 | 4.7 |
| border10 | celltypes | 1 | 0.089800 | 10777 | 10/39 | 2863 | 120011 | 4.1 |
| border30 | celltypes | 1 | 0.094625 | 11356 | 10/39 | 2863 | 120011 | 4.1 |
| dilate2 | celltypes | 1 | 0.143193 | 19485 | 10/39 | 2862 | 136075 | 4.2 |
| border30_recluster | quick_cluster | 1 | 0.117148 | 14059 | 8/37 | 2863 | 120011 | 4.6 |
| baysor | celltypes | 1 | 0.048090 | 6748 | 5/27 | 1762 | 140320 | 3.9 |

## Checks

- **border_degradation_monotone**: PASS — `{"chain": ["vendor", "border10", "border30"], "checked": true, "values": {"vendor": 0.07585138029014007, "border10": 0.08980010165734807, "border30": 0.0946246594062211}, "steps": [{"from": "vendor", "to": "border10", "rate_from": 0.07585138029014007, "rate_to": 0.08980010165734807, "increased": true}, {"from": "border10", "to": "border30", "rate_from": 0.08980010165734807, "rate_to": 0.0946246594062211, "increased": true}], "passed": true}`
- **dilation_worse_than_vendor**: PASS — `{"chain": ["vendor", "dilate2"], "checked": true, "values": {"vendor": 0.07585138029014007, "dilate2": 0.14319309204482822}, "steps": [{"from": "vendor", "to": "dilate2", "rate_from": 0.07585138029014007, "rate_to": 0.14319309204482822, "increased": true}], "passed": true}`
- **seed_stochasticity**: info — `{"rate_seed1": 0.07585138029014007, "rate_seed2": 0.07585138029014007, "abs_difference": 0.0, "relative_difference": 0.0, "note": "two full runs with different --seed (the seed reaches the NMF fit; the audit is factorization-independent) plus a fresh quick-clustering draw; this sets the comparison tolerance"}`
- **same_seed_deterministic**: info — `{"abs_difference": 0.0}`
- **transfer_vs_recluster_typing**: info — `{"transfer_rate": 0.0946246594062211, "recluster_rate": 0.11714759480380965, "vendor_rate": 0.07585138029014007}`
