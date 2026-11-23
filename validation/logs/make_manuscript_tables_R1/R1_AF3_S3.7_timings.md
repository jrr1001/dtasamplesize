Timings depend on the machine and on the load at the time of the run.

| scenario | method | B | grid | elapsed_s | note |
|---|---|---|---|---|---|
| P1 | exact (own enumeration up to N = 900) |  | N = 2..900 | 7.2144 | script's own qbeta enumeration, no legacy fields |
| P1 | exact (bam_sample_size defaults) | 5000 | automatic integer scan | 26.4200 | includes legacy per-arm Monte Carlo fields (B = 5000, n_range 20:500) |
| P1 | exact (bam_sample_size, B = 1000, n_range = 20:60) | 1000 | automatic integer scan | 9.7400 | minimal legacy fields; the headline exact result does not depend on B |
| P2 | exact (own enumeration up to N = 900) |  | N = 2..900 | 7.2522 | script's own qbeta enumeration, no legacy fields |
| P2 | exact (bam_sample_size defaults) | 5000 | automatic integer scan | 19.3700 | includes legacy per-arm Monte Carlo fields (B = 5000, n_range 20:500) |
| P2 | exact (bam_sample_size, B = 1000, n_range = 20:60) | 1000 | automatic integer scan | 9.2900 | minimal legacy fields; the headline exact result does not depend on B |
| P1 | monte_carlo automatic grid, seed 1 | 2000 | 372:1170 | 17.8700 | bam_sample_size defaults otherwise (includes legacy per-arm fields at this B) |
| P1 | monte_carlo explicit window, seed 1 | 2000 | 618:838 | 10.1900 | window around N*; includes legacy per-arm fields at this B |
| P1 | monte_carlo automatic grid, seed 1 | 10000 | 372:1125 | 82.4400 | bam_sample_size defaults otherwise (includes legacy per-arm fields at this B) |
| P1 | monte_carlo explicit window, seed 1 | 10000 | 618:838 | 42.9200 | window around N*; includes legacy per-arm fields at this B |
| P2 | monte_carlo automatic grid, seed 1 | 2000 | 171:1170 | 20.3300 | bam_sample_size defaults otherwise (includes legacy per-arm fields at this B) |
| P2 | monte_carlo explicit window, seed 1 | 2000 | 612:832 | 6.7800 | window around N*; includes legacy per-arm fields at this B |
| P2 | monte_carlo automatic grid, seed 1 | 10000 | 168:1125 | 86.2800 | bam_sample_size defaults otherwise (includes legacy per-arm fields at this B) |
| P2 | monte_carlo explicit window, seed 1 | 10000 | 612:832 | 24.3600 | window around N*; includes legacy per-arm fields at this B |
| P1 | monte_carlo explicit window, mean over declared seeds | 2000 |  | 8.3100 | mean elapsed per call |
| P2 | monte_carlo explicit window, mean over declared seeds | 2000 |  | 5.3040 | mean elapsed per call |
| P1 | monte_carlo explicit window, mean over declared seeds | 10000 |  | 41.6233 | mean elapsed per call |
| P2 | monte_carlo explicit window, mean over declared seeds | 10000 |  | 26.7133 | mean elapsed per call |
| P1 | monte_carlo explicit window, mean over declared seeds | 50000 |  | 184.4200 | mean elapsed per call |
| P2 | monte_carlo explicit window, mean over declared seeds | 50000 |  | 107.7733 | mean elapsed per call |
