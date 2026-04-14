# README

This is the artifact of the paper "Practical Silent Threshold Signatures and Silent Threshold Encryption for Dynamic Committees", submitted to CCS 2026 cycle B.

## Evaluation

Performance of epoch representation computation: `cargo bench --bench epoch`, corresponding to Table 1.

Performance of c-STS and c-STE schemes: `cargo run --example endtoend --release -- sign/encrypt <N> <n> <t>`, the program will output the time of each step.

- e.g., `cargo run --example endtoend --release -- encrypt 8192 4096 2048` will demonstrate the performance of c-STE under $N=8192$, $n=4096$ and $t=2048$.
- $N$ should be a power of 2.
- The time of "Aggregation" (when the first parameter is `sign`) corresponds to Table 2.

Performance of hints computation: `cargo bench --bench hints`, corresponding to Table 4.
