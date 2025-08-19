python qcds_full_gpt_22.py \
  --mode aer \
  --genome GRCh38 \
  --shots 1024 \
  --grover-max-iters 15 \
  --signal-mode rotate \
  --signal-seed 1337 \
  --evidence-bits 12 \
  --layers 32 \
  --iters-sched adaptive \
  --adaptive \
  --pilot-p0 \
  --pilot-p0-shots 4096 \
  --pilot-window 3 \
  --scout-shots 256 \
  --confirm-shots 4096 \
  --patience 3 \
  --drop-tol 0.03 \
  --multihit-threshold 0.90 \
  --oracle-type weighted \
  --oracle-weight-angle 2.0 \
  --oracle-power 3 \
  --print-hits \
  --layman \
  --out-dir results




python qcds_full_gpt_22.py \
  --mode ibm \
  --ibm-backend ibm_brisbane \
  --genome GRCh38 \
  --shots 1024 \
  --grover-max-iters 15 \
  --signal-mode rotate \
  --signal-seed 1337 \
  --evidence-bits 12 \
  --layers 16 \
  --iters-sched adaptive \
  --adaptive \
  --pilot-p0 \
  --pilot-p0-shots 2048 \
  --pilot-window 2 \
  --scout-shots 128 \
  --confirm-shots 2048 \
  --patience 3 \
  --drop-tol 0.03 \
  --multihit-threshold 0.90 \
  --oracle-type weighted \
  --oracle-weight-angle 1.4 \
  --oracle-power 2 \
  --print-hits \
  --layman \
  --out-dir results


python qcds_full_gpt_22.py \
  --mode aer \
  --genome GRCh38 \
  --shots 512 \
  --grover-max-iters 15 \
  --signal-mode rotate \
  --signal-seed 1337 \
  --evidence-bits 12 \
  --layers 16 \
  --iters-sched adaptive \
  --adaptive \
  --pilot-p0 \
  --pilot-p0-shots 2048 \
  --pilot-window 3 \
  --scout-shots 128 \
  --confirm-shots 2048 \
  --patience 2 \
  --drop-tol 0.05 \
  --multihit-threshold 0.90 \
  --oracle-type weighted \
  --oracle-weight-angle 1.6 \
  --oracle-power 2 \
  --layman \
  --out-dir results
