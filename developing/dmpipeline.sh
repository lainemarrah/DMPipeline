#!/bin/bash

conda activate dmpipeline

cd ${outdir}
if [-s ${name}.sorted]; then
  python scripts/pipeline.py [ARGS]
else
  python scripts/align_and_sort.py [ARGS]
  python scripts/pipeline.py
fi

if [-s ${name}.dmrpt]; then
  rm -r ${tempdir}
fi

conda deactivate
